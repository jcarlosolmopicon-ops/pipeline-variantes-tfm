#!/usr/bin/env python3
"""
Entrenamiento del modelo final (RF y XGBoost) sobre TODO TCGA-LUAD,
sin holdout, usando solo features extraibles tanto de TCGA-LUAD (MAF)
como del VCF anotado de HCC1395 (sin COSMIC, no disponible en VEP
offline para HCC1395 - ver seccion 4 del TFM).

Los modelos entrenados se serializan para su aplicacion posterior
sobre resultados/exp03/vep/annotated.vcf.gz
(scripts/apply_model_hcc1395.py).
"""

import re
import sys
import json
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
import xgboost as xgb
import joblib

MAF_PATH = "datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz"
OUT_DIR = "resultados/exp02_ml"

COLS = [
    "Hugo_Symbol", "Variant_Classification", "Variant_Type",
    "Tumor_Sample_Barcode", "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count",
    "SIFT", "PolyPhen",
    "Reference_Allele", "Tumor_Seq_Allele2",
]

RELEVANT_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site",
    "In_Frame_Del", "In_Frame_Ins"
}
NONRELEVANT_CLASSES = {
    "Silent", "Intron", "3'UTR", "5'UTR", "RNA", "3'Flank", "5'Flank"
}

# Features comunes a TCGA-LUAD (MAF) y HCC1395 (VCF anotado VEP), sin COSMIC
NUMERIC_FEATURES = [
    "SIFT_score", "PolyPhen_score", "SIFT_missing", "PolyPhen_missing",
    "VAF", "t_depth", "n_depth", "indel_len",
]
CATEGORICAL_FEATURES = ["Variant_Type"]
FEATURE_COLS = NUMERIC_FEATURES + CATEGORICAL_FEATURES


def parse_score(val):
    if pd.isna(val) or val in (".", ""):
        return np.nan
    m = re.search(r"\(([\d.]+)\)", str(val))
    if m:
        return float(m.group(1))
    return np.nan


def main():
    print(f"Leyendo MAF: {MAF_PATH}", file=sys.stderr)
    df = pd.read_csv(
        MAF_PATH, sep="\t", skiprows=1,
        usecols=COLS, low_memory=False,
        na_values=[".", ""],
    )
    print(f"Filas leidas: {len(df):,}", file=sys.stderr)

    df = df[df["Variant_Classification"].isin(RELEVANT_CLASSES | NONRELEVANT_CLASSES)].copy()
    df["label"] = df["Variant_Classification"].isin(RELEVANT_CLASSES).astype(int)

    df["SIFT_score"] = df["SIFT"].apply(parse_score)
    df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
    df["SIFT_missing"] = df["SIFT_score"].isna().astype(int)
    df["PolyPhen_missing"] = df["PolyPhen_score"].isna().astype(int)

    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df["n_depth"] = pd.to_numeric(df["n_depth"], errors="coerce")
    df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)

    df["ref_len"] = df["Reference_Allele"].astype(str).str.replace("-", "", regex=False).str.len()
    df["alt_len"] = df["Tumor_Seq_Allele2"].astype(str).str.replace("-", "", regex=False).str.len()
    df["indel_len"] = (df["alt_len"] - df["ref_len"]).abs()

    df["Variant_Type"] = df["Variant_Type"].fillna("UNKNOWN")

    df_model = df[FEATURE_COLS + ["label"]].copy()
    before = len(df_model)
    df_model = df_model.dropna(subset=["VAF", "t_depth"])
    print(f"Filas tras dropna: {len(df_model):,} (eliminadas {before - len(df_model):,})",
          file=sys.stderr)

    X = df_model[FEATURE_COLS]
    y = df_model["label"]

    preprocessor = ColumnTransformer(transformers=[
        ("num", SimpleImputer(strategy="median"), NUMERIC_FEATURES),
        ("cat", OneHotEncoder(handle_unknown="ignore"), CATEGORICAL_FEATURES),
    ])

    rf_pipe = Pipeline([
        ("prep", preprocessor),
        ("clf", RandomForestClassifier(
            n_estimators=200, max_depth=12, min_samples_leaf=20,
            n_jobs=-1, random_state=42, class_weight="balanced"
        )),
    ])

    xgb_pipe = Pipeline([
        ("prep", preprocessor),
        ("clf", xgb.XGBClassifier(
            n_estimators=300, max_depth=6, learning_rate=0.1,
            subsample=0.8, colsample_bytree=0.8,
            eval_metric="logloss", random_state=42, n_jobs=-1,
        )),
    ])

    os.makedirs(OUT_DIR, exist_ok=True)

    for name, pipe in [("RandomForest", rf_pipe), ("XGBoost", xgb_pipe)]:
        print(f"\nEntrenando {name} sobre todo el dataset ({len(X):,} filas)...", file=sys.stderr)
        pipe.fit(X, y)
        path = f"{OUT_DIR}/model_{name}.joblib"
        joblib.dump(pipe, path)
        print(f"Guardado en {path}", file=sys.stderr)

    with open(f"{OUT_DIR}/final_model_features.json", "w") as f:
        json.dump({
            "numeric_features": NUMERIC_FEATURES,
            "categorical_features": CATEGORICAL_FEATURES,
            "feature_cols": FEATURE_COLS,
            "n_training_rows": int(len(X)),
            "label_distribution": y.value_counts().to_dict(),
            "note": "Modelo final entrenado sobre todo TCGA-LUAD, sin COSMIC "
                    "(no disponible en VEP offline para HCC1395). "
                    "Validacion de generalizacion ya realizada en train_ml_module.py "
                    "(holdout por caso, AUC~0.96).",
        }, f, indent=2, default=str)

    print("\nListo.", file=sys.stderr)


if __name__ == "__main__":
    main()
