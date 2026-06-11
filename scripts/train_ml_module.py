#!/usr/bin/env python3
"""
Modulo ML - Clasificacion de impacto funcional de variantes somaticas (TCGA-LUAD)

Objetivo: clasificar variantes como funcionalmente relevantes (missense, nonsense,
frameshift, splice_site, nonstop, translation_start_site) vs. no relevantes
(silent, intron, UTR, flank, RNA) usando Random Forest y XGBoost.

Split por caso (Tumor_Sample_Barcode) para evitar data leakage entre
train/test del mismo paciente (Whalen et al., 2022).
"""

import re
import sys
import json
import gzip
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold, cross_validate
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import (roc_auc_score, precision_score, recall_score,
                              f1_score, classification_report, confusion_matrix)
import xgboost as xgb

MAF_PATH = "datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz"
OUT_DIR = "resultados/exp02_ml"

COLS = [
    "Hugo_Symbol", "Variant_Classification", "Variant_Type",
    "Tumor_Sample_Barcode", "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count",
    "SIFT", "PolyPhen", "IMPACT", "COSMIC",
    "Reference_Allele", "Tumor_Seq_Allele2",
]

# Clases consideradas funcionalmente relevantes (label = 1)
RELEVANT_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site",
    "In_Frame_Del", "In_Frame_Ins"
}
# Clases no relevantes (label = 0)
NONRELEVANT_CLASSES = {
    "Silent", "Intron", "3'UTR", "5'UTR", "RNA", "3'Flank", "5'Flank"
}


def parse_score(val):
    """Extrae el valor numerico de campos tipo 'deleterious(0.02)' o 'benign(0.99)'."""
    if pd.isna(val) or val in (".", ""):
        return np.nan
    m = re.search(r"\(([\d.]+)\)", str(val))
    if m:
        return float(m.group(1))
    return np.nan


def main():
    print(f"Leyendo MAF: {MAF_PATH}", file=sys.stderr)
    # El fichero tiene una linea de comentario inicial (#liftOver...)
    df = pd.read_csv(
        MAF_PATH, sep="\t", comment=None, skiprows=1,
        usecols=COLS, low_memory=False,
        na_values=[".", ""],
    )
    print(f"Filas leidas: {len(df):,}", file=sys.stderr)

    # --- Filtrar a clases consideradas (binario) ---
    df = df[df["Variant_Classification"].isin(RELEVANT_CLASSES | NONRELEVANT_CLASSES)].copy()
    print(f"Filas tras filtrar clases: {len(df):,}", file=sys.stderr)

    df["label"] = df["Variant_Classification"].isin(RELEVANT_CLASSES).astype(int)

    # --- Features derivadas ---
    df["SIFT_score"] = df["SIFT"].apply(parse_score)
    df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
    df["SIFT_missing"] = df["SIFT_score"].isna().astype(int)
    df["PolyPhen_missing"] = df["PolyPhen_score"].isna().astype(int)

    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df["n_depth"] = pd.to_numeric(df["n_depth"], errors="coerce")
    df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)

    df["COSMIC_present"] = (~df["COSMIC"].isin(["NONE", np.nan]) & df["COSMIC"].notna()).astype(int)

    df["ref_len"] = df["Reference_Allele"].astype(str).str.replace("-", "", regex=False).str.len()
    df["alt_len"] = df["Tumor_Seq_Allele2"].astype(str).str.replace("-", "", regex=False).str.len()
    df["indel_len"] = (df["alt_len"] - df["ref_len"]).abs()

    df["Variant_Type"] = df["Variant_Type"].fillna("UNKNOWN")

    numeric_features = [
        "SIFT_score", "PolyPhen_score", "SIFT_missing", "PolyPhen_missing",
        "VAF", "t_depth", "n_depth", "COSMIC_present", "indel_len",
    ]
    categorical_features = ["Variant_Type"]

    feature_cols = numeric_features + categorical_features
    df_model = df[feature_cols + ["label", "Tumor_Sample_Barcode"]].copy()

    # Eliminar filas con VAF/t_depth nulos (no calculables)
    before = len(df_model)
    df_model = df_model.dropna(subset=["VAF", "t_depth"])
    print(f"Filas tras dropna VAF/t_depth: {len(df_model):,} (eliminadas {before - len(df_model):,})",
          file=sys.stderr)

    print(f"Distribucion de label:\n{df_model['label'].value_counts()}", file=sys.stderr)
    print(f"Numero de casos (Tumor_Sample_Barcode): {df_model['Tumor_Sample_Barcode'].nunique()}",
          file=sys.stderr)

    # --- Split train/test por caso (80/20) ---
    groups = df_model["Tumor_Sample_Barcode"]
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    train_idx, test_idx = next(gss.split(df_model, df_model["label"], groups))

    X = df_model[feature_cols]
    y = df_model["label"]
    g = df_model["Tumor_Sample_Barcode"]

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    g_train = g.iloc[train_idx]

    n_train_cases = g_train.nunique()
    n_test_cases = g.iloc[test_idx].nunique()
    print(f"Train: {len(X_train):,} variantes, {n_train_cases} casos", file=sys.stderr)
    print(f"Test:  {len(X_test):,} variantes, {n_test_cases} casos", file=sys.stderr)

    # --- Preprocesador ---
    preprocessor = ColumnTransformer(transformers=[
        ("num", SimpleImputer(strategy="median"), numeric_features),
        ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
    ])

    # --- Random Forest ---
    rf_pipe = Pipeline([
        ("prep", preprocessor),
        ("clf", RandomForestClassifier(
            n_estimators=200, max_depth=12, min_samples_leaf=20,
            n_jobs=-1, random_state=42, class_weight="balanced"
        )),
    ])

    # --- XGBoost ---
    xgb_pipe = Pipeline([
        ("prep", preprocessor),
        ("clf", xgb.XGBClassifier(
            n_estimators=300, max_depth=6, learning_rate=0.1,
            subsample=0.8, colsample_bytree=0.8,
            eval_metric="logloss", random_state=42, n_jobs=-1,
        )),
    ])

    results = {}

    # --- Cross-validation estratificada por grupo (5-fold) sobre train ---
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    for name, pipe in [("RandomForest", rf_pipe), ("XGBoost", xgb_pipe)]:
        print(f"\n=== {name}: validacion cruzada (5-fold, stratified by group) ===", file=sys.stderr)
        cv_scores = cross_validate(
            pipe, X_train, y_train, groups=g_train, cv=cv,
            scoring=["roc_auc", "precision", "recall", "f1"],
            n_jobs=1,
        )
        cv_summary = {
            metric: {"mean": float(np.mean(cv_scores[f"test_{metric}"])),
                     "std": float(np.std(cv_scores[f"test_{metric}"]))}
            for metric in ["roc_auc", "precision", "recall", "f1"]
        }
        print(json.dumps(cv_summary, indent=2), file=sys.stderr)

        # --- Entrenar en todo train, evaluar en test held-out ---
        pipe.fit(X_train, y_train)
        y_pred = pipe.predict(X_test)
        y_proba = pipe.predict_proba(X_test)[:, 1]

        test_metrics = {
            "roc_auc": float(roc_auc_score(y_test, y_proba)),
            "precision": float(precision_score(y_test, y_pred)),
            "recall": float(recall_score(y_test, y_pred)),
            "f1": float(f1_score(y_test, y_pred)),
        }
        cm = confusion_matrix(y_test, y_pred).tolist()
        report = classification_report(y_test, y_pred, output_dict=True)

        print(f"\n{name} - Test set (held-out cases):", file=sys.stderr)
        print(json.dumps(test_metrics, indent=2), file=sys.stderr)

        results[name] = {
            "cv_train": cv_summary,
            "test": test_metrics,
            "confusion_matrix": cm,
            "classification_report": report,
        }

        # Feature importance (RF) / gain (XGB)
        if name == "RandomForest":
            ohe = pipe.named_steps["prep"].named_transformers_["cat"]
            cat_names = list(ohe.get_feature_names_out(categorical_features))
            all_names = numeric_features + cat_names
            importances = pipe.named_steps["clf"].feature_importances_
            results[name]["feature_importance"] = dict(
                sorted(zip(all_names, importances.tolist()), key=lambda x: -x[1])
            )

    # --- Guardar resultados ---
    import os
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(f"{OUT_DIR}/ml_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Resumen tabular
    summary_rows = []
    for name in results:
        t = results[name]["test"]
        summary_rows.append({
            "model": name,
            "test_roc_auc": round(t["roc_auc"], 4),
            "test_precision": round(t["precision"], 4),
            "test_recall": round(t["recall"], 4),
            "test_f1": round(t["f1"], 4),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(f"{OUT_DIR}/ml_summary.csv", index=False)
    print("\n=== RESUMEN FINAL ===", file=sys.stderr)
    print(summary_df.to_string(index=False), file=sys.stderr)

    # Metadata del dataset
    meta = {
        "n_total_variants": int(len(df)),
        "n_model_variants": int(len(df_model)),
        "n_train": int(len(X_train)),
        "n_test": int(len(X_test)),
        "n_train_cases": int(n_train_cases),
        "n_test_cases": int(n_test_cases),
        "label_distribution": df_model["label"].value_counts().to_dict(),
        "features": feature_cols,
    }
    with open(f"{OUT_DIR}/dataset_metadata.json", "w") as f:
        json.dump(meta, f, indent=2, default=str)

    print("\nResultados guardados en", OUT_DIR, file=sys.stderr)


if __name__ == "__main__":
    main()
