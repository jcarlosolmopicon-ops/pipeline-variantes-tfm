#!/usr/bin/env python3
"""
exp05_ablation_missing_flags.py

Experimento de ablación — cuantifica la contribución de los indicadores
estructurales SIFT_missing y PolyPhen_missing al AUC-ROC del módulo de
clasificación (sección 5.3.2 del TFM).

Reutiliza exactamente la misma lógica de carga, preprocesamiento, split
y validación cruzada que exp04 (train_ml_module.py), eliminando
únicamente las dos columnas estructurales del conjunto de features.

Salida: resultados/exp05/ (mismo formato que resultados/exp04/ para
permitir comparación directa).

Uso:
    python3 exp05_ablation_missing_flags.py

Requisitos: el mismo entorno conda/venv usado para exp04
(scikit-learn, xgboost, pandas, numpy).
"""

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import (
    GroupShuffleSplit,
    StratifiedGroupKFold,
    cross_validate,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBClassifier

# ─────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN — ajustar rutas si difieren en tu máquina
# ─────────────────────────────────────────────────────────────────────────
MAF_PATH = Path("datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz")
OUTDIR = Path("resultados/exp05")
OUTDIR.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
N_SPLITS_CV = 5
TEST_SIZE = 0.2

RELEVANT_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation",
    "Translation_Start_Site", "In_Frame_Del", "In_Frame_Ins",
}
NONRELEVANT_CLASSES = {
    "Silent", "Intron", "3'UTR", "5'UTR", "RNA", "3'Flank", "5'Flank",
}

COLS = [
    "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type",
    "SIFT", "PolyPhen", "t_depth", "t_alt_count", "n_depth",
    "COSMIC", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Reference_Allele",
]

# ─────────────────────────────────────────────────────────────────────────
# 1. CARGA — idéntica a train_ml_module.py (sección 4.5 del TFM)
# ─────────────────────────────────────────────────────────────────────────
print(f"[{time.strftime('%H:%M:%S')}] Cargando MAF...")
df = pd.read_csv(
    MAF_PATH, sep="\t", comment=None, skiprows=1,
    usecols=COLS, low_memory=False,
    na_values=[".", ""],
)
print(f"  {len(df):,} variantes cargadas")

df = df[df["Variant_Classification"].isin(RELEVANT_CLASSES | NONRELEVANT_CLASSES)].copy()
df["label"] = df["Variant_Classification"].isin(RELEVANT_CLASSES).astype(int)
print(f"  {len(df):,} variantes tras filtrar por clase")
print(f"  Relevantes: {df['label'].sum():,} ({df['label'].mean()*100:.1f}%)")


def parse_score(val):
    m = re.search(r"\(([\d.]+)\)", str(val))
    if m:
        return float(m.group(1))
    return np.nan


df["SIFT_score"] = df["SIFT"].apply(parse_score)
df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
df["SIFT_missing"] = df["SIFT_score"].isna().astype(int)
df["PolyPhen_missing"] = df["PolyPhen_score"].isna().astype(int)

df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)
df["COSMIC_present"] = df["COSMIC"].notna().astype(int)
df["indel_len"] = (
    df["Tumor_Seq_Allele2"].astype(str).str.len()
    - df["Reference_Allele"].astype(str).str.len()
).abs()

df = df.dropna(subset=["VAF", "t_depth"])
print(f"  {len(df):,} variantes tras dropna(VAF, t_depth)")

# ─────────────────────────────────────────────────────────────────────────
# 2. DEFINICIÓN DE FEATURES — ÚNICA DIFERENCIA RESPECTO A exp04
# ─────────────────────────────────────────────────────────────────────────
# exp04 (referencia): numeric_features incluía SIFT_missing y PolyPhen_missing
# exp05 (este experimento): se eliminan ambas columnas
numeric_features = [
    "SIFT_score", "PolyPhen_score",
    # "SIFT_missing", "PolyPhen_missing",   <-- ELIMINADAS EN exp05
    "VAF", "t_depth", "n_depth", "COSMIC_present", "indel_len",
]
categorical_features = ["Variant_Type"]

print(f"\n[ABLACIÓN] Features numéricas usadas ({len(numeric_features)}): {numeric_features}")
print("[ABLACIÓN] SIFT_missing y PolyPhen_missing EXCLUIDAS del modelo\n")

df_model = df.dropna(subset=["label"]).copy()
groups = df_model["Tumor_Sample_Barcode"]

# ─────────────────────────────────────────────────────────────────────────
# 3. SPLIT — idéntico a exp04 (mismo random_state -> mismas muestras en test)
# ─────────────────────────────────────────────────────────────────────────
gss = GroupShuffleSplit(n_splits=1, test_size=TEST_SIZE, random_state=RANDOM_STATE)
train_idx, test_idx = next(gss.split(df_model, df_model["label"], groups))

X_train = df_model.iloc[train_idx]
X_test = df_model.iloc[test_idx]
y_train = X_train["label"]
y_test = X_test["label"]
g_train = groups.iloc[train_idx]

print(f"Train: {len(X_train):,} variantes de {g_train.nunique():,} muestras")
print(f"Test:  {len(X_test):,} variantes de {groups.iloc[test_idx].nunique():,} muestras\n")

# ─────────────────────────────────────────────────────────────────────────
# 4. PREPROCESADOR Y MODELOS — mismos hiperparámetros que exp04
# ─────────────────────────────────────────────────────────────────────────
preprocessor = ColumnTransformer(transformers=[
    ("num", SimpleImputer(strategy="median"), numeric_features),
    ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
])

models = {
    "RandomForest": RandomForestClassifier(
        n_estimators=200, max_depth=12, min_samples_leaf=20,
        class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1,
    ),
    "XGBoost": XGBClassifier(
        n_estimators=300, max_depth=6, learning_rate=0.1,
        subsample=0.8, colsample_bytree=0.8,
        eval_metric="logloss", random_state=RANDOM_STATE, n_jobs=-1,
    ),
}

results = {"ablation": "SIFT_missing_PolyPhen_missing_removed", "models": {}}

for name, clf in models.items():
    print(f"[{time.strftime('%H:%M:%S')}] === {name} ===")
    pipe = Pipeline([("preprocessor", preprocessor), ("classifier", clf)])

    # --- Validación cruzada sobre train ---
    cv = StratifiedGroupKFold(n_splits=N_SPLITS_CV, shuffle=True, random_state=RANDOM_STATE)
    cv_scores = cross_validate(
        pipe, X_train, y_train, groups=g_train, cv=cv,
        scoring=["roc_auc", "precision", "recall", "f1"],
        n_jobs=1,
    )
    cv_summary = {
        "auc_roc": float(np.mean(cv_scores["test_roc_auc"])),
        "precision": float(np.mean(cv_scores["test_precision"])),
        "recall": float(np.mean(cv_scores["test_recall"])),
        "f1": float(np.mean(cv_scores["test_f1"])),
    }
    print(f"  CV (5-fold):  AUC={cv_summary['auc_roc']:.4f}  "
          f"Prec={cv_summary['precision']:.4f}  Rec={cv_summary['recall']:.4f}  "
          f"F1={cv_summary['f1']:.4f}")

    # --- Reentrenar sobre todo train, evaluar en test held-out ---
    pipe.fit(X_train, y_train)
    y_pred = pipe.predict(X_test)
    y_proba = pipe.predict_proba(X_test)[:, 1]

    test_summary = {
        "auc_roc": float(roc_auc_score(y_test, y_proba)),
        "precision": float(precision_score(y_test, y_pred)),
        "recall": float(recall_score(y_test, y_pred)),
        "f1": float(f1_score(y_test, y_pred)),
    }
    print(f"  Test held-out: AUC={test_summary['auc_roc']:.4f}  "
          f"Prec={test_summary['precision']:.4f}  Rec={test_summary['recall']:.4f}  "
          f"F1={test_summary['f1']:.4f}\n")

    results["models"][name] = {"cv": cv_summary, "test": test_summary}

    # Importancia de features (solo RF, para comparar con Tabla 4 del TFM)
    if name == "RandomForest":
        feat_names = (
            numeric_features
            + list(pipe.named_steps["preprocessor"]
                   .named_transformers_["cat"].get_feature_names_out(categorical_features))
        )
        importances = pipe.named_steps["classifier"].feature_importances_
        imp_dict = dict(sorted(zip(feat_names, importances.tolist()),
                                key=lambda x: -x[1]))
        results["models"][name]["feature_importance"] = imp_dict

# ─────────────────────────────────────────────────────────────────────────
# 5. GUARDAR RESULTADOS
# ─────────────────────────────────────────────────────────────────────────
results["n_train_rows"] = len(X_train)
results["n_test_rows"] = len(X_test)
results["n_train_groups"] = int(g_train.nunique())
results["n_test_groups"] = int(groups.iloc[test_idx].nunique())
results["numeric_features_used"] = numeric_features

out_json = OUTDIR / "ablation_results.json"
with open(out_json, "w") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print(f"[{time.strftime('%H:%M:%S')}] Resultados guardados en {out_json}")

# ─────────────────────────────────────────────────────────────────────────
# 6. COMPARATIVA DIRECTA CON exp04 (si el fichero existe)
# ─────────────────────────────────────────────────────────────────────────
exp04_results_path = Path("resultados/exp04/ml_results.json")
if exp04_results_path.exists():
    with open(exp04_results_path) as f:
        exp04 = json.load(f)

    print("\n" + "=" * 70)
    print("COMPARATIVA exp04 (con SIFT_missing/PolyPhen_missing) vs "
          "exp05 (sin ellas)")
    print("=" * 70)
    for model_name in ["RandomForest", "XGBoost"]:
        try:
            auc_exp04 = exp04[model_name]["test"]["auc_roc"] if "test" in exp04.get(model_name, {}) else exp04.get(model_name, {}).get("auc_roc")
        except Exception:
            auc_exp04 = None
        auc_exp05 = results["models"][model_name]["test"]["auc_roc"]
        print(f"{model_name:15s}  exp04 AUC: {auc_exp04}   "
              f"exp05 AUC: {auc_exp05:.4f}")
    print("\nNOTA: revisa manualmente el formato exacto de "
          "resultados/exp04/ml_results.json si la comparativa automática "
          "no imprime el valor de exp04 (puede variar la clave usada).")
else:
    print(f"\n[AVISO] No se encontró {exp04_results_path} — "
          "compara manualmente contra la Tabla 3 del TFM (AUC test: "
          "RF=0.9595, XGBoost=0.9596).")

print(f"\n[{time.strftime('%H:%M:%S')}] Experimento exp05 completado.")
