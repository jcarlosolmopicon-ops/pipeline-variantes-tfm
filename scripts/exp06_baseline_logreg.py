#!/usr/bin/env python3
"""
exp06_baseline_logreg.py

Línea base: regresión logística sobre el mismo conjunto de características,
el mismo split y la misma validación cruzada que exp04.

Objetivo: comprobar si el AUC-ROC de 0,96 obtenido por Random Forest y XGBoost
refleja capacidad real del modelo o si un clasificador lineal simple lo alcanza
igualmente. Es la comparación que faltaba para contextualizar los resultados.

Salida: resultados/exp06/baseline_results.json + una línea lista para pegar
en la Tabla 6 del TFM.

Uso (desde la raíz del repositorio):
    python3 scripts/exp06_baseline_logreg.py

Tiempo estimado: 5-15 min (la regresión logística sobre 3,6 M filas con
saga/liblinear puede tardar; se usa lbfgs con max_iter alto y escalado previo).
"""

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import (
    GroupShuffleSplit,
    StratifiedGroupKFold,
    cross_validate,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

# ─────────────────────────────────────────────────────────────────────────
MAF_PATH = Path("datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz")
OUTDIR = Path("resultados/exp06")
OUTDIR.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42

RELEVANT = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation",
    "Translation_Start_Site", "In_Frame_Del", "In_Frame_Ins",
}
NONRELEVANT = {"Silent", "Intron", "3'UTR", "5'UTR", "RNA", "3'Flank", "5'Flank"}

COLS = [
    "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type",
    "SIFT", "PolyPhen", "t_depth", "t_alt_count", "n_depth",
    "COSMIC", "Tumor_Seq_Allele2", "Reference_Allele",
]

# ─────────────────────────────────────────────────────────────────────────
print(f"[{time.strftime('%H:%M:%S')}] Cargando MAF...")
df = pd.read_csv(MAF_PATH, sep="\t", comment=None, skiprows=1,
                 usecols=COLS, low_memory=False, na_values=[".", ""])
print(f"  {len(df):,} variantes")

df = df[df["Variant_Classification"].isin(RELEVANT | NONRELEVANT)].copy()
df["label"] = df["Variant_Classification"].isin(RELEVANT).astype(int)


def parse_score(val):
    m = re.search(r"\(([\d.]+)\)", str(val))
    return float(m.group(1)) if m else np.nan


df["SIFT_score"] = df["SIFT"].apply(parse_score)
df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
df["SIFT_missing"] = df["SIFT_score"].isna().astype(int)
df["PolyPhen_missing"] = df["PolyPhen_score"].isna().astype(int)
df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)
df["COSMIC_present"] = df["COSMIC"].notna().astype(int)
df["indel_len"] = (df["Tumor_Seq_Allele2"].astype(str).str.len()
                   - df["Reference_Allele"].astype(str).str.len()).abs()
df = df.dropna(subset=["VAF", "t_depth"])
print(f"  {len(df):,} variantes tras dropna")

# MISMAS características que exp04 (incluidos los flags de ausencia)
numeric_features = ["SIFT_score", "PolyPhen_score", "SIFT_missing", "PolyPhen_missing",
                    "VAF", "t_depth", "n_depth", "COSMIC_present", "indel_len"]
categorical_features = ["Variant_Type"]

groups = df["Tumor_Sample_Barcode"]
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_STATE)
train_idx, test_idx = next(gss.split(df, df["label"], groups))

X_train, X_test = df.iloc[train_idx], df.iloc[test_idx]
y_train, y_test = X_train["label"], X_test["label"]
g_train = groups.iloc[train_idx]
print(f"Train: {len(X_train):,} | Test: {len(X_test):,}\n")

# StandardScaler es imprescindible para que la regresión logística converja
preprocessor = ColumnTransformer(transformers=[
    ("num", Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("sc", StandardScaler())]), numeric_features),
    ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
])

clf = LogisticRegression(max_iter=1000, class_weight="balanced",
                         random_state=RANDOM_STATE, n_jobs=-1)
pipe = Pipeline([("preprocessor", preprocessor), ("classifier", clf)])

print(f"[{time.strftime('%H:%M:%S')}] Validación cruzada (5-fold agrupada)...")
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
cvs = cross_validate(pipe, X_train, y_train, groups=g_train, cv=cv,
                     scoring=["roc_auc", "precision", "recall", "f1"], n_jobs=1)
cv_res = {k: float(np.mean(cvs[f"test_{k}"]))
          for k in ["roc_auc", "precision", "recall", "f1"]}
print(f"  CV:   AUC={cv_res['roc_auc']:.4f}  Prec={cv_res['precision']:.4f}  "
      f"Rec={cv_res['recall']:.4f}  F1={cv_res['f1']:.4f}")

print(f"[{time.strftime('%H:%M:%S')}] Entrenando modelo final y evaluando en test...")
pipe.fit(X_train, y_train)
y_pred = pipe.predict(X_test)
y_proba = pipe.predict_proba(X_test)[:, 1]
test_res = {
    "roc_auc": float(roc_auc_score(y_test, y_proba)),
    "precision": float(precision_score(y_test, y_pred)),
    "recall": float(recall_score(y_test, y_pred)),
    "f1": float(f1_score(y_test, y_pred)),
}
print(f"  Test: AUC={test_res['roc_auc']:.4f}  Prec={test_res['precision']:.4f}  "
      f"Rec={test_res['recall']:.4f}  F1={test_res['f1']:.4f}")

results = {
    "model": "LogisticRegression (baseline)",
    "features": numeric_features + categorical_features,
    "n_train_rows": len(X_train), "n_test_rows": len(X_test),
    "cv": cv_res, "test": test_res,
}
out = OUTDIR / "baseline_results.json"
json.dump(results, open(out, "w"), indent=2, ensure_ascii=False)

# ─────────────────────────────────────────────────────────────────────────

print("=" * 68)
print(f"Regresión logística | CV (train, 5-fold) | {cv_res['roc_auc']:.4f}* | "
      f"{cv_res['precision']:.4f}* | {cv_res['recall']:.4f}* | {cv_res['f1']:.4f}*")
print(f"Regresión logística | Test (held-out)    | {test_res['roc_auc']:.4f} | "
      f"{test_res['precision']:.4f} | {test_res['recall']:.4f} | {test_res['f1']:.4f}")
print("=" * 68)

delta_rf = 0.9595 - test_res["roc_auc"]
print(f"\nDiferencia de AUC frente a Random Forest (0,9595): {delta_rf:+.4f}")
if delta_rf < 0.02:
    print("→ La línea base alcanza un rendimiento muy similar. Interpretación:")
    print("  el problema está determinado por las características")
    print("  de anotación funcional, y los modelos de conjunto aportan poco margen.")
else:
    print("→ Los modelos de conjunto superan claramente a la línea base lineal,")
print(f"\nResultados guardados en {out}")
