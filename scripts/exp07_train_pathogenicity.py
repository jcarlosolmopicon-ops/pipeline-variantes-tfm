#!/usr/bin/env python3
"""
exp07 - Clasificacion de patogenicidad de variantes missense.

Entrena regresion logistica, Random Forest y XGBoost sobre las variantes
missense del MAF MC3 etiquetadas con ClinVar. La particion train/test se
agrupa por gen para impedir que el modelo memorice genes completos.

Se incluyen tres lineas base sin entrenamiento (SIFT, PolyPhen y su
combinacion) para poder interpretar la ganancia real de los modelos.

Salida: resultados/exp07/exp07_results.json
"""

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, confusion_matrix, f1_score,
                             precision_score, recall_score, roc_auc_score)
from sklearn.model_selection import (GroupShuffleSplit, StratifiedGroupKFold,
                                     cross_validate)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

DATA = "resultados/exp07/clinvar_missense_dataset.tsv"
OUT_DIR = "resultados/exp07"
SEED = 42

FEATURES = ["SIFT_score", "PolyPhen_score", "VAF_median", "t_depth_median",
            "n_depth_median", "COSMIC_present", "n_samples_MC3"]


def evaluate(name, y_true, score, pred, store):
    """Calcula las metricas de un clasificador y las imprime."""
    metrics = {
        "roc_auc": float(roc_auc_score(y_true, score)),
        "pr_auc": float(average_precision_score(y_true, score)),
        "precision": float(precision_score(y_true, pred, zero_division=0)),
        "recall": float(recall_score(y_true, pred, zero_division=0)),
        "f1": float(f1_score(y_true, pred, zero_division=0)),
        "confusion_matrix": confusion_matrix(y_true, pred).tolist(),
    }
    store[name] = metrics
    print(f"  {name:<26} AUC {metrics['roc_auc']:.4f}  PR-AUC {metrics['pr_auc']:.4f}  "
          f"F1 {metrics['f1']:.4f}  P {metrics['precision']:.4f}  R {metrics['recall']:.4f}")
    return metrics


def main():
    df = pd.read_csv(DATA, sep="\t")
    X, y, groups = df[FEATURES], df["label"].values, df["Hugo_Symbol"].values
    print(f"{len(df):,} variantes | {y.sum():,} patogenicas ({y.mean()*100:.1f} %) | "
          f"{df['Hugo_Symbol'].nunique():,} genes")

    # --- Particion agrupada por gen (80/20) ---
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    train_idx, test_idx = next(gss.split(X, y, groups))
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    g_train, g_test = groups[train_idx], groups[test_idx]
    assert not set(g_train) & set(g_test), "genes compartidos entre train y test"
    print(f"train: {len(train_idx):,} variantes / {pd.Series(g_train).nunique():,} genes | "
          f"test: {len(test_idx):,} / {pd.Series(g_test).nunique():,} genes")

    results = {
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "n_genes_train": int(pd.Series(g_train).nunique()),
        "n_genes_test": int(pd.Series(g_test).nunique()),
        "prevalence_train": round(float(y_train.mean()), 4),
        "prevalence_test": round(float(y_test.mean()), 4),
        "features": FEATURES,
        "baselines": {},
        "models": {},
    }

    # --- Lineas base sin entrenamiento ---
    # SIFT es inverso: valores bajos indican mayor impacto, de ahi el signo.
    print("\n=== Lineas base ===")
    sift, poly = X_test["SIFT_score"].values, X_test["PolyPhen_score"].values
    evaluate("SIFT", y_test, -sift, (sift < 0.05).astype(int), results["baselines"])
    evaluate("PolyPhen", y_test, poly, (poly > 0.85).astype(int), results["baselines"])
    evaluate("SIFT<0.05 y PolyPhen>0.85", y_test, (poly - sift + 1) / 2,
             ((sift < 0.05) & (poly > 0.85)).astype(int), results["baselines"])

    # --- Modelos ---
    models = {
        "LogisticRegression": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                       random_state=SEED)),
        ]),
        "RandomForest": RandomForestClassifier(
            n_estimators=400, max_depth=8, min_samples_leaf=10,
            class_weight="balanced", n_jobs=-1, random_state=SEED),
        "XGBoost": XGBClassifier(
            n_estimators=400, max_depth=4, learning_rate=0.05, subsample=0.8,
            colsample_bytree=0.8, eval_metric="logloss", random_state=SEED, n_jobs=-1),
    }

    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
    print("\n=== Modelos (test held-out; validacion cruzada agrupada por gen) ===")
    for name, model in models.items():
        scores = cross_validate(model, X_train, y_train, groups=g_train, cv=cv,
                                scoring=["roc_auc", "average_precision", "f1"], n_jobs=1)
        model.fit(X_train, y_train)
        proba = model.predict_proba(X_test)[:, 1]
        metrics = evaluate(name, y_test, proba, (proba > 0.5).astype(int), results["models"])
        metrics["cv"] = {
            m: {"mean": float(np.mean(scores[f"test_{m}"])),
                "std": float(np.std(scores[f"test_{m}"]))}
            for m in ("roc_auc", "average_precision", "f1")
        }
        if name == "RandomForest":
            metrics["feature_importance"] = dict(sorted(
                zip(FEATURES, model.feature_importances_.tolist()), key=lambda x: -x[1]))

    # --- Predicciones sobre el test, para las figuras ---
    # Las curvas ROC y de precision-sensibilidad se dibujan a partir de este
    # fichero, de modo que el script de figuras no reentrena nada.
    predicciones = pd.DataFrame({"label": y_test})
    predicciones["SIFT"] = -X_test["SIFT_score"].values
    predicciones["PolyPhen"] = X_test["PolyPhen_score"].values
    for name, model in models.items():
        predicciones[name] = model.predict_proba(X_test)[:, 1]
    predicciones.to_csv(f"{OUT_DIR}/test_predictions.csv", index=False)

    # --- Validacion externa sobre las variantes con clasificacion de oncogenicidad ---
    onc = df[df["ONC"].isin(["Oncogenic", "Likely_oncogenic"])]
    onc_test = onc[onc["Hugo_Symbol"].isin(set(g_test))]
    if len(onc_test):
        rf = models["RandomForest"]
        proba_onc = rf.predict_proba(onc_test[FEATURES])[:, 1]
        results["oncogenicity_check"] = {
            "n_total_oncogenic": int(len(onc)),
            "n_in_test_genes": int(len(onc_test)),
            "mean_proba": round(float(proba_onc.mean()), 4),
            "recall_at_0.5": round(float((proba_onc > 0.5).mean()), 4),
        }
        print(f"\nVariantes con ONC oncogenica en genes de test: {len(onc_test)} | "
              f"probabilidad media {proba_onc.mean():.3f} | "
              f"recall {(proba_onc > 0.5).mean():.3f}")

    with open(f"{OUT_DIR}/exp07_results.json", "w") as fh:
        json.dump(results, fh, indent=2)

    best_model = max(results["models"].items(), key=lambda kv: kv[1]["roc_auc"])
    best_base = max(results["baselines"].items(), key=lambda kv: kv[1]["roc_auc"])
    print(f"\nMejor modelo: {best_model[0]} AUC {best_model[1]['roc_auc']:.4f}")
    print(f"Mejor linea base: {best_base[0]} AUC {best_base[1]['roc_auc']:.4f}")
    print(f"Ganancia: {best_model[1]['roc_auc'] - best_base[1]['roc_auc']:+.4f}")
    print(f"\nGuardado en {OUT_DIR}/exp07_results.json")


if __name__ == "__main__":
    main()
