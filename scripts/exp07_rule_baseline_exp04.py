#!/usr/bin/env python3
"""
exp07 - Linea base trivial del planteamiento anterior (exp04).

Evalua una regla determinista sin entrenamiento sobre el mismo conjunto de
test que exp04, para poder comparar ambos planteamientos en igualdad de
condiciones. La regla marca como relevante toda variante que tenga
puntuacion SIFT o PolyPhen, o que sea una insercion o delecion.

Salida: resultados/exp07/exp04_rule_baseline.json
"""

import json
import re

import numpy as np
import pandas as pd
from sklearn.metrics import (confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score)
from sklearn.model_selection import GroupShuffleSplit

MAF = "datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz"
OUT_DIR = "resultados/exp07"
SEED = 42

RELEVANT_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site",
    "In_Frame_Del", "In_Frame_Ins",
}
NONRELEVANT_CLASSES = {
    "Silent", "Intron", "3'UTR", "5'UTR", "RNA", "3'Flank", "5'Flank",
}

COLS = ["Variant_Classification", "Variant_Type", "SIFT", "PolyPhen",
        "t_depth", "t_alt_count", "Tumor_Sample_Barcode"]


def parse_score(val):
    m = re.search(r"\(([\d.]+)\)", str(val))
    return float(m.group(1)) if m else np.nan


def main():
    df = pd.read_csv(MAF, sep="\t", skiprows=1, usecols=COLS,
                     low_memory=False, na_values=[".", ""])
    df = df[df["Variant_Classification"].isin(RELEVANT_CLASSES | NONRELEVANT_CLASSES)].copy()
    df["label"] = df["Variant_Classification"].isin(RELEVANT_CLASSES).astype(int)
    df["SIFT_score"] = df["SIFT"].apply(parse_score)
    df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["VAF", "t_depth"]).reset_index(drop=True)

    # Misma particion que exp04 para que las metricas sean comparables.
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    _, test_idx = next(gss.split(df, df["label"], df["Tumor_Sample_Barcode"]))
    test = df.iloc[test_idx]
    y = test["label"].values

    tiene_score = test["SIFT_score"].notna() | test["PolyPhen_score"].notna()
    regla = (tiene_score | test["Variant_Type"].isin(["DEL", "INS"])).astype(int).values

    # Cobertura de la implicacion "tiene score -> relevante" sobre todo el conjunto
    con_score = df["SIFT_score"].notna() | df["PolyPhen_score"].notna()
    precision_score_presente = float(df.loc[con_score, "label"].mean())

    resultados = {
        "descripcion": "regla determinista sin entrenamiento sobre el test de exp04",
        "n_test": int(len(test)),
        "precision": float(precision_score(y, regla)),
        "recall": float(recall_score(y, regla)),
        "f1": float(f1_score(y, regla)),
        "roc_auc_regla_binaria": float(roc_auc_score(y, regla)),
        "confusion_matrix": confusion_matrix(y, regla).tolist(),
        "implicacion_score_presente": {
            "n_variantes_con_score": int(con_score.sum()),
            "proporcion_relevantes": round(precision_score_presente, 6),
        },
    }

    with open(f"{OUT_DIR}/exp04_rule_baseline.json", "w") as fh:
        json.dump(resultados, fh, indent=2)

    print(f"test: {resultados['n_test']:,} variantes")
    print(f"regla trivial -> F1 {resultados['f1']:.4f}  "
          f"P {resultados['precision']:.4f}  R {resultados['recall']:.4f}")
    print(f"variantes con score: {con_score.sum():,}, de las cuales relevantes "
          f"{precision_score_presente*100:.4f} %")
    print(f"Guardado en {OUT_DIR}/exp04_rule_baseline.json")


if __name__ == "__main__":
    main()
