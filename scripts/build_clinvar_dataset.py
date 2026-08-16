#!/usr/bin/env python3
"""
exp07 - Construccion del conjunto de entrenamiento con etiqueta externa.

Cruza el MAF pan-cancer MC3 con ClinVar (GRCh38) por coincidencia exacta de
cromosoma, posicion y alelos, y retiene las variantes missense clasificadas
como patogenicas o benignas con criterios de revision declarados.

Salida: resultados/exp07/clinvar_missense_dataset.tsv
"""

import gzip
import json
import re

import numpy as np
import pandas as pd

CLINVAR = "datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz"
MAF = "datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz"
OUT_DIR = "resultados/exp07"

PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}

MAF_COLS = [
    "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
    "Variant_Type", "Variant_Classification", "Hugo_Symbol",
    "Tumor_Sample_Barcode", "SIFT", "PolyPhen",
    "t_depth", "t_alt_count", "n_depth", "COSMIC",
]


def info_field(info, key):
    m = re.search(rf"(?:^|;){key}=([^;]*)", info)
    return m.group(1) if m else ""


def parse_score(val):
    """Extrae el valor numerico de campos tipo 'deleterious(0.02)'."""
    m = re.search(r"\(([\d.]+)\)", str(val))
    return float(m.group(1)) if m else np.nan


def load_clinvar():
    """Devuelve {chrom:pos:ref:alt -> (CLNSIG, ONC, CLNREVSTAT)}."""
    cv = {}
    with gzip.open(CLINVAR, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            chrom, pos, ref, alt, info = f[0], f[1], f[3], f[4], f[7]
            if alt == "." or "," in alt:
                continue
            sig = info_field(info, "CLNSIG")
            if not sig:
                continue
            key = f"{chrom}:{pos}:{ref}:{alt}"
            cv[key] = (sig, info_field(info, "ONC"), info_field(info, "CLNREVSTAT"))
    return cv


def main():
    print("Leyendo ClinVar...")
    cv = load_clinvar()
    print(f"  {len(cv):,} registros con clasificacion")

    print("Cruzando con el MAF (SNP missense)...")
    hits = []
    n_missense = 0
    for chunk in pd.read_csv(MAF, sep="\t", skiprows=1, usecols=MAF_COLS,
                             low_memory=False, chunksize=500_000, dtype=str):
        sub = chunk[(chunk["Variant_Type"] == "SNP") &
                    (chunk["Variant_Classification"] == "Missense_Mutation")]
        n_missense += len(sub)
        keys = (sub["Chromosome"] + ":" + sub["Start_Position"] + ":" +
                sub["Reference_Allele"] + ":" + sub["Tumor_Seq_Allele2"])
        matched = keys.map(cv.get)
        found = matched.notna()
        if found.any():
            s = sub[found].copy()
            s["variant_key"] = keys[found]
            s["CLNSIG"] = [x[0] for x in matched[found]]
            s["ONC"] = [x[1] for x in matched[found]]
            s["CLNREVSTAT"] = [x[2] for x in matched[found]]
            hits.append(s)

    df = pd.concat(hits, ignore_index=True)
    print(f"  missense en el MAF: {n_missense:,}")
    print(f"  con registro en ClinVar: {len(df):,} filas "
          f"({df['variant_key'].nunique():,} variantes)")

    # --- Etiqueta binaria y filtro por calidad de la curacion ---
    df = df[df["CLNSIG"].isin(PATHOGENIC | BENIGN)].copy()
    df["label"] = df["CLNSIG"].isin(PATHOGENIC).astype(int)
    n_pre = df["variant_key"].nunique()
    df = df[df["CLNREVSTAT"].str.contains("criteria_provided", na=False)].copy()
    print(f"  etiquetables: {n_pre:,} -> con criterios declarados: "
          f"{df['variant_key'].nunique():,}")

    # --- Features ---
    df["SIFT_score"] = df["SIFT"].apply(parse_score)
    df["PolyPhen_score"] = df["PolyPhen"].apply(parse_score)
    for col in ("t_depth", "t_alt_count", "n_depth"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["VAF"] = (df["t_alt_count"] / df["t_depth"]).replace([np.inf, -np.inf], np.nan)
    df["COSMIC_present"] = (~df["COSMIC"].isin(["NONE"]) & df["COSMIC"].notna()).astype(int)

    # Una variante aparece en varias muestras y la etiqueta es de la variante:
    # se colapsa a una fila y la recurrencia pasa a ser una feature.
    print("Agregando a nivel de variante...")
    data = df.groupby("variant_key").agg(
        Chromosome=("Chromosome", "first"),
        Start_Position=("Start_Position", "first"),
        Reference_Allele=("Reference_Allele", "first"),
        Tumor_Seq_Allele2=("Tumor_Seq_Allele2", "first"),
        Hugo_Symbol=("Hugo_Symbol", "first"),
        CLNSIG=("CLNSIG", "first"),
        ONC=("ONC", "first"),
        CLNREVSTAT=("CLNREVSTAT", "first"),
        label=("label", "first"),
        SIFT_score=("SIFT_score", "first"),
        PolyPhen_score=("PolyPhen_score", "first"),
        VAF_median=("VAF", "median"),
        t_depth_median=("t_depth", "median"),
        n_depth_median=("n_depth", "median"),
        COSMIC_present=("COSMIC_present", "max"),
        n_samples_MC3=("Tumor_Sample_Barcode", "nunique"),
    ).reset_index()

    # Sin imputacion: se exigen ambos scores presentes.
    n_pre = len(data)
    data = data.dropna(subset=["SIFT_score", "PolyPhen_score"])
    print(f"  descartadas por falta de SIFT o PolyPhen: {n_pre - len(data):,}")

    data.to_csv(f"{OUT_DIR}/clinvar_missense_dataset.tsv", sep="\t", index=False)

    meta = {
        "clinvar_release": "GRCh38, fileDate 2026-08-08",
        "n_variants": int(len(data)),
        "n_pathogenic": int(data["label"].sum()),
        "n_benign": int((data["label"] == 0).sum()),
        "prevalence": round(float(data["label"].mean()), 4),
        "n_genes": int(data["Hugo_Symbol"].nunique()),
        "n_dropped_missing_scores": int(n_pre - len(data)),
        "n_with_oncogenicity": int((data["ONC"] != "").sum()),
        "clnsig_counts": data["CLNSIG"].value_counts().to_dict(),
    }
    with open(f"{OUT_DIR}/dataset_metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"\n{meta['n_variants']:,} variantes | {meta['n_pathogenic']:,} patogenicas "
          f"({meta['prevalence']*100:.1f} %) | {meta['n_genes']:,} genes")
    print(f"Guardado en {OUT_DIR}/clinvar_missense_dataset.tsv")


if __name__ == "__main__":
    main()
