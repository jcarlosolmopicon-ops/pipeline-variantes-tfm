#!/usr/bin/env python3
"""
Aplica los modelos RF/XGBoost entrenados sobre TCGA-LUAD
(scripts/train_final_model.py) a las variantes anotadas de HCC1395
(resultados/exp03/vep/annotated.vcf.gz).

Extrae para cada variante: SIFT, PolyPhen (del bloque CSQ de VEP,
transcrito CANONICAL si existe, si no el primero), VAF/t_depth/n_depth
(de los campos FORMAT del VCF de MuTect2), Variant_Type e indel_len
(de REF/ALT). No incluye COSMIC (no disponible en VEP offline para
HCC1395 - ver seccion 4 del TFM).

Salida: TSV con columnas originales + probabilidad de relevancia
funcional segun cada modelo.
"""

import re
import sys
import gzip
import json
import joblib
import numpy as np
import pandas as pd

VCF_PATH = "resultados/exp03/vep/annotated.vcf.gz"
OUT_DIR = "resultados/exp04"
MODEL_DIR = "resultados/exp04"

NUMERIC_FEATURES = [
    "SIFT_score", "PolyPhen_score", "SIFT_missing", "PolyPhen_missing",
    "VAF", "t_depth", "n_depth", "indel_len",
]
CATEGORICAL_FEATURES = ["Variant_Type"]
FEATURE_COLS = NUMERIC_FEATURES + CATEGORICAL_FEATURES


def parse_score(val):
    """Extrae el numero de 'deleterious(0.02)' / 'benign(0.99)'."""
    if not val:
        return np.nan
    m = re.search(r"\(([\d.]+)\)", val)
    if m:
        return float(m.group(1))
    return np.nan


def variant_type_and_indel_len(ref, alt):
    ref_len = len(ref.replace("-", ""))
    alt_len = len(alt.replace("-", ""))
    if ref_len == 1 and alt_len == 1:
        return "SNP", 0
    elif ref_len > alt_len:
        return "DEL", abs(alt_len - ref_len)
    elif ref_len < alt_len:
        return "INS", abs(alt_len - ref_len)
    else:
        return "ONP", abs(alt_len - ref_len)


def main():
    print(f"Leyendo VCF anotado: {VCF_PATH}", file=sys.stderr)

    csq_fields = None
    rows = []

    with gzip.open(VCF_PATH, "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: ([^"]+)"', line)
                csq_fields = m.group(1).split("|")
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info, fmt = fields[:9]
            tumor_sample = fields[9]

            # Solo variantes PASS (consistentes con somatic.filtered.vcf.gz)
            if filt != "PASS":
                continue

            # --- FORMAT (tumor sample) ---
            fmt_keys = fmt.split(":")
            tumor_vals = tumor_sample.split(":")
            fmt_dict = dict(zip(fmt_keys, tumor_vals))

            try:
                vaf = float(fmt_dict.get("AF", "nan").split(",")[0])
            except ValueError:
                vaf = np.nan
            try:
                t_depth = float(fmt_dict.get("DP", "nan"))
            except ValueError:
                t_depth = np.nan

            n_depth = np.nan
            if len(fields) > 10:
                normal_vals = fields[10].split(":")
                normal_dict = dict(zip(fmt_keys, normal_vals))
                try:
                    n_depth = float(normal_dict.get("DP", "nan"))
                except ValueError:
                    n_depth = np.nan

            # --- CSQ (VEP annotation) ---
            sift_score, polyphen_score = np.nan, np.nan
            consequence, symbol = None, None
            csq_match = re.search(r"CSQ=([^;\t]+)", info)
            if csq_match and csq_fields:
                transcripts = csq_match.group(1).split(",")
                chosen = None
                for t in transcripts:
                    vals = t.split("|")
                    d = dict(zip(csq_fields, vals))
                    if d.get("CANONICAL") == "YES":
                        chosen = d
                        break
                if chosen is None and transcripts:
                    chosen = dict(zip(csq_fields, transcripts[0].split("|")))
                if chosen:
                    sift_score = parse_score(chosen.get("SIFT", ""))
                    polyphen_score = parse_score(chosen.get("PolyPhen", ""))
                    consequence = chosen.get("Consequence")
                    symbol = chosen.get("SYMBOL")

            # --- Variant type / indel length ---
            vtype, indel_len = variant_type_and_indel_len(ref, alt)

            rows.append({
                "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt,
                "Gene": symbol, "Consequence": consequence,
                "SIFT_score": sift_score, "PolyPhen_score": polyphen_score,
                "SIFT_missing": int(np.isnan(sift_score)),
                "PolyPhen_missing": int(np.isnan(polyphen_score)),
                "VAF": vaf, "t_depth": t_depth, "n_depth": n_depth,
                "Variant_Type": vtype, "indel_len": indel_len,
            })

    df = pd.DataFrame(rows)
    print(f"Variantes PASS procesadas: {len(df):,}", file=sys.stderr)
    print(f"Variantes con VAF/t_depth utilizables: "
          f"{df[['VAF','t_depth']].notna().all(axis=1).sum():,}", file=sys.stderr)

    # Eliminar filas sin VAF/t_depth (no se pueden puntuar)
    df_model = df.dropna(subset=["VAF", "t_depth"]).copy()
    X = df_model[FEATURE_COLS]

    for name in ["RandomForest", "XGBoost"]:
        model_path = f"{MODEL_DIR}/model_{name}.joblib"
        print(f"Cargando modelo {name} desde {model_path}", file=sys.stderr)
        pipe = joblib.load(model_path)
        proba = pipe.predict_proba(X)[:, 1]
        df_model[f"proba_relevante_{name}"] = proba

    out_path = f"{OUT_DIR}/hcc1395_annotated_with_ml.tsv"
    df_model.to_csv(out_path, sep="\t", index=False)
    print(f"\nGuardado: {out_path} ({len(df_model):,} variantes)", file=sys.stderr)

    # Resumen
    summary = {
        "n_variants_pass": int(len(df)),
        "n_variants_scored": int(len(df_model)),
        "mean_proba_RandomForest": float(df_model["proba_relevante_RandomForest"].mean()),
        "mean_proba_XGBoost": float(df_model["proba_relevante_XGBoost"].mean()),
        "n_high_relevance_RF_gt_0.5": int((df_model["proba_relevante_RandomForest"] > 0.5).sum()),
        "n_high_relevance_XGB_gt_0.5": int((df_model["proba_relevante_XGBoost"] > 0.5).sum()),
        "top_genes_high_relevance_RF": (
            df_model[df_model["proba_relevante_RandomForest"] > 0.9]["Gene"]
            .value_counts().head(15).to_dict()
        ),
    }
    with open(f"{OUT_DIR}/hcc1395_ml_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== RESUMEN ===", file=sys.stderr)
    print(json.dumps(summary, indent=2, default=str), file=sys.stderr)


if __name__ == "__main__":
    main()
