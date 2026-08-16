#!/usr/bin/env python3
"""
Version parametrizada de apply_model_hcc1395.py para su uso como proceso de
Nextflow (modulo CLASSIFY, experimento exp08).

Aplica los modelos RF/XGBoost serializados a un VCF anotado por VEP. La
logica de extraccion es identica a la del script original: bloque CSQ con
transcrito CANONICAL (o el primero en su defecto), VAF y profundidades de
los campos FORMAT de MuTect2, tipo de variante y longitud de indel a partir
de REF/ALT. Solo cambian las rutas, que se reciben por linea de comandos.

Uso:
    python3 apply_model.py --vcf annotated.vcf.gz \
        --modelos <dir con model_*.joblib> --outdir <dir de salida>
"""

import argparse
import re
import sys
import gzip
import json
import os
import joblib
import numpy as np
import pandas as pd

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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="VCF anotado por VEP (bgzip)")
    parser.add_argument("--modelos", required=True,
                        help="Directorio con model_RandomForest.joblib y model_XGBoost.joblib")
    parser.add_argument("--outdir", required=True, help="Directorio de salida")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Leyendo VCF anotado: {args.vcf}", file=sys.stderr)

    csq_fields = None
    rows = []

    with gzip.open(args.vcf, "rt") as f:
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
        model_path = f"{args.modelos}/model_{name}.joblib"
        print(f"Cargando modelo {name} desde {model_path}", file=sys.stderr)
        pipe = joblib.load(model_path)
        # Prediccion secuencial: los modelos se serializaron con n_jobs=-1 y
        # la acumulacion multihilo de los arboles suma los flotantes en orden
        # no determinista, variando el ultimo decimal entre ejecuciones. Con
        # un solo hilo la salida es bit-identica en cada ejecucion.
        if hasattr(pipe, "named_steps") and "clf" in pipe.named_steps:
            clf = pipe.named_steps["clf"]
            if hasattr(clf, "n_jobs"):
                clf.n_jobs = 1
        proba = pipe.predict_proba(X)[:, 1]
        df_model[f"proba_relevante_{name}"] = proba

    out_path = f"{args.outdir}/hcc1395_annotated_with_ml.tsv"
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
    with open(f"{args.outdir}/hcc1395_ml_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== RESUMEN ===", file=sys.stderr)
    print(json.dumps(summary, indent=2, default=str), file=sys.stderr)


if __name__ == "__main__":
    main()
