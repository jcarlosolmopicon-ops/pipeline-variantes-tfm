# Liftover MAF TCGA-LUAD de GRCh37 a GRCh38

## Contexto
El MAF legacy de TCGA-LUAD descargado del GDC está en GRCh37 (hg19).
El pipeline usa GRCh38, por lo que se realizó liftover de coordenadas.

## Herramienta
CrossMap v0.7.3 — https://crossmap.readthedocs.io

## Cadena de conversión
hg19ToHg38.over.chain.gz (UCSC Genome Browser)
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

## Resultados
- Variantes originales (GRCh37): 3.600.963
- Variantes convertidas (GRCh38): 3.600.605
- Variantes omitidas: 358 (0.01%)

## Comandos ejecutados
```bash
# 1. Extraer coordenadas del MAF original
zcat TCGA-LUAD.mutect2.somatic.maf.gz \
  | grep -v "^#" | tail -n +2 \
  | awk 'BEGIN{OFS="\t"} {print $5, $6-1, $7, NR}' \
  > /tmp/luad_coords.bed

# 2. Liftover con CrossMap
CrossMap bed hg19ToHg38.over.chain.gz \
  /tmp/luad_coords.bed \
  /tmp/luad_coords_hg38.bed

# 3. Reconstruir MAF con coordenadas GRCh38
python3 scripts/liftover_maf.py
```

## Ficheros
- Input:  datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz (GRCh37)
- Output: datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz (GRCh38)
- Cadena: datos-raw/referencia/hg19ToHg38.over.chain.gz
