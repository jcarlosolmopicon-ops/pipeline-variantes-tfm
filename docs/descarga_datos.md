# Guía de descarga de datos

## Exp01 — SEQC2 / HCC1395 (validación del pipeline)

### Requisitos
- sra-tools >= 3.1.1 (`conda install -c bioconda sra-tools`)
- ~130 GB de espacio libre

### Descarga FASTQ completos
```bash
# Tumor — HCC1395 (WES, ~64 GB)
prefetch SRR7890824 --output-directory datos-raw/HCC1395 --max-size 100GB
fasterq-dump datos-raw/HCC1395/SRR7890824 --split-files --threads 4
gzip datos-raw/HCC1395/SRR7890824_*.fastq

# Normal — HCC1395BL (WES, ~50 GB)
prefetch SRR7890827 --output-directory datos-raw/HCC1395 --max-size 100GB
fasterq-dump datos-raw/HCC1395/SRR7890827 --split-files --threads 4
gzip datos-raw/HCC1395/SRR7890827_*.fastq
```

### Truth set (variantes de referencia)
```bash
# Descargar high-confidence call set del consorcio SEQC2
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/somatic/HCC1395_HCC1395BL/ \
     -r -l1 --no-parent -A "*.vcf.gz"
```

## Exp02 — TCGA-LUAD (módulo ML)

### Opción A: TCGAbiolinks en R (GRCh38, recomendado)
```r
library(TCGAbiolinks)
maf <- GDCquery_Maf("LUAD", pipelines = "mutect2")
write.table(maf,
  "datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf",
  sep="\t", quote=FALSE, row.names=FALSE)
```

### Opción B: wget directo (GRCh37, legacy)
```bash
wget "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" \
     -O datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh37.legacy.maf.gz
```

## Genoma de referencia GRCh38
```bash
mkdir -p datos-raw/referencia
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
     -O datos-raw/referencia/hg38.fa.gz
gunzip datos-raw/referencia/hg38.fa.gz

# Indexar para BWA-MEM2
bwa-mem2 index datos-raw/referencia/hg38.fa
samtools faidx datos-raw/referencia/hg38.fa
```

## Notas importantes
- Los ficheros FASTQ, BAM y MAF NO se versionan en Git (.gitignore)
- Los ficheros de referencia tampoco se versionan por tamaño
- Documentar fecha de descarga y MD5 en cada experimento
