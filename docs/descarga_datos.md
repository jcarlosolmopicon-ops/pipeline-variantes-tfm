# Guía de descarga de datos

Procedimiento para obtener los datos de partida del proyecto: el par tumor-normal
SEQC2/HCC1395, el MAF somático pan-cáncer del TCGA y los ficheros de referencia.

Las rutas son relativas a la raíz del repositorio. Los comandos se redactaron durante la
fase inicial del proyecto sobre WSL2 y se reescribieron a las rutas del repositorio tras la
migración a Ubuntu nativo (ver NB-007 del lab notebook).

---

## Requisitos previos

    sra-tools >= 3.1.1   # v2.9.6 falla con error SSL
    crossmap  >= 0.7.2
    conda activate tfm-variantes

---

## HCC1395 / HCC1395BL — validación del pipeline (exp03)

### Paso 1 — Descarga SRA con prefetch

    # Tumor HCC1395 (~64 GB)
    prefetch SRR7890824 --output-directory datos-raw/HCC1395 --max-size 100GB

    # Normal HCC1395BL (~70 GB)
    prefetch SRR7890827 --output-directory datos-raw/HCC1395 --max-size 100GB

`--max-size 100GB` es obligatorio: sin él, prefetch rechaza ficheros de más de 20 GB.
Si la descarga se interrumpe, relanzar el mismo comando — prefetch reanuda automáticamente.
Si queda un fichero `.sra.lock`, borrarlo antes de relanzar.

### Paso 2 — Conversión SRA a FASTQ con fasterq-dump

Los temporales alcanzan un pico de ~400 GB durante la fase de merge y se limpian solos al
terminar. Conviene comprobar el espacio libre antes de empezar.

    fasterq-dump datos-raw/HCC1395/SRR7890824/SRR7890824.sra \
      --split-files \
      --threads 4 \
      --mem 14GB \
      --outdir datos-raw/HCC1395

    fasterq-dump datos-raw/HCC1395/SRR7890827/SRR7890827.sra \
      --split-files \
      --threads 4 \
      --mem 14GB \
      --outdir datos-raw/HCC1395

### Paso 3 — Comprimir y verificar

    gzip datos-raw/HCC1395/SRR7890824_1.fastq
    gzip datos-raw/HCC1395/SRR7890824_2.fastq

    gzip -t datos-raw/HCC1395/SRR7890824_1.fastq.gz && echo "R1 OK"
    gzip -t datos-raw/HCC1395/SRR7890824_2.fastq.gz && echo "R2 OK"

    # Borrar el .sra original tras verificar (libera ~65 GB)
    rm -rf datos-raw/HCC1395/SRR7890824/

### Paso 4 — Descargar el truth set de SEQC2

    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/somatic/HCC1395_HCC1395BL/ \
      -r -l1 --no-parent -A "*.vcf.gz" \
      -P datos-raw/HCC1395/truth_set/

Se obtienen `highconf_sSNV.vcf.gz`, `highconf_sINDEL.vcf.gz` y `highconf_regions.bed`.
La adaptación del truth set para hap.py se documenta en NB-008 del lab notebook.

---

## MAF somático pan-cáncer del TCGA — módulo ML (exp04)

El fichero se descarga bajo la denominación del proyecto TCGA-LUAD, pero su contenido real
corresponde al conjunto pan-cáncer MC3 del TCGA (10.295 muestras tumorales de 33 tipos de
cáncer). Ver sección 5.3.1 de la memoria.

### Descarga del MAF GRCh37 legacy

    wget "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" \
      -O datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz

TCGAbiolinks se descartó por conflicto de dependencias con Conda.

### Liftover GRCh37 a GRCh38 con CrossMap v0.7.2

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O datos-raw/referencia/hg19ToHg38.over.chain.gz

    CrossMap maf \
      datos-raw/referencia/hg19ToHg38.over.chain.gz \
      datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz \
      datos-raw/referencia/GRCh38.fa \
      GRCh38 \
      datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf

De los 3.600.963 registros del fichero original se convierten 3.598.760 (99,94 %); los
2.203 restantes quedan en el fichero `.unmap`. Detalle en `docs/liftover_GRCh37_to_GRCh38.md`.

---

## Genoma de referencia GRCh38

    mkdir -p datos-raw/referencia

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
      -O datos-raw/referencia/GRCh38.fa.gz
    gunzip datos-raw/referencia/GRCh38.fa.gz

    # Índice BWA (~5 GB, 1-2 h). Se emplea BWA 0.7.17 clásico: BWA-MEM2 dio errores
    # reproducibles de memoria insuficiente en el hardware disponible (ver NB-007).
    bwa index datos-raw/referencia/GRCh38.fa

    # Índices para samtools y GATK
    samtools faidx datos-raw/referencia/GRCh38.fa
    gatk CreateSequenceDictionary -R datos-raw/referencia/GRCh38.fa

---

## Caché de VEP GRCh38 v110

    vep_install -a cf -s homo_sapiens -y GRCh38 --CACHEDIR ~/.vep

Ocupa ~15 GB. El pipeline la monta dentro del contenedor mediante `containerOptions`
(ver `nextflow.config`).

---

## dbSNP138 para BQSR

Se emplea dbSNP138 del bundle de recursos de GATK, por compatibilidad directa con las
herramientas del Broad Institute y con las coordenadas GRCh38 del resto del flujo.

    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf \
      -O datos-raw/referencia/dbsnp138.vcf

    gatk IndexFeatureFile -I datos-raw/referencia/dbsnp138.vcf

---

## ClinVar (etiquetas de patogenicidad, exp07)

Distribución en GRCh38, de acceso libre y sin registro. Se conserva el nombre con la fecha de
la versión para que quede trazable qué release se empleó, porque ClinVar se actualiza
semanalmente y las clasificaciones cambian.

    mkdir -p datos-raw/clinvar

    curl -L -o datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz \
      https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

    curl -L -o datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz.tbi \
      https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

    # dejar constancia de la version exacta
    zcat datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz | head -20 \
      | grep -E '^##(fileDate|source|reference)' > resultados/exp07/clinvar_version.txt

El fichero pesa 185 MB y está en `.gitignore`. La versión empleada tiene
`fileDate=2026-08-08`.


---

## Notas

- Los ficheros FASTQ, BAM, MAF y de referencia no se versionan en Git (ver `.gitignore`).
- CADD v1.7 (~80 GB) no llegó a instalarse: la anotación final no incluye CADD ni SpliceAI.
- Verificar siempre la integridad con `gzip -t` antes de borrar los originales.
