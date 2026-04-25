# Guia de descarga de datos
## Actualizado: 25 de abril de 2026

---

## Requisitos previos

    # Versiones requeridas
    sra-tools >= 3.1.1   # IMPORTANTE: v2.9.6 falla con error SSL
    crossmap >= 0.7.3
    conda activate tfm-variantes

    # Montar disco externo K: antes de cualquier descarga
    sudo mkdir -p /mnt/k
    sudo mount -t drvfs K: /mnt/k

    # Directorio temporal en C: para fasterq-dump
    mkdir -p ~/tmp-fasterq

---

## Exp01 — SEQC2 / HCC1395 (validacion del pipeline)

### Paso 1 — Descarga SRA con prefetch

    # Tumor HCC1395 (~64 GB)
    prefetch SRR7890824 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

    # Normal HCC1395BL (~70 GB) — despues del tumor
    prefetch SRR7890827 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

Nota: --max-size 100GB es obligatorio. Sin el, prefetch rechaza ficheros >20 GB.
Si se interrumpe, relanzar el mismo comando — prefetch reanuda automaticamente.
Si aparece un .sra.lock, borrarlo con rm *.sra.lock antes de relanzar.

### Paso 2 — Conversion SRA a FASTQ con fasterq-dump

IMPORTANTE: Los temporales (~400 GB) deben ir a C: (~/tmp-fasterq/) y
la salida final a K:. Si ambos van al mismo disco, el proceso tarda ~17h
por contension de I/O.

    # Tumor
    START=$(date +%s)
    echo "Inicio: $(date '+%Y-%m-%d %H:%M:%S')"

    fasterq-dump /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/SRR7890824.sra \
      --split-files \
      --threads 4 \
      --mem 14GB \
      --temp ~/tmp-fasterq/ \
      --outdir /mnt/k/TFM-bioinformatica/datos-raw/HCC1395

    END=$(date +%s) && ELAPSED=$((END-START))
    echo "Fin: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Tiempo: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"

    # Normal (mismo comando cambiando SRR7890824 por SRR7890827)
    fasterq-dump /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890827/SRR7890827.sra \
      --split-files \
      --threads 4 \
      --mem 14GB \
      --temp ~/tmp-fasterq/ \
      --outdir /mnt/k/TFM-bioinformatica/datos-raw/HCC1395

Monitor de progreso (segunda terminal):

    watch -n 30 'df -h /mnt/k && echo "---" && df -h ~ && echo "---" \
      && du -sh ~/tmp-fasterq/ 2>/dev/null && echo "---" \
      && du -sh /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/'

Watchdog de espacio (tercera terminal):

    while true; do
      FREE=$(df ~ | awk 'NR==2 {print $4}')
      FREE_GB=$((FREE / 1024 / 1024))
      echo "$(date '+%H:%M:%S') — C: libre: ${FREE_GB} GB"
      if [ "$FREE" -lt 52428800 ]; then
        echo "ESPACIO CRITICO — matando fasterq-dump"
        pkill -f fasterq-dump
        break
      fi
      sleep 60
    done

### Paso 3 — Comprimir FASTQs

NOTA: En discos montados via drvfs usar gzip -c para evitar errores
de permisos al renombrar ficheros.

    # Comprimir en paralelo
    gzip -c /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq > \
      /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq.gz &

    gzip -c /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq > \
      /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq.gz &

    wait && echo "Compresion completada"

### Paso 4 — Verificar integridad y limpiar

    # Verificar integridad de los .fastq.gz
    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq.gz && echo "R1 OK"
    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq.gz && echo "R2 OK"

    # Borrar .fastq sin comprimir tras verificacion
    rm /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq
    rm /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq

    # Borrar .sra original tras verificacion (opcional, libera ~65 GB)
    rm -rf /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/

### Paso 5 — Descargar truth set SEQC2

    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/somatic/HCC1395_HCC1395BL/ \
      -r -l1 --no-parent -A "*.vcf.gz" \
      -P /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/

---

## Exp02 — TCGA-LUAD (modulo ML)

### Descarga MAF GRCh37 legacy (wget directo)

    wget "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" \
      -O /mnt/k/TFM-bioinformatica/datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh37.legacy.maf.gz

Nota: TCGAbiolinks descartado por conflicto de dependencias con Conda.

### Liftover GRCh37 a GRCh38 con CrossMap v0.7.3

    # Descargar chain file
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg19ToHg38.over.chain.gz

    # Ejecutar liftover
    CrossMap.py maf \
      /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg19ToHg38.over.chain.gz \
      /mnt/k/TFM-bioinformatica/datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh37.legacy.maf.gz \
      /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg38.fa \
      /mnt/k/TFM-bioinformatica/datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf

Resultado: 3.600.605 variantes (perdida 0.01% en liftover).

---

## Genoma de referencia GRCh38

    mkdir -p /mnt/k/TFM-bioinformatica/datos-raw/referencia

    # Descargar referencia
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
      -O /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg38.fa.gz

    gunzip /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg38.fa.gz

    # Indexar para BWA-MEM2 (~60 GB, puede tardar 1-2h)
    bwa-mem2 index /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg38.fa

    # Indexar para samtools/GATK
    samtools faidx /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg38.fa

---

## Cache VEP GRCh38 v110

    vep_install -a cf -s homo_sapiens -y GRCh38 \
      --CACHEDIR /mnt/k/TFM-bioinformatica/datos-raw/referencia/vep_cache

Tamano: ~15 GB. Puede tardar varias horas segun el ancho de banda.

---

## dbSNP b156 (para GATK BQSR)

    wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz \
      -O /mnt/k/TFM-bioinformatica/datos-raw/referencia/dbsnp_156.vcf.gz

    # Indexar
    tabix -p vcf /mnt/k/TFM-bioinformatica/datos-raw/referencia/dbsnp_156.vcf.gz

---

## Notas importantes

- Los ficheros FASTQ, BAM y MAF NO se versionan en Git (.gitignore).
- Los ficheros de referencia tampoco se versionan por tamano.
- CADD v1.7 (80 GB) — usar API online en lugar de descarga local.
- Documentar fecha de descarga y verificacion en el lab notebook.
- Siempre verificar integridad con gzip -t antes de borrar originales.
