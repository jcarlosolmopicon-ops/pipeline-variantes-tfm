# Lab Notebook — Pipeline Variantes Somaticas
## Registro tecnico de sesiones de trabajo
## Repo: ~/pipeline-variantes | Maquina: DESKTOP-QK2T4QJ (WSL2 Ubuntu)

---

## Entrada NB-001 — Setup inicial del entorno de trabajo

**Fecha:** 24 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Instalacion y configuracion del stack completo de herramientas bioinformaticas
en Windows 11 con WSL2 (Ubuntu). Docker Desktop ya estaba instalado y configurado
con integracion WSL2.

### Regla de trabajo en WSL2

Todo el proyecto reside en el home de WSL2 (~/) y nunca en /mnt/c/.
Los archivos en /mnt/c/ son lentos y pueden causar errores en herramientas
bioinformaticas que asumen inodos POSIX.

### Pasos ejecutados

    # Verificacion Docker
    docker run hello-world

    # Instalacion Java 17
    sudo apt install -y openjdk-17-jdk curl wget git unzip

    # Instalacion Nextflow
    curl -s https://get.nextflow.io | bash
    mkdir -p ~/bin && mv nextflow ~/bin/
    nextflow -version

    # Instalacion Miniconda (Linux x86_64 dentro de WSL2)
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

    # Configuracion de canales conda
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    # Creacion del entorno tfm-variantes
    conda env create -f entorno/environment.yml
    conda activate tfm-variantes

### Versiones instaladas

| Herramienta | Version |
|-------------|---------|
| Nextflow    | 25.10.4 |
| Python      | 3.11    |
| Java        | 17      |
| Docker      | Desktop integrado con WSL2 |

### Decisiones tecnicas tomadas

- GRCh38 como genoma de referencia (GRCh37 obsoleto).
- Nextflow DSL2 sobre Snakemake: integracion nativa Docker/Singularity
  y disponibilidad de nf-core como referencia.
- BWA-MEM2 sobre BWA-MEM: mismos resultados con velocidad 2x.
- VEP v110 sobre ANNOVAR: sistema de plugins mas extensible (CADD, SpliceAI).
- RF + XGBoost sobre deep learning: interpretabilidad clinica y robustez
  ante datos desbalanceados.
- Conda (bioconda) sobre pip+venv: gestiona dependencias no Python (C libs, Java).

### Problemas encontrados y soluciones

- Docker no disponible en WSL: Docker Desktop no estaba en ejecucion.
  Solucion: abrir Docker Desktop primero.
- Permission denied con Docker:
  sudo usermod -aG docker $USER && newgrp docker

---

## Entrada NB-002 — Inicializacion Git y configuracion Nextflow

**Fecha:** 24 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Inicializacion del repositorio Git en ~/pipeline-variantes, construccion
de la estructura de directorios del proyecto y configuracion de Nextflow
con perfiles local (Docker) y conda_local. Verificacion con smoke test.

### Comandos ejecutados

    # Inicializacion del repo
    mkdir -p ~/pipeline-variantes && cd ~/pipeline-variantes
    git init
    git config user.name "Juan Carlos Olmo Picon"
    git config user.email "tu@email.com"

    # Estructura de directorios
    mkdir -p config datos-intermedios/nextflow-work datos-raw/{HCC1395,TCGA-LUAD,referencia,tests}
    mkdir -p docs entorno resultados/{exp01,exp02} scripts
    find . -type d -empty -exec touch {}/.gitkeep \;

    # Smoke test Nextflow
    nextflow run hello -profile local

### Commits realizados

    494cc0d config: add nextflow.config with local and conda_local profiles
    0121785 chore: add .gitignore for nextflow, data files and python cache
    a347968 chore: add conda environment.yml v1.0 with bioinformatics and ML deps
    22cd5b1 config: add nextflow.config profiles and hpc.config for slurm
    ebacf1a chore: init repo directory structure with .gitkeep placeholders
    956ab11 docs: add README.md with project description and repo structure

### nextflow.config (perfil local y conda_local)

    profiles {
        local {
            process.executor  = 'local'
            docker.enabled    = true
            docker.runOptions = '--platform linux/amd64'
        }
        conda_local {
            process.executor = 'local'
            conda.enabled    = true
        }
    }
    workDir        = 'datos-intermedios/nextflow-work'
    conda.cacheDir = "$HOME/.conda-nf-cache"
    includeConfig  'config/hpc.config'

### Problemas encontrados y soluciones

- Error "Unknown configuration profile: local": nextflow.config no existia
  en la raiz del repo. Solucion: crear nextflow.config en ~/pipeline-variantes/
- Fichero "reword" creado por accidente durante git rebase -i.
  Solucion: rm reword antes del commit.
- Mensajes de commit sin prefijo en los primeros 3 commits.
  Solucion: git rebase -i HEAD~3 usando reword para corregirlos.

---

## Entrada NB-003 — Seleccion de dataset y descarga de datos

**Fecha:** 24 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Seleccion del dataset de trabajo y descarga de los ficheros SRA y MAF.
Se cambio el enfoque de variantes germinales a variantes somaticas (tumor vs
normal) por mayor relevancia clinica y disponibilidad de truth set de referencia.

### Decision de dataset

- Exp01 (validacion pipeline): SEQC2/HCC1395 — par tumor/normal con truth set
  de alta confianza (39.536 SNVs + 2.020 INDELs HighConf).
- Exp02 (modulo ML): TCGA-LUAD MAF — 585 casos de adenocarcinoma de pulmon.
- Descartado: TCGA para datos crudos BAM/FASTQ por requerir acceso dbGaP.

### Comandos de descarga SRA

    # Actualizacion sra-tools (v2.9.6 fallaba con error SSL)
    conda install sra-tools=3.1.1

    # Descarga tumor HCC1395 (SRR7890824, ~64 GB)
    prefetch SRR7890824 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

    # Descarga normal HCC1395BL (SRR7890827, ~70 GB)
    prefetch SRR7890827 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

### Descarga y liftover TCGA-LUAD MAF

    # Descarga MAF GRCh37 legacy
    wget "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" \
      -O /mnt/k/TFM-bioinformatica/datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh37.legacy.maf.gz

    # Descarga chain file para liftover
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg19ToHg38.over.chain.gz

    # Liftover GRCh37 -> GRCh38 con CrossMap v0.7.3
    CrossMap.py maf hg19ToHg38.over.chain.gz \
      TCGA-LUAD.mutect2.GRCh37.legacy.maf.gz \
      /mnt/k/.../referencia/hg38.fa \
      TCGA-LUAD.mutect2.GRCh38.maf

### Resultado liftover

| Metrica | Valor |
|---------|-------|
| Variantes entrada (GRCh37) | 3.600.963 |
| Variantes salida (GRCh38)  | 3.600.605 |
| Perdidas en liftover       | 358 (0.01%) |

### Problemas encontrados y soluciones

- SRA-tools v2.9.6 error SSL (mbedtls_ssl_handshake -9984).
  Solucion: conda install sra-tools=3.1.1
- prefetch rechazaba ficheros >20 GB.
  Solucion: anadir --max-size 100GB
- MAF GDC era GRCh37 legacy. GDC GRCh38 ya no distribuye MAF agregado.
  Solucion: liftover con CrossMap v0.7.3
- TCGAbiolinks: conflicto de dependencias con Conda.
  Solucion: descartado en favor de CrossMap.
- Lock file SRR7890827.sra.lock al cerrar terminal.
  Solucion: rm *.sra.lock y relanzar prefetch (reanuda automaticamente).

### Commits generados

    004ca6d docs: add dataset README and data download guide for HCC1395 and TCGA-LUAD
    a8c5f7b docs: add liftover procedure GRCh37 to GRCh38 for TCGA-LUAD MAF
    fe67618 docs: add storage strategy for large files on external drive K
    7e1b21e chore: update conda environment to sra-tools=3.1.1 and export explicit lock

---

## Entrada NB-004 — Estrategia de almacenamiento en disco externo K:

**Fecha:** 24 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

El disco C: (WSL2) solo disponia de ~67 GB libres, insuficiente para los
FASTQs descomprimidos (~422 GB). Se adopto una estrategia de almacenamiento
en disco externo K: (1.9 TB) para todos los datos pesados del proyecto.

### Inventario de discos

| Disco | Montaje WSL2 | Libre | Uso |
|-------|-------------|-------|-----|
| C: (interno) | ~ (home WSL2) | ~67 GB | Codigo, conda, Nextflow workDir |
| D: (interno) | /mnt/d | ~319 GB | Reserva |
| K: (externo USB) | /mnt/k | ~890 GB | Datos pesados del TFM |

### Montar disco K: en WSL2

    sudo mkdir -p /mnt/k
    sudo mount -t drvfs K: /mnt/k

Este comando debe ejecutarse en cada inicio de WSL2 ya que el montaje
no es persistente entre sesiones.

### Estructura de datos en K:

    /mnt/k/TFM-bioinformatica/
    ├── datos-raw/
    │   ├── HCC1395/          # SRA y FASTQs tumor/normal
    │   ├── TCGA-LUAD/        # MAF GRCh37 y GRCh38
    │   └── referencia/       # GRCh38, indices, VEP cache, chain files
    ├── datos-intermedios/    # Temporales fasterq-dump (se limpian solos)
    └── resultados/           # Outputs del pipeline

### Problemas documentados

- FASTQs parciales (60 GB x2) generados por fasterq-dump en C:.
  Solucion: eliminados con rm; relanzar con --outdir apuntando a K:.
- Temporales fasterq-dump de 315 GB consumieron espacio de K:.
  Solucion: rm -rf datos-intermedios/* y relanzar.

---

## Entrada NB-005 — Conversion FASTQ tumor HCC1395

**Fecha:** 24-25 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Conversion de SRR7890824.sra (tumor HCC1395, 65 GB) a formato FASTQ mediante
fasterq-dump. El disco externo K: fue identificado como cuello de botella
principal de I/O al tener temporales y salida en el mismo disco.

### Comando ejecutado

    fasterq-dump /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/SRR7890824.sra \
      --split-files \
      --threads 2 \
      --mem 8GB \
      --temp /mnt/k/TFM-bioinformatica/datos-intermedios/ \
      --outdir /mnt/k/TFM-bioinformatica/datos-raw/HCC1395

### Verificacion e integridad

    # Verificar ficheros generados
    ls -lh /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824*.fastq*

    # Test de integridad
    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq.gz && echo "R1 OK"
    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq.gz && echo "R2 OK"

    # Borrar SRA original tras verificacion
    rm -rf /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/

### Metricas de ejecucion

| Metrica              | Valor                            |
|----------------------|----------------------------------|
| Inicio               | 24/04/2026 ~18:00                |
| Fin                  | 25/04/2026 ~10:44                |
| Duracion total       | ~16h 44min                       |
| SRA entrada          | 65 GB                            |
| FASTQs sin comprimir | ~232 GB (R1: 116 GB, R2: 116 GB) |
| FASTQs comprimidos   | ~95 GB (R1: 45 GB, R2: 50 GB)   |
| Temporales maximos   | ~401 GB                          |
| Ratio compresion     | ~2.4x                            |

### Problemas encontrados y soluciones

1. Cuello de botella I/O: temporales y salida en el mismo disco K:.
   Solucion: mover temporales a C: en futuras conversiones (ver NB-006).

2. gzip Operation not permitted en disco K: montado via drvfs.
   Solucion: usar gzip -c con redireccion stdout al fichero .gz.
   En la practica el primer gzip funciono y borro los originales.

3. Temporales fasterq-dump: techo de ~401 GB durante fase merge.
   Se limpian automaticamente al finalizar.

### Artefactos generados

- SRR7890824_1.fastq.gz — 45 GB — Verificado con gzip -t
- SRR7890824_2.fastq.gz — 50 GB — Verificado con gzip -t
- SRR7890824.sra borrado tras verificacion

---

## Entrada NB-006 — Optimizacion almacenamiento: compactacion VHDX WSL2

**Fecha:** 25 de abril de 2026
**Maquina:** DESKTOP-QK2T4QJ (Windows 11 + WSL2)
**Estado:** Completado

### Contexto

El disco virtual de WSL2 (ext4.vhdx) habia crecido hasta 583 GB en C:
por los ficheros temporales de fasterq-dump, dejando solo 67 GB libres.
Se compacto el VHDX y se adopto nueva arquitectura de discos separando
temporales (C:) de datos finales (K:).

### Diagnostico previo

    # Uso real del filesystem WSL2
    df -h /
    # Resultado: 42 GB usados de 1007 GB

    # Mayor consumidor
    du -sh ~/*
    # Resultado: 15 GB miniconda3

### Limpieza previa a la compactacion

    conda clean --all -y          # libero 3.38 GB
    docker system prune -f        # libero 1.115 GB
    rm ~/Miniconda3-latest-Linux-x86_64.sh   # 155 MB
    rm ~/ncbi_error_report.txt

    # Rellenar espacio libre con ceros para maximizar compactacion
    sudo dd if=/dev/zero of=/zeros.tmp bs=1M status=progress 2>/dev/null
    sudo rm /zeros.tmp

### Compactacion VHDX desde PowerShell Admin

    # Ruta del VHDX encontrada en:
    # C:\Users\Juan Carlos\AppData\Local\wsl\{4559b619-56e6-48c1-b26b-a75229854874}\ext4.vhdx

    wsl --shutdown

    diskpart
    select vdisk file="C:\Users\Juan Carlos\AppData\Local\wsl\{4559b619-56e6-48c1-b26b-a75229854874}\ext4.vhdx"
    attach vdisk readonly
    compact vdisk
    detach vdisk
    exit

    # Verificacion resultado
    Get-PSDrive C | Select-Object @{N='Libre(GB)';E={[math]::Round($_.Free/1GB,0)}}

### Metricas

| Metrica            | Antes  | Despues |
|--------------------|--------|---------|
| ext4.vhdx tamano   | 583 GB | ~30 GB  |
| C: libre           | 67 GB  | 619 GB  |
| Espacio recuperado | —      | +552 GB |

### Problemas encontrados y soluciones

- Optimize-VHD no disponible: requiere Hyper-V (no disponible en Windows Home).
  Solucion: usar diskpart con compact vdisk.
- diskpart "archivo en uso": WSL2 no estaba completamente apagado.
  Solucion: reiniciar Windows completamente antes de ejecutar diskpart.

### Nueva arquitectura de discos adoptada

| Disco | Uso | Libre |
|-------|-----|-------|
| C: interno | Temporales fasterq-dump (~/tmp-fasterq/), conda, codigo | 619 GB |
| K: externo | FASTQs .gz, referencias, resultados permanentes | ~544 GB |

La separacion de I/O entre discos reduce el tiempo estimado de conversion
del normal en ~50% respecto al tumor.

### Proximo paso

Conversion SRR7890827 (muestra normal, 70 GB) con temporales en C::

    mkdir -p ~/tmp-fasterq

    fasterq-dump /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890827/SRR7890827.sra \
      --split-files \
      --threads 4 \
      --mem 14GB \
      --temp ~/tmp-fasterq/ \
      --outdir /mnt/k/TFM-bioinformatica/datos-raw/HCC1395

---

## NB-XX — Validacion con hap.py contra truth set SEQC2 (HCC1395)

### Objetivo
Evaluar precision/recall del pipeline (MuTect2 + filtros) frente al
truth set HighConf/MedConf de SEQC2 (39.536 SNVs + 2.020 INDELs).

### Entorno
- hap.py no disponible via conda/pkrusche (imagen Docker obsoleta,
  manifest v1 no soportado por containerd >= 2.1)
- Usada imagen alternativa: jmcdani20/hap.py:v0.3.12

### Preparacion del truth set
El truth set SEQC2 (highconf_sSNV.vcf.gz + highconf_sINDEL.vcf.gz)
no es directamente compatible con hap.py por dos motivos:

1. **Sin columna FORMAT/sample**: hap.py requiere genotipos.
   Solucion: anadir sample dummy "TRUTH" con GT.
2. **FILTER invalido**: el campo FILTER contiene valores compuestos
   ("PASS;HighConf", "PASS;MedConf"), no validos en VCF estandar.
   hap.py interpreta estos registros como "no PASS" y los descarta
   (TRUTH.TOTAL=0 en primer intento).

Pasos aplicados:
    # 1. Concatenar SNV + INDEL truth sets
    bcftools concat -a highconf_sSNV.vcf.gz highconf_sINDEL.vcf.gz \
      -O z -o truth.vcf.gz

    # 2. Normalizar FILTER a PASS (todas las variantes son
    #    HighConf o MedConf, ambas validas)
    zcat truth.vcf.gz | awk 'BEGIN{OFS="\t"}
      /^#/{print; next}
      {$7="PASS"; print}' | bgzip > truth.fixed.vcf.gz

    # 3. Reconstruir cabecera: contigs completos (tomados del VCF
    #    de mutect2, que ya tiene ##contig de GRCh38) + ##FORMAT GT
    zcat somatic.filtered.vcf.gz | grep "^##contig" > contigs.txt
    {
      grep "^##fileformat" truth_header_noctig.txt
      cat contigs.txt
      grep "^##" truth_header_noctig.txt | grep -v "^##fileformat"
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    } > new_header.txt

    # 4. Anadir sample TRUTH con GT 0/1 (heterocigoto, esperado
    #    para variantes somaticas; GT 1/1 dio TP=0)
    {
      cat new_header.txt
      zcat truth.fixed.vcf.gz | awk 'BEGIN{OFS="\t"}
        /^#CHROM/{print $0,"FORMAT","TRUTH"; next}
        /^#/{next}
        {print $0,"GT","0/1"}'
    } | bgzip > truth.gt.vcf.gz
    tabix -p vcf truth.gt.vcf.gz

### Comando hap.py
    docker run -it --rm \
      -v ~/pipeline-variantes:/data \
      jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
      /data/resultados/exp03/happy/truth.gt.vcf.gz \
      /data/resultados/exp03/mutect2/somatic.filtered.vcf.gz \
      -f /data/datos-raw/HCC1395/truth_set/highconf_regions.bed \
      -r /data/datos-raw/referencia/GRCh38.fa \
      -o /data/resultados/exp03/happy/seqc2_eval \
      --pass-only

### Resultados (somatic.filtered.vcf.gz vs SEQC2 truth set)
| Tipo  | TRUTH.TOTAL | TP    | FN   | Recall | Precision | F1    |
|-------|-------------|-------|------|--------|-----------|-------|
| SNP   | 39447       | 35453 | 3994 | 0.8988 | 0.9781    | 0.9367|
| INDEL | 1625        | 1412  | 213  | 0.8689 | 0.6054    | 0.7136|

### Interpretacion
- SNVs: alta precision (97.8%) y buen recall (89.9%) -> resultado
  solido, comparable a benchmarks publicados de MuTect2 en HCC1395.
- INDELs: recall aceptable (86.9%) pero precision baja (60.5%),
  esperado en MuTect2 sin filtros especificos de indels adicionales
  (ej. realineacion local agresiva, panel of normals).
- Frac_NA alta (54-79%) refleja variantes del query fuera de las
  regiones highconf_regions.bed (zona no evaluable, no FP reales).

### Salidas generadas
resultados/exp03/happy/
  seqc2_eval.summary.csv
  seqc2_eval.extended.csv
  seqc2_eval.vcf.gz (+ .tbi)
  seqc2_eval.roc.*.csv.gz
  seqc2_eval.runinfo.json
  truth.gt.vcf.gz (+ .tbi)  # truth set corregido, reproducible

### Proximo paso
Redactar seccion 5 (Resultados) y 6 (Conclusiones) del TFM con
estas metricas. Decidir alcance del modulo ML (RF/XGBoost, Exp02)
segun tiempo disponible hasta deposito (20 junio).
