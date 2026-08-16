# Lab Notebook — Pipeline de variantes somáticas

Registro técnico de las sesiones de trabajo del proyecto. El repositorio vive en
`~/pipeline-variantes`. El proyecto se desarrolló en dos entornos sucesivos: la fase inicial
sobre WSL2 (NB-001 a NB-006) y el resto sobre Ubuntu nativo tras la migración descrita en
NB-007.

---

## Entrada NB-001 — Setup inicial del entorno de trabajo

**Fecha:** 24 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Instalación y configuración del stack de herramientas bioinformáticas en Windows 11 con
WSL2 (Ubuntu). Docker Desktop ya estaba instalado con integración WSL2.

### Regla de trabajo en WSL2

Todo el proyecto reside en el home de WSL2 (`~/`) y nunca en `/mnt/c/`. Los archivos en
`/mnt/c/` son lentos y pueden causar errores en herramientas que asumen inodos POSIX.

### Pasos ejecutados

    # Verificación de Docker
    docker run hello-world

    # Java 17
    sudo apt install -y openjdk-17-jdk curl wget git unzip

    # Nextflow
    curl -s https://get.nextflow.io | bash
    mkdir -p ~/bin && mv nextflow ~/bin/
    nextflow -version

    # Miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

    # Canales de conda
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    # Entorno del proyecto
    conda env create -f entorno/environment.yml
    conda activate tfm-variantes

### Versiones instaladas

| Herramienta | Versión |
|-------------|---------|
| Nextflow    | 25.10.4 |
| Python      | 3.11    |
| Java        | 17      |
| Docker      | Desktop integrado con WSL2 |

### Decisiones técnicas tomadas

- GRCh38 como genoma de referencia (GRCh37 obsoleto).
- Nextflow DSL2 sobre Snakemake: integración nativa con Docker/Singularity y disponibilidad
  de nf-core como referencia.
- BWA-MEM2 sobre BWA-MEM: mismos resultados con velocidad 2x.
  **Revertido en NB-007** — BWA-MEM2 se descartó por consumo de memoria poco predecible.
- VEP v110 sobre ANNOVAR: sistema de plugins más extensible (CADD, SpliceAI).
  **Parcialmente revertido** — CADD y SpliceAI requieren ficheros de plugin externos que no
  llegaron a instalarse; la anotación final no los incluye.
- RF + XGBoost sobre deep learning: interpretabilidad y robustez ante datos desbalanceados.
- Conda (bioconda) sobre pip+venv: gestiona dependencias no Python (librerías C, Java).

### Problemas encontrados y soluciones

- Docker no disponible en WSL: Docker Desktop no estaba en ejecución.
  Solución: abrir Docker Desktop primero.
- `Permission denied` con Docker:
  `sudo usermod -aG docker $USER && newgrp docker`

---

## Entrada NB-002 — Inicialización de Git y configuración de Nextflow

**Fecha:** 24 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Inicialización del repositorio en `~/pipeline-variantes`, construcción de la estructura de
directorios y configuración de Nextflow con los perfiles `local` (Docker) y `conda_local`.
Verificación con smoke test.

### Comandos ejecutados

    # Inicialización del repositorio
    mkdir -p ~/pipeline-variantes && cd ~/pipeline-variantes
    git init
    git config user.name "Juan Carlos Olmo Picon"
    git config user.email "jcarlosolmopicon@gmail.com"

    # Estructura de directorios
    mkdir -p config datos-intermedios/nextflow-work datos-raw/{HCC1395,TCGA-LUAD,referencia,tests}
    mkdir -p docs entorno resultados/{exp01,exp02} scripts
    find . -type d -empty -exec touch {}/.gitkeep \;

    # Smoke test de Nextflow
    nextflow run hello -profile local

### Commits realizados

    494cc0d config: add nextflow.config with local and conda_local profiles
    0121785 chore: add .gitignore for nextflow, data files and python cache
    a347968 chore: add conda environment.yml v1.0 with bioinformatics and ML deps
    22cd5b1 config: add nextflow.config profiles and hpc.config for slurm
    ebacf1a chore: init repo directory structure with .gitkeep placeholders
    956ab11 docs: add README.md with project description and repo structure

### nextflow.config (versión inicial)

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

- Error `Unknown configuration profile: local`: `nextflow.config` no existía en la raíz del
  repositorio. Solución: crearlo en `~/pipeline-variantes/`.
- Fichero `reword` creado por accidente durante un `git rebase -i`.
  Solución: `rm reword` antes del commit.
- Mensajes de commit sin prefijo en los tres primeros commits.
  Solución: `git rebase -i HEAD~3` usando `reword`.

---

## Entrada NB-003 — Selección de dataset y descarga de datos

**Fecha:** 24 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Selección del dataset de trabajo y descarga de los ficheros SRA y MAF. Se cambió el enfoque
de variantes germinales a variantes somáticas (tumor frente a normal) por mayor relevancia
clínica y por la disponibilidad de un truth set de referencia.

### Decisión de dataset

- **Validación del pipeline:** SEQC2/HCC1395, par tumor-normal con truth set público.
- **Módulo ML:** MAF somático del TCGA descargado del GDC, en coordenadas GRCh37.
- Descartado: TCGA para datos crudos BAM/FASTQ, por requerir acceso dbGaP.

  **Corrección posterior (agosto de 2026).** El MAF se descargó bajo la denominación del
  proyecto TCGA-LUAD y así se registró en su momento, pero su contenido real corresponde al
  conjunto pan-cáncer MC3 del TCGA: 3.598.760 variantes de 10.295 muestras tumorales de 33
  tipos de cáncer, de las que el proyecto LUAD aporta el 6,1 %. La discrepancia se detectó
  durante la verificación final de los datos y se documenta en la sección 5.3.1 de la memoria.

  El recuento del truth set de SEQC2 registrado aquí (39.536 SNV + 2.020 INDEL) procede de
  la publicación de Fang et al. (2021). Los ficheros efectivamente descargados contienen
  39.560 SNV y 1.922 INDEL, de los que 37.398 y 1.754 son HighConf y el resto MedConf.

### Comandos de descarga SRA

    # Actualización de sra-tools (v2.9.6 fallaba con error SSL)
    conda install sra-tools=3.1.1

    # Tumor HCC1395 (SRR7890824, ~64 GB)
    prefetch SRR7890824 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

    # Normal HCC1395BL (SRR7890827, ~70 GB)
    prefetch SRR7890827 \
      --output-directory /mnt/k/TFM-bioinformatica/datos-raw/HCC1395 \
      --max-size 100GB

### Descarga y liftover del MAF

    # MAF GRCh37 legacy
    wget "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" \
      -O /mnt/k/TFM-bioinformatica/datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz

    # Chain file
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O /mnt/k/TFM-bioinformatica/datos-raw/referencia/hg19ToHg38.over.chain.gz

    # Liftover GRCh37 -> GRCh38
    CrossMap maf hg19ToHg38.over.chain.gz \
      TCGA-LUAD.mutect2.somatic.maf.gz \
      GRCh38.fa GRCh38 \
      TCGA-LUAD.mutect2.GRCh38.maf

### Resultado del liftover

| Métrica | Valor |
|---------|-------|
| Registros en el MAF original (GRCh37) | 3.600.963 |
| No convertidos (fichero `.unmap`) | 2.203 (0,061 %) |
| Convertidos a GRCh38 | 3.598.760 (99,94 %) |

La versión de CrossMap empleada en el liftover definitivo fue la **0.7.2**, registrada en la
línea de comentario que la propia herramienta antepone al fichero de salida. Ver
`docs/liftover_GRCh37_to_GRCh38.md`.

### Problemas encontrados y soluciones

- sra-tools v2.9.6 daba error SSL (`mbedtls_ssl_handshake -9984`).
  Solución: `conda install sra-tools=3.1.1`.
- prefetch rechazaba ficheros de más de 20 GB.
  Solución: añadir `--max-size 100GB`.
- El MAF del GDC era GRCh37 legacy; el GDC ya no distribuye el MAF agregado en GRCh38.
  Solución: liftover con CrossMap.
- TCGAbiolinks: conflicto de dependencias con Conda. Solución: descartado.
- Fichero `SRR7890827.sra.lock` al cerrar la terminal.
  Solución: `rm *.sra.lock` y relanzar prefetch, que reanuda automáticamente.

### Commits generados

    004ca6d docs: add dataset README and data download guide for HCC1395 and TCGA-LUAD
    a8c5f7b docs: add liftover procedure GRCh37 to GRCh38 for TCGA-LUAD MAF
    fe67618 docs: add storage strategy for large files on external drive K
    7e1b21e chore: update conda environment to sra-tools=3.1.1 and export explicit lock

---

## Entrada NB-004 — Estrategia de almacenamiento en disco externo

**Fecha:** 24 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado — superado por NB-007

### Contexto

El disco C: (WSL2) solo disponía de ~67 GB libres, insuficientes para los FASTQ
descomprimidos (~422 GB). Se adoptó una estrategia de almacenamiento en el disco externo K:
(1,9 TB) para todos los datos pesados del proyecto.

### Inventario de discos

| Disco | Montaje WSL2 | Libre | Uso |
|-------|-------------|-------|-----|
| C: (interno) | `~` (home WSL2) | ~67 GB | Código, conda, workDir de Nextflow |
| D: (interno) | `/mnt/d` | ~319 GB | Reserva |
| K: (externo USB) | `/mnt/k` | ~890 GB | Datos pesados del TFM |

### Montaje del disco K:

    sudo mkdir -p /mnt/k
    sudo mount -t drvfs K: /mnt/k

El montaje no es persistente: hay que repetirlo en cada inicio de WSL2.

### Estructura de datos en K:

    /mnt/k/TFM-bioinformatica/
    ├── datos-raw/
    │   ├── HCC1395/          # SRA y FASTQ tumor/normal
    │   ├── TCGA-LUAD/        # MAF GRCh37 y GRCh38
    │   └── referencia/       # GRCh38, índices, caché de VEP, chain files
    ├── datos-intermedios/    # Temporales de fasterq-dump
    └── resultados/           # Salidas del pipeline

### Problemas documentados

- FASTQ parciales (60 GB x2) generados por fasterq-dump en C:.
  Solución: eliminados; relanzar con `--outdir` apuntando a K:.
- Temporales de fasterq-dump de 315 GB consumieron espacio de K:.
  Solución: `rm -rf datos-intermedios/*` y relanzar.

---

## Entrada NB-005 — Conversión a FASTQ del tumor HCC1395

**Fecha:** 24-25 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (WSL2 Ubuntu)
**Estado:** Completado

### Contexto

Conversión de `SRR7890824.sra` (tumor HCC1395, 65 GB) a FASTQ mediante fasterq-dump. El
disco externo K: resultó ser el cuello de botella de E/S al tener temporales y salida en el
mismo disco.

### Comando ejecutado

    fasterq-dump /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/SRR7890824.sra \
      --split-files \
      --threads 2 \
      --mem 8GB \
      --temp /mnt/k/TFM-bioinformatica/datos-intermedios/ \
      --outdir /mnt/k/TFM-bioinformatica/datos-raw/HCC1395

### Verificación

    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_1.fastq.gz && echo "R1 OK"
    gzip -t /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824_2.fastq.gz && echo "R2 OK"

    # Borrar el SRA original tras verificar
    rm -rf /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/SRR7890824/

### Métricas de ejecución

| Métrica              | Valor                            |
|----------------------|----------------------------------|
| Inicio               | 24/04/2026 ~18:00                |
| Fin                  | 25/04/2026 ~10:44                |
| Duración total       | ~16 h 44 min                     |
| SRA de entrada       | 65 GB                            |
| FASTQ sin comprimir  | ~232 GB (R1: 116 GB, R2: 116 GB) |
| FASTQ comprimidos    | ~95 GB (R1: 45 GB, R2: 50 GB)    |
| Temporales máximos   | ~401 GB                          |
| Ratio de compresión  | ~2,4x                            |

### Problemas encontrados y soluciones

- Cuello de botella de E/S: temporales y salida en el mismo disco K:.
  Solución: mover los temporales a C: en las conversiones siguientes (ver NB-006).
- `gzip: Operation not permitted` en el disco K: montado vía drvfs.
  Solución: usar `gzip -c` con redirección al fichero `.gz`.
- Los temporales de fasterq-dump alcanzaron un techo de ~401 GB durante la fase de merge.
  Se limpian automáticamente al finalizar.

### Artefactos generados

- `SRR7890824_1.fastq.gz` — 45 GB — verificado con `gzip -t`
- `SRR7890824_2.fastq.gz` — 50 GB — verificado con `gzip -t`
- `SRR7890824.sra` borrado tras la verificación

---

## Entrada NB-006 — Compactación del disco virtual de WSL2

**Fecha:** 25 de abril de 2026
**Máquina:** DESKTOP-QK2T4QJ (Windows 11 + WSL2)
**Estado:** Completado — superado por NB-007

### Contexto

El disco virtual de WSL2 (`ext4.vhdx`) había crecido hasta 583 GB en C: por los ficheros
temporales de fasterq-dump, dejando solo 67 GB libres. Se compactó el VHDX y se adoptó una
nueva arquitectura de discos separando temporales (C:) de datos finales (K:).

### Diagnóstico previo

    df -h /          # 42 GB usados de 1007 GB
    du -sh ~/*       # mayor consumidor: 15 GB de miniconda3

### Limpieza previa a la compactación

    conda clean --all -y                     # liberó 3,38 GB
    docker system prune -f                   # liberó 1,115 GB
    rm ~/Miniconda3-latest-Linux-x86_64.sh   # 155 MB
    rm ~/ncbi_error_report.txt

    # Rellenar el espacio libre con ceros para maximizar la compactación
    sudo dd if=/dev/zero of=/zeros.tmp bs=1M status=progress 2>/dev/null
    sudo rm /zeros.tmp

### Compactación desde PowerShell con permisos de administrador

    wsl --shutdown

    diskpart
    select vdisk file="C:\Users\...\AppData\Local\wsl\{...}\ext4.vhdx"
    attach vdisk readonly
    compact vdisk
    detach vdisk
    exit

### Métricas

| Métrica            | Antes  | Después |
|--------------------|--------|---------|
| Tamaño de ext4.vhdx| 583 GB | ~30 GB  |
| C: libre           | 67 GB  | 619 GB  |
| Espacio recuperado | —      | +552 GB |

### Problemas encontrados y soluciones

- `Optimize-VHD` no disponible: requiere Hyper-V, ausente en Windows Home.
  Solución: usar `diskpart` con `compact vdisk`.
- `diskpart` daba "archivo en uso": WSL2 no estaba completamente apagado.
  Solución: reiniciar Windows antes de ejecutar diskpart.

### Nueva arquitectura de discos adoptada

| Disco | Uso | Libre |
|-------|-----|-------|
| C: interno | Temporales de fasterq-dump (`~/tmp-fasterq/`), conda, código | 619 GB |
| K: externo | FASTQ comprimidos, referencias, resultados permanentes | ~544 GB |

Separar la E/S entre discos redujo el tiempo de conversión de la muestra normal en
aproximadamente un 50 % respecto al tumor.

---

## Entrada NB-007 — Migración de WSL2 a Ubuntu nativo (arranque dual)

**Fecha:** 20 de mayo de 2026
**Máquina de origen:** DESKTOP-QK2T4QJ (Windows 11 + WSL2 Ubuntu)
**Máquina de destino:** olmop-MS-7E26 (Ubuntu nativo, arranque dual junto a Windows)
**Estado:** Completado

### Contexto

La migración se decidió durante la ejecución completa del pipeline sobre el par
HCC1395/HCC1395BL, no antes. La corrida se interrumpió de forma repetida por agotamiento de
memoria y por inestabilidad del entorno, lo que impedía completar las etapas más exigentes
(alineamiento y MuTect2) en una sola ejecución.

### Problemas acumulados en WSL2

1. **Gestión de memoria.** WSL2 opera sobre una máquina virtual ligera con memoria asignada
   dinámicamente. Bajo carga sostenida, BWA-MEM2 y GATK alcanzaban el límite disponible y
   morían por OOM, sin que el ajuste de `.wslconfig` lo resolviera de forma fiable.

2. **Rendimiento de E/S sobre drvfs.** Los datos pesados residían en el disco externo K:,
   montado mediante drvfs. El acceso a ficheros grandes a través de esa capa de traducción
   es sustancialmente más lento que sobre un sistema de ficheros nativo, y varias
   herramientas (sra-tools entre ellas) fallaban al intentar escribir en esa ruta.

3. **Crecimiento del disco virtual.** `ext4.vhdx` crece de forma dinámica pero no se compacta
   solo, lo que obligaba a intervenciones manuales periódicas (ver NB-006).

### Decisión

Instalar Ubuntu en arranque dual junto a Windows, en lugar de seguir ajustando la
configuración de WSL2.

### Pasos ejecutados

    # 1. Instalación de Ubuntu en partición propia (arranque dual)
    # 2. Reinstalación del stack: Java 17, Nextflow, Miniconda, Docker
    # 3. Recreación del entorno Conda
    conda env create -f entorno/environment.yml
    conda activate tfm-variantes

    # 4. Clonado del repositorio
    git clone https://github.com/jcarlosolmopicon-ops/pipeline-variantes-tfm.git

    # 5. Traslado de los datos pesados al almacenamiento nativo (ext4 sobre NVMe)
    #    y verificación de integridad tras la copia
    gzip -t datos-raw/HCC1395/*.fastq.gz

    # 6. Relanzamiento del pipeline
    nextflow run main.nf -profile local --step all -resume

### Resultado

La ejecución completa terminó sin interrupciones. Los tiempos mejoraron de forma apreciable
al eliminar la capa drvfs, aunque MuTect2 siguió siendo el cuello de botella del flujo (más
de 10 horas sobre el par completo).

### Consecuencias para el resto del proyecto

- La estrategia de almacenamiento de NB-004 y NB-006 (montaje de K:, temporales en C:,
  compactación del VHDX) deja de aplicar. Todos los datos residen a partir de aquí en el
  sistema de ficheros nativo.
- Se descartó definitivamente **BWA-MEM2** en favor de BWA 0.7.17 clásico. Aunque el cambio
  de entorno alivió la presión de memoria, BWA-MEM2 seguía mostrando un perfil de consumo
  poco predecible para el volumen de datos de este trabajo.
- Todos los resultados reportados en la memoria (exp03 a exp06) se obtuvieron en este entorno.

### Estrategia de almacenamiento vigente

Todo el proyecto reside en el sistema de ficheros nativo (ext4 sobre NVMe), sin capas de
traducción intermedias.

| Contenido | Ubicación |
|-----------|-----------|
| Código y entorno | `~/pipeline-variantes/`, `~/miniconda3/` |
| Datos de partida (FASTQ, referencias, MAF) | `datos-raw/` |
| Intermedios de Nextflow (workDir) | `datos-intermedios/` |
| Resultados | `resultados/expNN/` |

En Git se versionan el código, la configuración, la documentación y los ficheros de
resultados ligeros (métricas, informes, modelos serializados). Quedan excluidos por
`.gitignore` los ficheros pesados: FASTQ, BAM alineados y recalibrados, referencias e
índices, y el MAF de entrenamiento.

### Configuración final

| Elemento | Valor |
|----------|-------|
| Sistema | Ubuntu (arranque dual junto a Windows 11) |
| Hostname | olmop-MS-7E26 |
| CPUs asignadas al pipeline | 12 |
| Memoria asignada al pipeline | 28 GB |
| Almacenamiento | NVMe, ext4 nativo |
| Perfil de ejecución | `local` (Docker) |

---

## Entrada NB-008 — Validación con hap.py frente al truth set de SEQC2

**Fecha:** 11 de junio de 2026
**Máquina:** olmop-MS-7E26 (Ubuntu nativo)
**Estado:** Completado

### Objetivo

Evaluar precisión y recall del pipeline (MuTect2 + FilterMutectCalls) frente al truth set de
SEQC2 para el par HCC1395/HCC1395BL.

### Entorno

hap.py no está disponible vía conda y la imagen Docker de pkrusche está obsoleta (manifest
v1, no soportado por containerd >= 2.1). Se usó la imagen alternativa
`jmcdani20/hap.py:v0.3.12`.

### Preparación del truth set

El truth set de SEQC2 (`highconf_sSNV.vcf.gz` + `highconf_sINDEL.vcf.gz`) no es directamente
compatible con hap.py por dos motivos:

1. **Sin columna FORMAT/sample.** hap.py requiere genotipos.
   Solución: añadir una muestra ficticia `TRUTH` con campo GT.
2. **FILTER no estándar.** El campo contiene valores compuestos (`PASS;HighConf`,
   `PASS;MedConf`) que hap.py interpreta como "no PASS" y descarta, dando `TRUTH.TOTAL=0`
   en el primer intento.

Pasos aplicados:

    # 1. Concatenar los truth sets de SNV e INDEL
    bcftools concat -a highconf_sSNV.vcf.gz highconf_sINDEL.vcf.gz \
      -O z -o truth.vcf.gz

    # 2. Normalizar FILTER a PASS (HighConf y MedConf son ambas válidas)
    zcat truth.vcf.gz | awk 'BEGIN{OFS="\t"}
      /^#/{print; next}
      {$7="PASS"; print}' | bgzip > truth.fixed.vcf.gz

    # 3. Reconstruir la cabecera: contigs de GRCh38 (tomados del VCF de MuTect2)
    #    más la definición de FORMAT/GT
    zcat somatic.filtered.vcf.gz | grep "^##contig" > contigs.txt
    {
      grep "^##fileformat" truth_header_noctig.txt
      cat contigs.txt
      grep "^##" truth_header_noctig.txt | grep -v "^##fileformat"
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    } > new_header.txt

    # 4. Añadir la muestra TRUTH con GT 0/1 (heterocigoto, coherente con la naturaleza
    #    somática de las variantes; con GT 1/1 el resultado fue TP=0)
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

### Resultados

| Tipo  | TRUTH.TOTAL | TP    | FN   | Recall | Precision | F1     |
|-------|-------------|-------|------|--------|-----------|--------|
| SNP   | 39.447      | 35.453| 3.994| 0,8988 | 0,9781    | 0,9367 |
| INDEL | 1.625       | 1.412 | 213  | 0,8689 | 0,6054    | 0,7136 |

Nota operativa: `Frac_NA` es alta (54 % en SNV, 79 % en INDEL) porque recoge las variantes
del VCF de consulta situadas fuera de `highconf_regions.bed`, es decir, zona no evaluable.
No son falsos positivos. La interpretación de estas métricas se desarrolla en la sección 5.2
de la memoria.

### Salidas generadas

    resultados/exp03/happy/
      seqc2_eval.summary.csv
      seqc2_eval.extended.csv
      seqc2_eval.vcf.gz (+ .tbi)
      seqc2_eval.roc.*.csv.gz
      seqc2_eval.runinfo.json
      truth.gt.vcf.gz (+ .tbi)   # truth set adaptado, reproducible

---

## Entrada NB-009 — Replanteamiento del módulo ML con etiqueta externa (exp07)

**Fecha:** 16 de agosto de 2026
**Máquina:** olmop-MS-7E26 (Ubuntu nativo)
**Estado:** Completado

### Motivo

Al revisar exp04 antes de cerrar la memoria aparecieron dos problemas en el planteamiento, no
en la ejecución:

1. **La etiqueta era predecible sin entrenar.** Una regla determinista ("tiene puntuación SIFT
   o PolyPhen, o es un indel") evaluada sobre el mismo test de exp04 da F1 = 0,9407, frente a
   0,9405 de Random Forest y 0,9406 de la regresión logística. Solo XGBoost la supera, y por
   nueve diezmilésimas. Script: `scripts/exp07_rule_baseline_exp04.py`.

2. **La ablación de exp05 no eliminaba lo que decía eliminar.** Al quitar `SIFT_missing` y
   `PolyPhen_missing` se mantiene el `SimpleImputer(strategy="median")`, de modo que el valor
   imputado sigue marcando qué filas carecían de puntuación:

   | Columna | Mediana (train) | Filas imputadas | Filas reales con ese valor | Pureza |
   |---|---|---|---|---|
   | `PolyPhen_score` | 0,600 | 1.386.273 | 348 | 99,97 % |
   | `SIFT_score` | 0,030 | 1.504.495 | 63.782 | 95,93 % |

   Un árbol que particione en esos valores reconstruye el indicador. La caída de AUC de 0,0010
   es, por tanto, compatible con que la señal siguiera disponible.

La causa común es la definición del problema: con la etiqueta derivada de
`Variant_Classification`, tener puntuación SIFT implica `label = 1` en el 100,0000 % de las
1.913.573 variantes que la tienen. Restringir la evaluación al subconjunto con puntuación,
que es la comprobación habitual, tampoco sirve aquí porque ese subconjunto contiene una sola
clase.

### Decisión

Conservar exp04 a exp06 como están y añadir un planteamiento nuevo con etiqueta externa, en
lugar de reajustar el anterior.

### Datos

ClinVar en GRCh38, `fileDate=2026-08-08`, descargado el 16 de agosto de 2026:

    mkdir -p datos-raw/clinvar
    curl -L -o datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz \
      https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    curl -L -o datos-raw/clinvar/clinvar_GRCh38_20260810.vcf.gz.tbi \
      https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

El fichero (185 MB) está en `.gitignore`. La versión queda registrada en
`resultados/exp07/clinvar_version.txt`.

### Cruce y filtrado

Coincidencia exacta `cromosoma:posición:ref:alt`, solo SNV: los indels se representan de forma
distinta en MAF y VCF y no se pueden emparejar sin normalizar. De los 3.425.534 SNV del MAF,
323.626 tienen registro en ClinVar (242.638 variantes distintas).

    python3 scripts/build_clinvar_dataset.py

| Paso | Variantes |
|---|---|
| Missense con clasificación P/LP o B/LB | 19.042 |
| Con `criteria_provided` en `CLNREVSTAT` | 18.462 |
| Con SIFT y PolyPhen disponibles | **15.734** |

Las 2.728 descartadas por falta de puntuación tienen una prevalencia de patogénicas del 22,3 %,
frente al 29,3 % de las retenidas. Se descartan en vez de imputarse: imputar es lo que
invalidó exp05.

### Tres decisiones de diseño

1. **Solo missense.** Sobre todas las clases, la consecuencia determina la etiqueta de ClinVar
   casi por completo (48.194 de 48.270 silenciosas son benignas; 6.474 de 6.527 nonsense son
   patogénicas). Entrenar sobre el conjunto completo reproduciría el problema de exp04. En
   missense conviven las dos clases: 13.396 benignas frente a 5.646 patogénicas.
2. **Agrupación por gen**, no por muestra. Aquí el eje de fuga es el gen: sin agrupar, el
   modelo memoriza que las variantes de TP53 tienden a ser patogénicas. El script comprueba
   con un `assert` que ningún gen aparece en las dos particiones.
3. **Sin imputación**, por lo dicho arriba.

### Resultados

    python3 scripts/exp07_train_pathogenicity.py

15.734 variantes, 4.605 patogénicas (29,3 %), 6.603 genes. Train 12.686 / 5.282 genes, test
3.048 / 1.321 genes.

| Modelo | AUC-ROC | PR-AUC | Precisión | Recall | F1 |
|---|---|---|---|---|---|
| SIFT (línea base) | 0,9132 | 0,7380 | 0,6547 | 0,8866 | 0,7532 |
| PolyPhen (línea base) | 0,9233 | 0,8033 | 0,7167 | 0,7847 | 0,7492 |
| SIFT < 0,05 y PolyPhen > 0,85 | 0,9323 | 0,8273 | 0,7728 | 0,7442 | 0,7583 |
| Regresión logística | 0,9330 | 0,8197 | 0,6763 | 0,8947 | 0,7703 |
| Random Forest | 0,9406 | 0,8427 | 0,7090 | 0,8854 | 0,7874 |
| XGBoost | 0,9412 | 0,8453 | 0,7699 | 0,7940 | 0,7818 |

La ganancia sobre la mejor línea base es de 0,0089 en AUC-ROC y 0,029 en F1. Es pequeña, pero
se mantiene en las tres métricas y coincide con la validación cruzada (RF: 0,9411 ± 0,0053).

Importancias del Random Forest: PolyPhen 0,493 y SIFT 0,447; las otras cinco suman 0,059. La
recurrencia en las 10.295 muestras del MAF aporta 0,009, previsiblemente porque la mayoría de
las variantes del conjunto aparece en una o dos muestras.

Comprobación sobre las variantes que ClinVar clasifica además como oncogénicas: solo 33 caen
en genes del conjunto de test; probabilidad media 0,843 y recall 0,879. Sin negativos, no es
una validación.

### Limitaciones anotadas para la memoria

- `CLNSIG` es significado clínico germinal, no oncogenicidad somática. Entre los genes con más
  patogénicas aparecen SCN1A, FBN1 y COL4A5, que son de enfermedad mendeliana.
- Los criterios ACMG admiten evidencia computacional (PP3/BP4), así que SIFT y PolyPhen pueden
  haber intervenido en el etiquetado. Las métricas son una cota superior.
- Solo SNV, una sola base de datos, y la comprobación de oncogenicidad se apoya en 33 casos.

### Salidas generadas

    resultados/exp07/
      clinvar_missense_dataset.tsv    # conjunto de modelado (15.734 filas)
      dataset_metadata.json
      exp07_results.json              # métricas, CV, importancias, comprobación ONC
      exp04_rule_baseline.json        # regla determinista sobre el test de exp04
      clinvar_version.txt
      test_predictions.csv            # etiqueta y puntuacion de cada modelo en test
      build_log.txt / train_log.txt / rule_baseline_log.txt
