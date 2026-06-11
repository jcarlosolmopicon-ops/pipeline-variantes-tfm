# Pipeline de análisis de variantes somáticas

TFM — Máster Universitario en Bioinformática, UAX

## Descripción

Pipeline reproducible en Nextflow DSL2 para la detección de variantes somáticas a partir de datos de secuenciación de exoma completo (WES) tumor-normal, con un módulo de clasificación patogénica mediante Random Forest y XGBoost entrenado sobre TCGA-LUAD (en desarrollo).

## Estructura del repositorio

```
pipeline-variantes/
├── config/                    # Perfiles de configuración (local/Docker, HPC/SLURM)
├── modules/                   # Módulos Nextflow DSL2 (fastqc, bwa, bqsr, mutect2, vep, etc.)
├── datos-intermedios/         # workDir de Nextflow (en .gitignore)
├── datos-raw/                 # FASTQs y referencias (en .gitignore)
│   └── HCC1395/truth_set/     # Truth set SEQC2 (SNV + INDEL HighConf)
├── docs/                       # Lab notebook, ADRs y documentación técnica
├── entorno/                    # environment.yml (Conda)
├── resultados/
│   └── exp03_run_completa/     # Pipeline completo sobre HCC1395 (FastQC→VEP→hap.py)
├── main.nf                     # Pipeline principal
└── nextflow.config             # Configuración principal (perfiles local/conda_local)
```
## Requisitos

- Nextflow >= 25.10.4
- Docker (perfil `local`, activo en producción)
- Conda / Miniconda con canales bioconda y conda-forge (perfil `conda_local`, alternativo)

## Instalación del entorno

```bash
conda env create -f entorno/environment.yml
conda activate tfm-variantes
```

## Ejecución

```bash
nextflow run main.nf -profile local \
  --input "datos-raw/HCC1395/*_{1,2}.fastq.gz" \
  --outdir resultados/exp03_run_completa \
  --step all \
  -resume
```

## Experimentos

| ID | Descripción | Estado |
|----|-------------|--------|
| exp03_run_completa | Pipeline completo sobre par tumor-normal SEQC2/HCC1395 (FastQC, BWA, BQSR, MuTect2, VEP) + validación hap.py contra truth set SEQC2 | Completado |
| Módulo ML (TCGA-LUAD) | Clasificación de patogenicidad con Random Forest/XGBoost sobre variantes de TCGA-LUAD (liftover GRCh38 completado, 99.99%) | Pendiente |

### Resultados de validación (hap.py vs. truth set SEQC2/HCC1395)

| Tipo | Recall | Precisión | F1-score |
|------|--------|-----------|----------|
| SNV | 0.899 | 0.978 | 0.937 |
| INDEL | 0.869 | 0.605 | 0.714 |

## Autor

Juan Carlos Olmo Picón — jcarlosolmopicon@gmail.com  
Tutor: Beatriz Magán Pinto  
UAX, 2026
