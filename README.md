# Pipeline de análisis de variantes genéticas

TFM — Máster Universitario en Bioinformática, UAX

## Descripción

Pipeline reproducible en Nextflow DSL2 para el análisis de variantes genéticas (SNVs e indels) a partir de datos de secuenciación de exoma completo (WES), con clasificación patogénica automática mediante Random Forest y XGBoost entrenados sobre ClinVar.

## Estructura del repositorio
pipeline-variantes/
├── config/              # Perfiles de configuración (HPC, cloud)
├── datos-intermedios/   # BAMs, VCFs intermedios (en .gitignore)
├── datos-raw/           # FASTQs y referencias (en .gitignore)
│   └── tests/           # Datos reducidos para pruebas
├── docs/                # Documentación técnica
├── entorno/             # environment.yml y Dockerfile
├── resultados/
│   ├── exp01/           # Validación NA12878 (GIAB)
│   └── exp02/           # Clasificación ML sobre ClinVar
├── scripts/             # Scripts Python/R y notebooks
├── main.nf              # Pipeline principal (pendiente)
└── nextflow.config      # Configuración principal
## Requisitos

- Nextflow >= 25.10.4
- Docker Desktop (con integración WSL2 activada)
- Conda / Miniconda con canales bioconda y conda-forge

## Instalación del entorno

```bash
conda env create -f entorno/environment.yml
conda activate tfm-variantes
```

## Ejecución

```bash
# Prueba local (datos de test)
nextflow run main.nf -profile local --input datos-raw/tests/

# Cluster HPC (SLURM)
nextflow run main.nf -profile hpc --input datos-raw/
```

## Experimentos

| ID | Descripción | Estado |
|----|-------------|--------|
| exp01 | Validación pipeline sobre NA12878 chr20 (GIAB) | Pendiente |
| exp02 | Entrenamiento y evaluación modelos ML (ClinVar) | Pendiente |

## Autor

Juan Carlos Olmo Picon — jcarlosolmopicon@gmail.com  
Tutor: Beatriz Magan Pinto  
UAX, 2026
