# Pipeline de análisis de variantes somáticas

TFM — Máster Universitario en Bioinformática, UAX

## Descripción

Pipeline reproducible en Nextflow DSL2 para la detección de variantes somáticas a partir de datos de secuenciación de exoma completo (WES) tumor-normal, con un módulo de clasificación de relevancia funcional mediante Random Forest y XGBoost entrenado sobre el MAF somático pan-cáncer MC3 del TCGA.

El flujo va de las lecturas FASTQ crudas hasta un fichero TSV en el que cada variante detectada lleva asociada una probabilidad de relevancia funcional:

```
FASTQ → QC → alineamiento → preprocesamiento BAM → MuTect2 → VEP ─┬→ validación hap.py
                                                                  └→ clasificación ML
```

## Estructura del repositorio

```
pipeline-variantes/
├── config/                    # Perfiles de configuración (local/Docker, HPC/SLURM)
├── modules/                   # Módulos Nextflow DSL2 (fastqc, bwa, bqsr, mutect2, vep, etc.)
├── scripts/                   # Módulo ML: entrenamiento, aplicación y experimentos exp05/exp06
├── datos-intermedios/         # workDir de Nextflow (en .gitignore)
├── datos-raw/                 # FASTQs y referencias (en .gitignore)
│   └── HCC1395/truth_set/     # Truth set SEQC2 (SNV + INDEL, HighConf y MedConf)
├── docs/                      # Lab notebook (NB-001 a NB-006) y documentación técnica
├── entorno/                   # environment.yml (Conda)
├── resultados/
│   ├── exp03/                 # Pipeline completo sobre HCC1395 (FastQC→VEP→hap.py)
│   ├── exp04/                 # Módulo ML: entrenamiento, evaluación y aplicación
│   ├── exp05/                 # Experimento de ablación
│   ├── exp06/                 # Línea base de regresión logística
│   ├── exp01/                 # ⚠️ Pruebas tempranas sobre subconjuntos (ver nota)
│   └── exp02/                 # ⚠️ Pruebas tempranas sobre subconjuntos (ver nota)
├── main.nf                    # Pipeline principal
└── nextflow.config            # Configuración principal (perfiles local/conda_local)
```

> **Nota sobre `exp01/` y `exp02/`.** Corresponden a ejecuciones de prueba tempranas sobre
> subconjuntos reducidos de datos y no guardan relación con los experimentos descritos en la
> memoria. Se conservan solo por trazabilidad histórica.

> **Nomenclatura.** En el texto de la memoria los dos experimentos principales se denominan
> *exp01* (validación sobre HCC1395) y *exp02* (módulo ML). Por la evolución incremental del
> proyecto, en este repositorio corresponden a `exp03/` y `exp04/` respectivamente.

## Requisitos

- Nextflow >= 25.10.4
- Docker (perfil `local`, activo en producción)
- Conda / Miniconda con canales bioconda y conda-forge (perfil `conda_local`, alternativo)

El perfil `local` con Docker es el **entorno canónico de ejecución**: las versiones efectivamente
empleadas en cada etapa del pipeline son las declaradas en las imágenes de los módulos.
`environment.yml` se mantiene como apoyo para el desarrollo interactivo y para el módulo de
aprendizaje automático.

## Instalación del entorno

```bash
conda env create -f entorno/environment.yml
conda activate tfm-variantes
```

## Ejecución

```bash
nextflow run main.nf -profile local \
  --input "datos-raw/HCC1395/*_{1,2}.fastq.gz" \
  --outdir resultados/exp03 \
  --step all \
  -resume
```

El parámetro `--step` admite `qc`, `trim`, `align`, `call`, `annotate` o `all`, lo que permite
ejecutar etapas concretas durante el desarrollo.

### Módulo de aprendizaje automático

```bash
python3 scripts/train_ml_module.py              # exp04: entrenamiento y evaluación
python3 scripts/train_final_model.py            # modelo final sobre el conjunto completo
python3 scripts/apply_model_hcc1395.py          # aplicación sobre las variantes de HCC1395
python3 scripts/exp05_ablation_missing_flags.py # exp05: ablación de indicadores estructurales
python3 scripts/exp06_baseline_logreg.py        # exp06: línea base de regresión logística
```

## Experimentos

| ID | Descripción | Estado |
|----|-------------|--------|
| exp03 | Pipeline completo sobre par tumor-normal SEQC2/HCC1395 (FastQC, BWA, BQSR, MuTect2, VEP) + validación hap.py contra truth set SEQC2 | Completado |
| exp04 | Módulo ML: entrenamiento y evaluación de Random Forest y XGBoost sobre el MAF pan-cáncer MC3, y aplicación sobre las variantes de HCC1395 | Completado |
| exp05 | Ablación: cuantifica cuánto del rendimiento depende de los indicadores `SIFT_missing` y `PolyPhen_missing` | Completado |
| exp06 | Línea base: regresión logística sobre las mismas características y particiones, para contextualizar el AUC de los modelos de conjunto | Completado |

### Resultados de validación (hap.py vs. truth set SEQC2/HCC1395)

| Tipo | Recall | Precisión | F1-score |
|------|--------|-----------|----------|
| SNV | 0.899 | 0.978 | 0.937 |
| INDEL | 0.869 | 0.605 | 0.714 |

La menor precisión en INDELs es coherente con la ausencia de un Panel of Normals dedicado y con
las dificultades conocidas de representación de inserciones y deleciones en el alineamiento.

### Resultados del módulo ML (test held-out, agrupado por muestra)

| Modelo | AUC-ROC | Precisión | Recall | F1-score |
|--------|---------|-----------|--------|----------|
| Random Forest | 0.9595 | 0.9887 | 0.8968 | 0.9405 |
| XGBoost | 0.9596 | 0.9852 | 0.9017 | 0.9416 |
| Regresión logística (línea base) | 0.9565 | 0.9780 | 0.9059 | 0.9406 |

La línea base lineal queda a tres milésimas de AUC de los modelos de conjunto y los supera en
recall, lo que indica que el problema está determinado en gran medida por las características de
anotación funcional disponibles y no por la complejidad del clasificador. El experimento de
ablación (exp05) apunta en la misma dirección: eliminar `SIFT_missing` y `PolyPhen_missing`
reduce el AUC en apenas 0.0010 (RF) y 0.0014 (XGBoost).

Aplicados sobre las 126.014 variantes PASS del VCF anotado de HCC1395, Random Forest clasifica
7.232 (5.7 %) como funcionalmente relevantes y XGBoost 11.930 (9.5 %), con umbral 0.5.

## Datos

**Validación del pipeline.** Par tumor-normal SEQC2/HCC1395 (línea celular de cáncer de mama
triple negativo): `SRR7890824` (tumor, ~100x) y `SRR7890827` (normal, ~80x), WES en Illumina
HiSeq X. BioProject `SRP162370`. Truth set con 39.560 SNVs y 1.922 INDELs.

**Entrenamiento del clasificador.** MAF somático **pan-cáncer MC3 del TCGA** (Ellrott et al.,
2018): 3.598.760 variantes de **10.295 muestras tumorales** de **33 tipos de cáncer**, en
coordenadas GRCh38 tras liftover con CrossMap 0.7.2.

> ⚠️ **Denominación del fichero.** El MAF se descargó del GDC bajo el nombre del proyecto
> TCGA-LUAD, pero su contenido real corresponde al conjunto pan-cáncer MC3: el proyecto LUAD
> aporta solo el 6.1 % de las variantes (219.815 de 509 muestras). La incidencia se detectó
> durante la verificación final de los datos y se documenta en la sección 5.3.1 de la memoria.
> No invalida los resultados, pero obliga a interpretar el clasificador como un **modelo
> pan-cáncer**, no específico de adenocarcinoma de pulmón.

Los ficheros pesados (FASTQ, BAM alineados y recalibrados, referencias y MAF) están excluidos
mediante `.gitignore`, siguiendo la práctica habitual en pipelines genómicos. El pipeline los
regenera de forma determinista a partir de los datos de partida.

## Limitaciones conocidas

- No se empleó un **Panel of Normals** dedicado, lo que explica parte de la baja precisión en INDELs.
- La anotación no incluye **CADD ni SpliceAI** (requieren ficheros de plugin externos no
  instalados) ni las frecuencias de **gnomAD**, **COSMIC** o **ClinVar**, que sí están en la caché
  offline pero no se solicitaron en la invocación de VEP.
- La variable objetivo del clasificador es un **proxy de impacto funcional** derivado de
  `Variant_Classification`, no patogenicidad clínica en sentido estricto.
- La separación entrenamiento/prueba se agrupó por muestra tumoral y no por paciente: de las
  10.295 muestras, 71 pacientes aportan más de una, de modo que la fuga potencial afecta como
  máximo al 1.4 % de los casos.

## Documentación

| Fichero | Contenido |
|---------|-----------|
| `docs/lab_notebook.md` | Registro técnico de sesiones (entradas NB-001 a NB-006) |
| `docs/descarga_datos.md` | Procedimiento de descarga de los datos de partida y las referencias |
| `docs/almacenamiento_datos.md` | Estrategia de almacenamiento de la fase inicial del proyecto |
| `docs/liftover_GRCh37_to_GRCh38.md` | Conversión de coordenadas del MAF con CrossMap |

> **Entornos de ejecución.** El proyecto se desarrolló en dos entornos sucesivos. La fase
> inicial (entradas NB-001 a NB-006) se llevó a cabo sobre Windows 11 con WSL2 y un disco
> externo montado vía drvfs. Durante la ejecución completa del pipeline se migró a **Ubuntu
> nativo en arranque dual** por problemas de memoria y de rendimiento de E/S, migración que
> se documenta en la entrada NB-007 del lab notebook. Todos los resultados descritos en la
> memoria (exp03 a exp06) se obtuvieron en ese segundo entorno: 12 CPUs, 28 GB de RAM y
> almacenamiento NVMe con ext4 nativo.

## Autor

Juan Carlos Olmo Picón — jcarlosolmopicon@gmail.com  
Tutora: Beatriz Magán Pinto  
UAX, 2026
