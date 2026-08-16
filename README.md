# Pipeline de análisis de variantes somáticas

TFM — Máster Universitario en Bioinformática, UAX

## Descripción

Pipeline reproducible en Nextflow DSL2 para la detección de variantes somáticas a partir de datos de secuenciación de exoma completo (WES) tumor-normal, con un módulo de clasificación mediante Random Forest y XGBoost entrenado sobre el MAF somático pan-cáncer MC3 del TCGA.

El módulo de clasificación tiene dos planteamientos. El primero (exp04) predice relevancia funcional a partir de una etiqueta derivada de la propia anotación; su evaluación crítica (exp05 a exp07) muestra que una regla determinista sin entrenamiento iguala a los modelos, de modo que el rendimiento procede de la definición de las clases. El segundo (exp07) sustituye esa etiqueta por las clasificaciones de patogenicidad de ClinVar y se restringe a variantes missense, donde la pregunta no está resuelta de antemano.

El flujo va de las lecturas FASTQ crudas hasta un fichero TSV en el que cada variante detectada lleva asociada una probabilidad:

```
FASTQ → QC → alineamiento → preprocesamiento BAM → MuTect2 → VEP ─┬→ validación hap.py
                                                                  └→ clasificación ML
```

## Estructura del repositorio

```
pipeline-variantes/
├── config/                    # Perfiles de configuración (local/Docker, HPC/SLURM)
├── modules/                   # Módulos Nextflow DSL2 (fastqc, bwa, bqsr, mutect2, vep, etc.)
├── scripts/                   # Módulo ML: entrenamiento, aplicación y experimentos exp05 a exp07
├── datos-intermedios/         # workDir de Nextflow (en .gitignore)
├── datos-raw/                 # FASTQs, referencias y ClinVar (en .gitignore)
│   ├── HCC1395/truth_set/     # Truth set SEQC2 (SNV + INDEL, HighConf y MedConf)
│   └── clinvar/               # VCF de ClinVar GRCh38 (etiquetas de patogenicidad)
├── docs/                      # Lab notebook (NB-001 a NB-009) y documentación técnica
├── entorno/                   # environment.yml (Conda)
├── resultados/
│   ├── exp03/                 # Pipeline completo sobre HCC1395 (FastQC→VEP→hap.py)
│   ├── exp04/                 # Módulo ML: entrenamiento, evaluación y aplicación
│   ├── exp05/                 # Experimento de ablación
│   ├── exp06/                 # Línea base de regresión logística
│   └── exp07/                 # Clasificación de patogenicidad con etiqueta de ClinVar
├── main.nf                    # Pipeline principal
└── nextflow.config            # Configuración principal (perfiles local/conda_local)
```

> **Nomenclatura de los experimentos.** En el texto de la memoria los dos experimentos
> principales se denominan *exp01* (validación sobre HCC1395) y *exp02* (módulo ML). Por la
> evolución incremental del proyecto, en este repositorio corresponden a `exp03/` y `exp04/`.
> Las carpetas locales `exp01/` y `exp02/`, que contienen ejecuciones de prueba tempranas
> sobre subconjuntos reducidos, no se versionan.

## Requisitos

- Nextflow >= 24.0 (los resultados publicados se obtuvieron con 25.10.4)
- Docker (perfil `local`, activo en producción)
- Conda / Miniconda con canales bioconda y conda-forge (perfil `conda_local`, alternativo)

El perfil `local` con Docker es el **entorno canónico de ejecución**: las versiones efectivamente
empleadas en cada etapa del pipeline son las declaradas en las imágenes de los módulos.
`environment.yml` se mantiene como apoyo para el desarrollo interactivo y para el módulo de
aprendizaje automático.

## Instalación del entorno

Para reproducir los resultados publicados, usar la exportación del entorno con el que se
obtuvieron (Ubuntu 26.04, mayo de 2026):

```bash
conda env create -f entorno/environment_ubuntu_2026-05.yml
conda activate tfm-variantes
```

`entorno/environment.yml` declara las dependencias sin fijar versiones y sirve para levantar un
entorno de trabajo actualizado; `entorno/environment_from_history.yml` es la exportación de la
fase inicial sobre WSL2 y se conserva por trazabilidad histórica.

Las herramientas del pipeline no proceden de ninguno de estos ficheros: se ejecutan en
contenedores Docker con etiqueta fija declarada en cada `modules/*.nf`.

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
python3 scripts/exp07_rule_baseline_exp04.py    # regla determinista sobre el test de exp04
python3 scripts/build_clinvar_dataset.py        # exp07: cruce del MAF con ClinVar
python3 scripts/exp07_train_pathogenicity.py    # exp07: entrenamiento y líneas base
```

Requieren el entorno `tfm-variantes` (`conda activate tfm-variantes`), que es donde están
scikit-learn y XGBoost.

## Experimentos

| ID | Descripción |
|----|-------------|
| exp03 | Pipeline completo sobre par tumor-normal SEQC2/HCC1395 (FastQC, BWA, BQSR, MuTect2, VEP) + validación hap.py contra truth set SEQC2 |
| exp04 | Módulo ML: entrenamiento y evaluación de Random Forest y XGBoost sobre el MAF pan-cáncer MC3, y aplicación sobre las variantes de HCC1395 |
| exp05 | Ablación: retira los indicadores `SIFT_missing` y `PolyPhen_missing`. Ver la nota de más abajo sobre su interpretación |
| exp06 | Línea base: regresión logística sobre las mismas características y particiones, para contextualizar el AUC de los modelos de conjunto |
| exp07 | Clasificación de patogenicidad con etiqueta externa de ClinVar, restringida a variantes missense y con partición agrupada por gen. Incluye la regla determinista evaluada sobre el test de exp04 |

### Resultados de validación (hap.py vs. truth set SEQC2/HCC1395)

| Tipo | Recall | Precisión | F1-score |
|------|--------|-----------|----------|
| SNV | 0.899 | 0.978 | 0.937 |
| INDEL | 0.869 | 0.605 | 0.714 |

### Primer planteamiento — relevancia funcional (exp04, test held-out agrupado por muestra)

| Modelo | AUC-ROC | Precisión | Recall | F1-score |
|--------|---------|-----------|--------|----------|
| Random Forest | 0.9595 | 0.9887 | 0.8968 | 0.9405 |
| XGBoost | 0.9596 | 0.9852 | 0.9017 | 0.9416 |
| Regresión logística (línea base) | 0.9565 | 0.9780 | 0.9059 | 0.9406 |
| **Regla determinista (sin entrenamiento)** | — | 0.9781 | 0.9060 | **0.9407** |

La regla marca como relevante toda variante con puntuación SIFT o PolyPhen y todo indel. Iguala
a Random Forest y a la regresión logística sobre la misma partición, lo que sitúa el
rendimiento en la definición de las clases: tener puntuación SIFT implica la etiqueta positiva
en el 100.0000 % de las 1.913.573 variantes que la tienen.

La ablación de exp05 reduce el AUC en solo 0.0010 (RF) y 0.0013 (XGBoost), pero ese resultado
no acredita independencia de la señal estructural: al retirar los indicadores se mantiene la
imputación por la mediana, y el valor imputado sigue marcando qué filas carecían de puntuación
(pureza del 99.97 % en PolyPhen y del 95.93 % en SIFT).

Aplicados sobre las 126.014 variantes PASS del VCF anotado de HCC1395, Random Forest clasifica
7.232 (5.7 %) como funcionalmente relevantes y XGBoost 11.930 (9.5 %), con umbral 0.5.

### Segundo planteamiento — patogenicidad con etiqueta de ClinVar (exp07)

15.734 variantes missense de 6.603 genes, 29.3 % patogénicas. Partición y validación cruzada
agrupadas por gen; test de 3.048 variantes de 1.321 genes.

| Modelo | AUC-ROC | PR-AUC | Precisión | Recall | F1-score |
|--------|---------|--------|-----------|--------|----------|
| SIFT (línea base) | 0.9132 | 0.7380 | 0.6547 | 0.8866 | 0.7532 |
| PolyPhen (línea base) | 0.9233 | 0.8033 | 0.7167 | 0.7847 | 0.7492 |
| SIFT < 0.05 y PolyPhen > 0.85 | 0.9323 | 0.8273 | 0.7728 | 0.7442 | 0.7583 |
| Regresión logística | 0.9330 | 0.8197 | 0.6763 | 0.8947 | 0.7703 |
| Random Forest | 0.9406 | 0.8427 | 0.7090 | 0.8854 | 0.7874 |
| XGBoost | 0.9412 | 0.8453 | 0.7699 | 0.7940 | 0.7818 |

Sobre la mejor línea base, los modelos ganan 0.0089 de AUC-ROC y 0.029 de F1. Es poco, pero se
mantiene en las tres métricas y coincide con la validación cruzada (RF: 0.9411 ± 0.0053). El F1
absoluto baja de 0.94 a 0.79 respecto de exp04 porque el problema deja de estar resuelto por la
categoría de consecuencia.

## Datos

**Validación del pipeline.** Par tumor-normal SEQC2/HCC1395 (línea celular de cáncer de mama
triple negativo): `SRR7890824` (tumor, ~100x) y `SRR7890827` (normal, ~80x), WES en Illumina
HiSeq X. BioProject `SRP162370`. Truth set con 39.560 SNVs y 1.922 INDELs.

**Entrenamiento del clasificador.** MAF somático **pan-cáncer MC3 del TCGA** (Ellrott et al.,
2018): 3.598.760 variantes de **10.295 muestras tumorales** de **33 tipos de cáncer**, en
coordenadas GRCh38 tras liftover con CrossMap 0.7.2.

> **Denominación del fichero.** El MAF se descargó del GDC bajo el nombre del proyecto
> TCGA-LUAD, pero su contenido real corresponde al conjunto pan-cáncer MC3: el proyecto LUAD
> aporta solo el 6.1 % de las variantes (219.815 de 509 muestras). La incidencia se detectó
> durante la verificación final de los datos y se documenta en la sección 5.3.1 de la memoria.
> No invalida los resultados, pero obliga a interpretar el clasificador como un **modelo
> pan-cáncer**, no específico de adenocarcinoma de pulmón.

**Etiquetas de patogenicidad.** ClinVar en GRCh38 (`fileDate` 2026-08-08), de acceso libre y sin
registro. Se cruza con el MAF por coincidencia exacta de cromosoma, posición y alelos,
restringido a SNV. De los 3.425.534 SNV del MAF, 323.626 tienen registro en ClinVar; el
subconjunto missense con clasificación patogénica o benigna y criterios de revisión declarados
son 18.462 variantes, de las que 15.734 conservan SIFT y PolyPhen.

Los ficheros pesados (FASTQ, BAM alineados y recalibrados, referencias, MAF y el VCF de
ClinVar) están excluidos mediante `.gitignore`. Los detalles de descarga están en `docs/descarga_datos.md`.

## Limitaciones conocidas

- No se empleó un **Panel of Normals** dedicado, lo que explica parte de la baja precisión en INDELs.
- La anotación de HCC1395 no incluye **CADD ni SpliceAI** (requieren ficheros de plugin externos
  no instalados) ni las frecuencias de **gnomAD**, **COSMIC** o **ClinVar**, que sí están en la
  caché offline pero no se solicitaron en la invocación de VEP. Las etiquetas de exp07 no salen
  de esa caché, sino del VCF de ClinVar descargado aparte.
- La variable objetivo de exp04 es un **proxy de impacto funcional** derivado de
  `Variant_Classification`, no patogenicidad clínica en sentido estricto, y además es
  predecible sin modelo (ver los resultados de arriba).
- La separación entrenamiento/prueba de exp04 se agrupó por muestra tumoral y no por paciente:
  las 10.295 muestras corresponden a 10.224 pacientes y 71 pacientes aportan más de una
  muestra, de modo que la fuga potencial afecta como máximo al 1.4 % de las muestras.
- En exp07, `CLNSIG` expresa **significado clínico germinal**, no oncogenicidad somática: entre
  los genes con más variantes patogénicas aparecen SCN1A, FBN1 o COL4A5, de enfermedad
  mendeliana. El campo `ONC` de ClinVar, que sí es oncogenicidad, solo cruza con 544 variantes
  y ninguna benigna, insuficiente para entrenar.
- Los criterios ACMG admiten evidencia computacional (PP3, BP4), de modo que SIFT y PolyPhen
  pueden haber intervenido en el etiquetado de ClinVar. Las métricas de exp07 son una cota
  superior.
- El cruce con ClinVar se limita a SNV y descarta 2.728 variantes sin puntuación funcional,
  cuya prevalencia de patogénicas (22.3 %) difiere de la del conjunto retenido (29.3 %).

## Documentación

| Fichero | Contenido |
|---------|-----------|
| `docs/lab_notebook.md` | Registro técnico de sesiones (entradas NB-001 a NB-009) |
| `docs/descarga_datos.md` | Procedimiento de descarga de los datos de partida y las referencias |
| `docs/liftover_GRCh37_to_GRCh38.md` | Conversión de coordenadas del MAF con CrossMap |
| `datos-raw/README.md` | Descripción de los datos de partida, del truth set y de ClinVar |
| `entorno/environment_ubuntu_2026-05.yml` | Entorno conda con el que se obtuvieron los resultados |

El pipeline se ejecutó sobre Ubuntu 26.04 LTS con Nextflow 25.10.4 y Docker 29.1.3, en una
estación de trabajo con 12 CPUs, 28 GB de RAM y almacenamiento NVMe. La fase inicial del proyecto
se desarrolló sobre WSL2; la migración se documenta en NB-007.

## Autor

Juan Carlos Olmo Picón   
Tutora: Beatriz Magán Pinto  
UAX, 2026
