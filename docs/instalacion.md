# Instalación y despliegue en otro equipo

Guía para clonar el repositorio y ejecutar el pipeline en una máquina distinta de la que
produjo los resultados publicados. Todas las rutas de `nextflow.config` son relativas al
directorio del proyecto (`${projectDir}`) o al `HOME` del usuario, de modo que no hay nada
que editar tras el clonado: basta con colocar los datos donde se indica.

---

## 1. Requisitos

| Componente | Versión | Notas |
|---|---|---|
| Linux x86-64 | — | Los resultados publicados se obtuvieron sobre Ubuntu 26.04 LTS |
| Docker | ≥ 24 | Motor de contenedores del perfil `local` |
| Java | ≥ 17 | Requisito de Nextflow |
| Nextflow | ≥ 24.0 | Los resultados publicados se obtuvieron con 25.10.4 |
| Miniconda | — | Para el módulo de clasificación (perfil `conda_local`) |

Espacio en disco orientativo para el par completo HCC1395/HCC1395BL: **~250 GB** de datos de
partida (FASTQ, referencia, MAF) y **~350 GB** de intermedios en
`datos-intermedios/nextflow-work/`. El modo de prueba (`-profile test`) trabaja con unos
pocos MB.

```bash
# Nextflow
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/

# Comprobación
nextflow -version && docker --version && java -version
```

## 2. Clonado y estructura esperada

```bash
git clone https://github.com/jcarlosolmopicon-ops/pipeline-variantes-tfm.git
cd pipeline-variantes-tfm
```

Los datos pesados no se versionan. El pipeline los espera bajo `datos-raw/`:

```
datos-raw/
├── HCC1395/            # FASTQ del par tumor-normal + truth_set/
├── TCGA-LUAD/          # MAF pan-cáncer MC3 (GRCh38)
├── clinvar/            # VCF de ClinVar (solo para exp07)
├── referencia/         # GRCh38.fa (+ índices), dbsnp138.vcf.gz
└── tests/              # FASTQ reducidos — SÍ versionados
```

La caché offline de VEP se instala en `~/.vep` (parámetro `params.vep_cache`).

## 3. Descarga de los datos

Todos los comandos de descarga, con sumas de verificación y notas de resolución de
problemas, están en [`descarga_datos.md`](descarga_datos.md):

| Sección | Datos |
|---|---|
| HCC1395 / HCC1395BL | FASTQ del SRA y truth set de SEQC2 |
| MAF pan-cáncer del TCGA | Fichero del GDC y liftover a GRCh38 |
| Genoma de referencia GRCh38 | FASTA e índices |
| Caché de VEP GRCh38 v110 | `vep_install … --CACHEDIR ~/.vep` |
| dbSNP138 | Sitios conocidos para BQSR |
| ClinVar | Etiquetas de patogenicidad (solo exp07) |

## 4. Prueba de humo (sin datos pesados)

Con solo el clon, los FASTQ de `datos-raw/tests/` permiten comprobar la instalación:

```bash
nextflow run main.nf -profile test --step qc
```

Debe terminar con FastQC y MultiQC en `resultados/test/`. Los pasos posteriores del perfil
`test` (`--step all`) ya requieren el genoma de referencia, dbSNP y la caché de VEP.

## 5. Ejecución completa

```bash
nextflow run main.nf -profile local \
  --input "datos-raw/HCC1395/*_{1,2}.fastq.gz" \
  --outdir resultados/exp03 \
  --step all \
  -resume
```

`--step` admite `qc`, `trim`, `align`, `call`, `annotate`, `classify` o `all`. El paso
`classify` aplica los modelos serializados de `resultados/exp04/` al VCF anotado del
`--outdir` indicado y publica `classify/hcc1395_annotated_with_ml.tsv` y su resumen JSON:

```bash
# Solo la clasificación, sobre un VCF ya anotado
nextflow run main.nf -profile conda_local -c config/exp08.config \
  --step classify --outdir resultados/exp03

# Equivalente como workflow independiente (experimento exp08)
nextflow run exp08.nf -profile conda_local -c config/exp08.config \
  --vcf resultados/exp03/vep/annotated.vcf.gz \
  --modelos resultados/exp04 \
  --outdir resultados/exp08
```

La inferencia es **determinista**: dos ejecuciones producen ficheros bit-idénticos (la
predicción se fuerza a un solo hilo al cargar los modelos; véase NB-010 del lab notebook).

## 6. Entornos conda

No hay que crear entornos a mano para ejecutar el pipeline: el proceso `CLASSIFY` declara
`entorno/environment_ubuntu_2026-05.yml` y Nextflow lo construye y cachea en
`~/.conda-nf-cache` la primera vez. Para trabajar interactivamente con los scripts del
módulo de clasificación:

```bash
conda env create -f entorno/environment_ubuntu_2026-05.yml
conda activate tfm-variantes
```

## 7. Problemas conocidos

**Recolección de métricas de Nextflow y coreutils de uutils.** En sistemas donde
`/usr/bin/date` es el de uutils-coreutils ≥ 0.8 (Ubuntu 25.10+ actualizado), `date +%s%3N`
devuelve nanosegundos completos en lugar de milisegundos. La función `nxf_date` del wrapper
de Nextflow 25.10.4 no lo tolera y la recolección de métricas de `trace`/`report`/`timeline`
**aborta la tarea después de completarse el cómputo** (error
`.command.run: … variable sin asignar`). Los resultados del proceso son correctos; solo
falla la contabilidad. Solución: cargar `config/exp08.config`, que desactiva esos informes:

```bash
nextflow run main.nf -profile local -c config/exp08.config …
```

Comprobación rápida del sistema: `date +%s%3N` debe devolver **13 dígitos**; si devuelve 19,
el sistema está afectado.

**Plataforma de las imágenes Docker.** El perfil `local` fija `--platform linux/amd64`; en
hosts ARM las imágenes se ejecutan bajo emulación, con la penalización correspondiente.
