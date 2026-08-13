# Estrategia de almacenamiento
## Actualizado: 25 de abril de 2026 (revisado en agosto de 2026)

> **Este documento describe la estrategia inicial, ya superada.** Durante la ejecucion
> completa del pipeline el proyecto migro de WSL2 a **Ubuntu nativo en arranque dual**
> (maquina olmop-MS-7E26), por problemas de memoria y de rendimiento de E/S sobre el
> disco externo montado via drvfs. Ver NB-007 del lab notebook.
>
> A partir de esa migracion **nada de lo que sigue aplica**: no hay montaje de K:, no
> hay separacion entre temporales en C: y datos en el disco externo, y no hay VHDX que
> compactar. La estrategia vigente se resume al final del documento.

---

## Arquitectura de discos

| Disco | Montaje WSL2 | Libre | Uso |
|-------|-------------|-------|-----|
| C: (interno) | ~ (/dev/sdd) | ~866 GB | Temporales fasterq-dump, conda, codigo, Nextflow workDir |
| D: (interno) | /mnt/d | ~319 GB | Reserva / overflow |
| K: (externo USB) | /mnt/k | ~924 GB | Datos permanentes del TFM |

## Regla de oro

- Temporales (fasterq-dump --temp) → siempre en C: via ~/tmp-fasterq/
- Datos finales (FASTQs .gz, BAMs, VCFs, referencias) → siempre en K:
- Codigo y entorno → siempre en C: (~/pipeline-variantes/, ~/miniconda3/)
- Nunca usar /mnt/c/ para datos del proyecto (lento, inodos no POSIX)

## Montar disco K: en WSL2

    sudo mkdir -p /mnt/k
    sudo mount -t drvfs K: /mnt/k

Ejecutar en cada inicio de WSL2 — el montaje no es persistente entre sesiones.

## Comandos de inicio de sesion

    sudo mount -t drvfs K: /mnt/k
    conda activate tfm-variantes
    cd ~/pipeline-variantes
    git status
    df -h ~ && df -h /mnt/k

## Estado de datos en K: (25/04/2026)

| Fichero | Tamano | Estado |
|---------|--------|--------|
| SRR7890824_1.fastq.gz (tumor R1) | 45 GB | Verificado con gzip -t |
| SRR7890824_2.fastq.gz (tumor R2) | 50 GB | Verificado con gzip -t |
| SRR7890827.sra (normal) | 70 GB | Conversion en curso |
| MAF pan-cancer GRCh38 (descargado como TCGA-LUAD) | ~718 MB | Liftover completado |
| hg19ToHg38.over.chain.gz | — | Descargado |

## Pendiente de descargar

| Fichero | Tamano estimado | Prioridad |
|---------|----------------|-----------|
| GRCh38 referencia (hg38.fa.gz) | ~3 GB | Alta |
| Indices BWA | ~5 GB | Alta |
| VEP cache GRCh38 v110 | ~15 GB | Alta |
| dbSNP138 (para BQSR, bundle GATK) | ~20 GB | Alta |
| gnomAD v4 exomas | ~30 GB | No descargado (la cache VEP ya incluye gnomAD) |
| CADD v1.7 | 0-80 GB | Descartado, no instalado |
| SEQC2 truth set VCF | ~50 MB | Alta |

## Balance de almacenamiento estimado al completar el TFM

| Concepto | Tamano | Disco |
|----------|--------|-------|
| FASTQs tumor comprimidos | 95 GB | K: |
| FASTQs normal comprimidos | ~95 GB | K: |
| Referencia + indices BWA | ~8 GB | K: |
| VEP cache + dbSNP | ~35 GB | K: |
| gnomAD exomas | ~30 GB | K: |
| BAMs alineados tumor+normal | ~100 GB | K: (borrar tras VCF) |
| VCFs intermedios y finales | ~5 GB | K: |
| MAF pan-cancer MC3 del TCGA | ~1 GB | K: |
| Total estimado K: | ~424 GB | K: libre actual: 924 GB |

## Nota sobre compactacion VHDX

El disco virtual ext4.vhdx crece dinamicamente pero no se compacta solo.
Compactar con diskpart si C: baja de 200 GB libres. Ver NB-006 del lab
notebook para el procedimiento completo.

    # Verificar espacio antes de operaciones grandes
    df -h ~ && df -h /mnt/k

---

## Estrategia vigente (tras la migracion a Ubuntu nativo)

Todo el proyecto reside en el sistema de ficheros nativo (ext4 sobre NVMe), sin capas
de traduccion intermedias. Esto elimina de raiz los tres problemas que motivaron la
migracion: los limites de memoria de la maquina virtual de WSL2, la lentitud del acceso
via drvfs y el crecimiento descontrolado del disco virtual.

| Contenido | Ubicacion |
|-----------|-----------|
| Codigo y entorno | ~/pipeline-variantes/, ~/miniconda3/ |
| Datos de partida (FASTQ, referencias, MAF) | datos-raw/ |
| Intermedios de Nextflow (workDir) | datos-intermedios/ |
| Resultados | resultados/expNN/ |

### Que se versiona y que no

En Git se versionan el codigo, la configuracion, la documentacion y los ficheros de
resultados ligeros (metricas, informes, modelos serializados). Quedan excluidos por
.gitignore los ficheros pesados: FASTQ, BAM alineados y recalibrados, referencias e
indices, y el MAF de entrenamiento. El pipeline los regenera de forma determinista a
partir de los datos de partida, siguiendo la practica habitual en pipelines genomicos.

### Comandos de inicio de sesion

    conda activate tfm-variantes
    cd ~/pipeline-variantes
    git status
    df -h .
