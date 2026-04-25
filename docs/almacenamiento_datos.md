# Estrategia de almacenamiento
## Actualizado: 25 de abril de 2026

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
| TCGA-LUAD.mutect2.GRCh38.maf.gz | ~718 MB | Liftover completado |
| hg19ToHg38.over.chain.gz | — | Descargado |

## Pendiente de descargar

| Fichero | Tamano estimado | Prioridad |
|---------|----------------|-----------|
| GRCh38 referencia (hg38.fa.gz) | ~3 GB | Alta |
| Indices BWA-MEM2 | ~60 GB | Alta |
| VEP cache GRCh38 v110 | ~15 GB | Alta |
| dbSNP b156 (para BQSR) | ~20 GB | Alta |
| gnomAD v4 exomas | ~30 GB | Media |
| CADD v1.7 (o usar API) | 0-80 GB | Media |
| SEQC2 truth set VCF | ~50 MB | Alta |

## Balance de almacenamiento estimado al completar el TFM

| Concepto | Tamano | Disco |
|----------|--------|-------|
| FASTQs tumor comprimidos | 95 GB | K: |
| FASTQs normal comprimidos | ~95 GB | K: |
| Referencia + indices BWA-MEM2 | ~63 GB | K: |
| VEP cache + dbSNP | ~35 GB | K: |
| gnomAD exomas | ~30 GB | K: |
| BAMs alineados tumor+normal | ~100 GB | K: (borrar tras VCF) |
| VCFs intermedios y finales | ~5 GB | K: |
| TCGA-LUAD MAF | ~1 GB | K: |
| Total estimado K: | ~424 GB | K: libre actual: 924 GB |

## Nota sobre compactacion VHDX

El disco virtual ext4.vhdx crece dinamicamente pero no se compacta solo.
Compactar con diskpart si C: baja de 200 GB libres. Ver NB-006 del lab
notebook para el procedimiento completo.

    # Verificar espacio antes de operaciones grandes
    df -h ~ && df -h /mnt/k
