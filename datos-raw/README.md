# Datos crudos del proyecto

## Exp01 — Validación del pipeline: SEQC2 / HCC1395

| Muestra | Tipo | Accesión SRA | Cobertura | Tecnología |
|---------|------|-------------|-----------|-----------|
| HCC1395 | Tumor (TNBC) | SRR7890824 | ~100x WES | Illumina HiSeq X |
| HCC1395BL | Normal matched | SRR7890827 | ~80x WES | Illumina HiSeq X |

- **BioProject:** SRP162370
- **Consorcio:** SEQC2 Somatic Mutation Working Group
- **Publicación:** Zhao et al., Sci Data 2021 / Fang et al., Nat Biotechnol 2021
- **Truth set:** https://sites.google.com/view/seqc2
- **Referencia:** GRCh38
- **Tipo de variantes:** SNVs somáticas (~39.536 HighConf) + INDELs (~2.020)

## Exp02 — Módulo ML: TCGA-LUAD

| Dataset | Tipo | Acceso | Muestras |
|---------|------|--------|---------|
| TCGA-LUAD MAF | Mutaciones somáticas WXS (MuTect2) | Open access (GDC) | 585 casos |

- **Proyecto GDC:** TCGA-LUAD
- **Tipo de cáncer:** Adenocarcinoma de pulmón
- **Referencia:** GRCh38
- **Acceso:** Sin dbGaP — open access
- **URL:** https://portal.gdc.cancer.gov/projects/TCGA-LUAD

## Notas de descarga

- Los FASTQ de HCC1395 NO se versionan en Git (en .gitignore)
- El MAF de TCGA-LUAD tampoco se versiona (tamaño)
- Ver comandos de descarga en docs/descarga_datos.md
