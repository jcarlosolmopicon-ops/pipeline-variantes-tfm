# Datos de partida del proyecto

Ninguno de los ficheros descritos aquí se versiona en Git por tamaño (ver `.gitignore`).
Los comandos de descarga están en `docs/descarga_datos.md`.

## Validación del pipeline (exp03) — SEQC2 / HCC1395

| Muestra | Tipo | Accesión SRA | Cobertura | Tecnología |
|---------|------|--------------|-----------|------------|
| HCC1395 | Tumor (cáncer de mama triple negativo) | SRR7890824 | ~100x WES | Illumina HiSeq X |
| HCC1395BL | Normal pareado (linfocitos B de la misma donante) | SRR7890827 | ~80x WES | Illumina HiSeq X |

- **BioProject:** SRP162370
- **Consorcio:** SEQC2 Somatic Mutation Working Group
- **Publicaciones:** Fang et al., *Nat Biotechnol* 2021; Zhao et al., *Sci Data* 2021
- **Truth set:** https://sites.google.com/view/seqc2
- **Referencia:** GRCh38
- **Acceso:** abierto, sin dbGaP

### Truth set descargado

| Fichero | Variantes | HighConf | MedConf |
|---------|-----------|----------|---------|
| `highconf_sSNV.vcf.gz` | 39.560 | 37.398 | 2.162 |
| `highconf_sINDEL.vcf.gz` | 1.922 | 1.754 | 168 |

Se acompaña de `highconf_regions.bed`, que delimita las regiones evaluables. La validación
con hap.py emplea ambas categorías de confianza; el procedimiento de adaptación del truth
set está en NB-008 del lab notebook.

## Módulo de clasificación (exp04) — MAF somático pan-cáncer MC3 del TCGA

| Dataset | Tipo | Acceso | Contenido |
|---------|------|--------|-----------|
| MAF pan-cáncer MC3 | Mutaciones somáticas WXS, consenso de siete algoritmos | Abierto (GDC) | 3.598.760 variantes de 10.295 muestras tumorales, 33 tipos de cáncer |

- **Proyecto MC3:** Ellrott et al., *Cell Systems* 2018
- **Referencia:** GRCh38 tras liftover desde GRCh37 con CrossMap 0.7.2
- **Acceso:** abierto, sin dbGaP — los MAF de nivel abierto están enmascarados para
  variantes germinales
- **URL del portal:** https://portal.gdc.cancer.gov

**Denominación del fichero.** El MAF se descargó del GDC bajo el nombre del proyecto
TCGA-LUAD y las rutas del repositorio conservan esa denominación, pero su contenido real
corresponde al conjunto pan-cáncer MC3: el proyecto LUAD aporta solo el 6,1 % de las
variantes (219.815, de 509 muestras). La incidencia se detectó durante la verificación final
de los datos y se documenta en la sección 5.3.1 de la memoria. Obliga a interpretar el
clasificador como un modelo pan-cáncer, no específico de adenocarcinoma de pulmón.

## Genoma de referencia y recursos auxiliares

| Fichero | Uso |
|---------|-----|
| `referencia/GRCh38.fa` (+ `.fai`, `.dict`, índices BWA) | Alineamiento y llamada de variantes |
| `referencia/dbsnp138.vcf.gz` (+ `.tbi`) | Sitios conocidos para BQSR |
| `referencia/hg19ToHg38.over.chain.gz` | Liftover del MAF |
| Caché de VEP `~/.vep/homo_sapiens/110_GRCh38` | Anotación funcional en modo offline |

## Datos de prueba

`tests/` contiene un subconjunto reducido de lecturas (tumor y normal) para ejecutar el
pipeline de extremo a extremo en pocos minutos con el perfil `test`. Estos sí se versionan.
