# Liftover del MAF somático de GRCh37 a GRCh38

## Contexto

El MAF legacy descargado del GDC está en coordenadas GRCh37 (hg19). El pipeline trabaja
sobre GRCh38, por lo que fue necesario convertir las coordenadas antes de usar el fichero
como conjunto de entrenamiento del módulo de clasificación.

El MAF se descargó bajo la denominación del proyecto TCGA-LUAD, pero su contenido real
corresponde al conjunto pan-cáncer MC3 del TCGA (10.295 muestras tumorales de 33 tipos de
cáncer); el proyecto LUAD aporta el 6,1 % de las variantes. La incidencia se detectó durante
la verificación final de los datos y se documenta en la sección 5.3.1 de la memoria.

## Herramienta

CrossMap v0.7.2 — https://crossmap.readthedocs.io

Cadena de conversión: `hg19ToHg38.over.chain.gz` (UCSC Genome Browser),
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

## Procedimiento ejecutado

El módulo `maf` de CrossMap convierte el fichero en una sola pasada: reescribe las
coordenadas, actualiza el campo `NCBI_Build` y deja en un fichero `.unmap` aparte los
registros cuya posición no tiene correspondencia en la cadena de conversión.

```bash
CrossMap maf \
  datos-raw/referencia/hg19ToHg38.over.chain.gz \
  datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz \
  datos-raw/referencia/GRCh38.fa \
  GRCh38 \
  datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf

bgzip datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf
```

CrossMap antepone al fichero de salida una línea de comentario con los parámetros de la
ejecución, que es la razón por la que los scripts del módulo ML cargan el MAF con
`skiprows=1`:

```
#liftOver: Program=CrossMapv0.7.2, Time=May13,2026, ChainFile=.../hg19ToHg38.over.chain.gz, NewRefGenome=.../GRCh38.fa
```

## Resultados

| Etapa | Variantes | % |
|-------|-----------|---|
| Registros en el MAF original (GRCh37) | 3.600.963 | 100 % |
| No convertidos (fichero `.unmap`) | 2.203 | 0,061 % |
| Convertidos a GRCh38 y empleados en el modelado | 3.598.760 | 99,94 % |

Las tres cifras son verificables directamente sobre los ficheros:

```bash
# Registros del MAF original (la primera línea es la cabecera de columnas)
zcat datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz | tail -n +2 | wc -l

# Registros no convertidos (una cabecera + los registros descartados)
wc -l datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.unmap

# Registros del MAF convertido (línea 1: comentario de CrossMap; línea 2: cabecera)
zcat datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz | tail -n +3 | wc -l
```

## Ficheros

- **Entrada:** `datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz` (GRCh37)
- **Salida:** `datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz` (GRCh38)
- **Descartes:** `datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.unmap`
- **Cadena:** `datos-raw/referencia/hg19ToHg38.over.chain.gz`

Ninguno se versiona en Git por tamaño (ver `.gitignore`).
