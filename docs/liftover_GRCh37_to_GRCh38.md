# Liftover del MAF somático de GRCh37 a GRCh38

## Contexto

El MAF legacy descargado del GDC está en coordenadas GRCh37 (hg19). El pipeline trabaja
sobre GRCh38, por lo que fue necesario convertir las coordenadas antes de usar el fichero
como conjunto de entrenamiento del módulo de clasificación.

> **Nota sobre el contenido del fichero.** El MAF se descargó bajo la denominación del
> proyecto TCGA-LUAD, pero su contenido real corresponde al conjunto pan-cáncer MC3 del
> TCGA (10.295 muestras tumorales de 33 tipos de cáncer). El proyecto LUAD aporta solo el
> 6,1 % de las variantes. La incidencia se detectó durante la verificación final de los
> datos y se documenta en la sección 5.3.1 de la memoria.

## Herramienta

CrossMap v0.7.2 — https://crossmap.readthedocs.io

## Cadena de conversión

`hg19ToHg38.over.chain.gz` (UCSC Genome Browser)
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

## Resultados

| Etapa | Variantes |
|-------|-----------|
| Originales (GRCh37) | 3.600.963 |
| Convertidas por CrossMap (GRCh38) | 3.600.605 |
| Omitidas en la conversión | 358 (0,01 %) |
| En el MAF final reconstruido | 3.598.760 |

La diferencia de 1.845 registros entre las variantes convertidas y las del MAF final se
produce durante la reconstrucción del fichero con las nuevas coordenadas: se descartan las
filas cuyas coordenadas convertidas no pudieron reasociarse de forma unívoca al registro
original. El conjunto empleado en el modelado es el de 3.598.760 variantes.

## Procedimiento ejecutado

```bash
# 1. Extraer coordenadas del MAF original a formato BED
#    La cuarta columna (NR) actúa como identificador de fila para el reensamblado posterior
zcat TCGA-LUAD.mutect2.somatic.maf.gz \
  | grep -v "^#" | tail -n +2 \
  | awk 'BEGIN{OFS="\t"} {print $5, $6-1, $7, NR}' \
  > /tmp/maf_coords.bed

# 2. Liftover de las coordenadas
CrossMap bed \
  datos-raw/referencia/hg19ToHg38.over.chain.gz \
  /tmp/maf_coords.bed \
  /tmp/maf_coords_hg38.bed

# 3. Reconstruir el MAF sustituyendo las coordenadas originales por las convertidas,
#    emparejando por el identificador de fila de la cuarta columna del BED
```

> **Pendiente de consolidar.** El paso 3 se realizó de forma interactiva durante el
> desarrollo. Para reproducibilidad completa conviene extraer ese código a
> `scripts/liftover_maf.py` y versionarlo, tal como se recoge entre las tareas pendientes
> del proyecto.

## Ficheros

- **Entrada:** `datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.somatic.maf.gz` (GRCh37)
- **Salida:** `datos-raw/TCGA-LUAD/TCGA-LUAD.mutect2.GRCh38.maf.gz` (GRCh38)
- **Cadena:** `datos-raw/referencia/hg19ToHg38.over.chain.gz`

## Referencia

Zhao, H., Sun, Z., Wang, J., Huang, H., Kocher, J. P., y Wang, L. (2014). CrossMap: a
versatile tool for coordinate conversion between genome assemblies. *Bioinformatics*,
30(7), 1006-1007. https://doi.org/10.1093/bioinformatics/btt730
