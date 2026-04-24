## Problema identificado
El disco C: (WSL2, /dev/sdd) tiene capacidad limitada (~66 GB libres en Windows).
Los FASTQs de HCC1395 descomprimidos ocupan ~422 GB (tumor) + ~400 GB (normal).

## Solución adoptada
Disco duro externo K: (1.9 TB, ~890 GB libres) para datos pesados.

## Estructura en disco externo (K:)
/mnt/k/TFM-bioinformatica/
├── datos-raw/
│   ├── HCC1395/
│   │   ├── SRR7890824/   # Tumor .sra (65 GB)
│   │   └── SRR7890827/   # Normal .sra (70 GB)
│   ├── TCGA-LUAD/        # MAF GRCh37 + GRCh38
│   └── referencia/       # hg19ToHg38.over.chain.gz
└── datos-intermedios/    # Temporales fasterq-dump

## Montaje en WSL2
```bash
sudo mkdir -p /mnt/k
sudo mount -t drvfs K: /mnt/k
```

## Notas
- El montaje no es persistente — hay que ejecutarlo cada vez que se reinicia WSL2
- Los FASTQs finales irán también a /mnt/k/TFM-bioinformatica/datos-raw/HCC1395/
- El procesamiento (Nextflow workDir) permanece en disco interno ~/pipeline-variantes/datos-intermedios/
