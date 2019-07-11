# Drop-seq
Repository for Drop-seq Analysis

## Usage
1. Install
* Enter your information from [this link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
* Install cellranger **using your own token&id**
```
mkdir tools
cd tools

// Modify this code
wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=(put your token & id)"
tar exzvf cellranger-3.0.2.tar.gz
```
* If you want to get more references, see [this link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

```
qsub -l os7 -cwd path/to/preparation.sh
```

2. Run
```
// Example
qsub -cwd -l os7 drop-seq.sh -1 path/to/barcode.fastq.gz -2 path/to/transcript.fastq.gz -i ID_NAME -r REFERENCE -c EC
```

## Directory
```
.
├── README.md
├── data
│   ├── barcode.fastq.gz
│   └── transcript.fastq.gz
├── drop-seq.sh
├── output
│   └── sample1
├── preparation.sh
├── reference
│   ├── GRCh38
│   ├── Homo_sapiens
|   └── etc.
└── tools
    ├── 2.7.1a.tar.gz
    ├── homer
    └── etc.
```










