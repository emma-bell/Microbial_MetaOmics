# 16S rRNA gene amplicon analysis

For this workshop you will be assigned a dataset amplicon libraries sequenced on an Illumina MiSeq.

Input = Illumina paired end reads (R1 and R2 fastq)

## Requirements
* Install conda and an environment with cutadapt
* Install R Studio and DADA2

## Trimming primers

For this step you need to know the length of your reads (e.g., 2x250bp) and the primers used for amplification. We'll then use cutadapt to remove them.

### 1. Run cutadapt

In terminal activate the conda environment with cutadapt installed.

```
conda activate cutadapt
```

Then run the cutadapt script.
```
bash cutadapt_515_806.sh
```

Column 2 shows the fraction of reads retained in each sample.
Column 3 shows the fraction of bps retained in each sample.

**Q: How many reads were retained?**

## Processing amplicon data with DADA2
We'll process the data in R Studio using [DADA2](https://benjjneb.github.io/dada2/index.html).