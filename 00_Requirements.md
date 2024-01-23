# Hands-on bioinformatics for microbial meta-omics (ENV-621).

In the hands-on part of this course, we'll learn how to go from raw Illumina short reads (shotgun metagenomics) to metagenome-assembled genomes annotated with taxonomy and function.

# Requirements

This course is for beginners in microbial meta-omics and there are no prerequisites. However, basic unix skills will be beneficial. If you are not familiar with unix I recommend running through [this resource](https://astrobiomike.github.io/unix/unix-intro) before the class.


To be able to participate
* You will need a SCITAS account and be able to connect via SSH following [these instructions](https://scitas-doc.epfl.ch/user-guide/using-clusters/connecting-to-the-clusters/).
* You will need to install the below tools on your laptop.

# Bioinformatics tools
## Installed on your laptop
* [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
* [R/RStudio](https://posit.co/download/rstudio-desktop/)
* To install packages in R, open RStudio and run the following commands in the console window:
    * [DADA2](https://benjjneb.github.io/dada2/index.html)

    To install Dada2 open RStudio and run the following command in the console window:
        ```
        if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install("dada2")
        ```

        Check the installation by loading Dada2:
        ```
        library(dada2)
        ```
        And checking which version is installed:
        ```
        packageVersion("dada2")
        ```
    * [ampvis2](https://kasperskytte.github.io/ampvis2/)
        ```
        install.packages("remotes")
        remotes::install_github("kasperskytte/ampvis2")
        ```
    * [Tidyverse](https://tidyverse.tidyverse.org)
        ```
        install.packages("tidyverse")
        ```
    * [Complex Heatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
        ```
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        BiocManager::install("ComplexHeatmap")
        ```

## Tools we'll use that are installed on SCITAS
* [FastQC](https://anaconda.org/bioconda/fastqc)
* [MultiQC](https://anaconda.org/bioconda/multiqc)
* [Fastp](https://anaconda.org/bioconda/fastp)
* [MEGAHIT](https://anaconda.org/bioconda/megahit)
* [SPAdes](https://github.com/ablab/spades)
* [seqkit](https://anaconda.org/bioconda/seqkit)
* [BBMap](https://anaconda.org/bioconda/bbmap)
* [Strobealign](https://github.com/ksahlin/strobealign)
* [CONCOCT](https://github.com/BinPro/CONCOCT)
* [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/)
* [DASTool](https://github.com/cmks/DAS_Tool)
* [DRep](https://github.com/MrOlm/drep)
* [CheckM](https://github.com/Ecogenomics/CheckM)
* [CheckM2](https://github.com/chklovski/CheckM2)
* [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
* [METABOLIC](https://github.com/AnantharamanLab/METABOLIC)
* [CoverM](https://github.com/wwood/CoverM)

# Set-up your working directory
After logging in to SCITAS set up your working directory.

```
#move into your working directory
cd /scratch/USERNAME

#copy the relevant "samples.txt" file to your working directory
cp samples.txt .

#make a folder for the raw dataset and move into it
mkdir 00_raw
cd 00_raw

#create a symbolic link from the raw sequence files to your working directory
ln -s YOUR_DATASET/*.fq.gz .
cd ..
```

# Overview of the metagenomic workflow
* Quality trim short reads and remove adaptors with fastp, perform basic QC with FastQC and look at multiple samples with MultiQC
* Assemble reads with MEGAHIT and SPAdes and compare the quality of each assembly with Quast
* Map reads with Strobealign and Samtools
* Perform binning with MetaBAT2 and CONCOCT
* Check the quality of genome bins with CheckM and CheckM2
* Refine bins with DAS Tool
* Dereplicate bins with dRep
* Assign taxonomy with GTDB-Tk
* Annotate function with METABOLIC
* Calculate relative abundance with CoverM