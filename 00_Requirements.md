# Hands-on bioinformatics for microbial meta-omics (ENV-621).

In the hands-on part of this course, we'll learn how to go from raw Illumina short reads (shotgun metagenomics) to metagenome-assembled genomes annotated with taxonomy and function.

# Requirements

This course is for beginners in microbial meta-omics and there are no prerequisites. However, basic unix skills will be beneficial. If you are not familiar with unix I recommend running through [this resource](https://astrobiomike.github.io/unix/unix-intro) before the class.


To be able to participate
* You will need a SCITAS account and be able to connect via SSH following [these instructions](https://scitas-doc.epfl.ch/user-guide/using-clusters/connecting-to-the-clusters/). 
* You will need to install [R/RStudio](https://posit.co/download/rstudio-desktop/) on your laptop.

# Bioinformatics tools
## Installed on your laptop
* [R/RStudio](https://posit.co/download/rstudio-desktop/)

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
cd your_working_directory

#copy the relevant "samples.txt" file to your working directory
cp <your_dataset_sample_list.txt> <.>

#make a folder for the raw dataset and move into it
mkdir 00_raw
cd 00_raw

#create a symbolic link from the raw sequence files to your working directory
ln -s <path_to_your_dataset> <your_working_directory>
cd ..
```

# Overview of the workflow
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