# Metatranscriptomics and Metaproteomics

## Introduction

In this section, we will be using a new dataset present here:

```
ls /work/eml-course_bioinfo_metaomics/datasets/arsenic
```

It is composed of samples from soil-derived cultures having both:
* metagenomes (DNA reads) 
* metatranscriptomes (RNA reads) 
* metaproteomes (pre-processed protein spectral data)

The combination of these datasets is essential for gaining insights into the genomic potential and active gene expression within microbial communities.

## Previous Steps

Genome assembly, binning, and annotation have all been already done similarly to what we did in the previous section. 

In addition, metatranscriptomic reads have already been aligned to our metagenomes as follows:

1. **Alignment with Bowtie2:** Metatranscriptome reads are aligned to metagenome assemblies using Bowtie2. This step generates SAM files.

2. **Sorting with Sambamba:** SAM files are sorted into BAM files using Sambamba, facilitating downstream analyses.

In the previous sections of this course, you have completed similar alignment and sorting steps. So here we will start from the sorted alignments files to directly have the raw number of RNA reads mapping to the genes present in the metagenome assembly.

# Preparation

To avoid being lost between the datasets, you can create a new one: 

```
mkdir /scratch/$(whoami)/arsenic && cd /scratch/$(whoami)/arsenic
```

And go inside:

```
cd /scratch/$(whoami)/arsenic
```

## Feature Counting

To quantify the number of reads aligned to specific genomic features, such as coding sequences (CDS), we will use featureCounts.

# Prerequisites

Let's prepare our scripts:

```
nano 07_featurecounts.sh
```
```
#!/bin/bash

cd /data

input1="/data2/datasets/arsenic/metatranscriptome/reads"

input2="/data2/datasets/arsenic/metagenome/assemblies/prodigal"

input3="/data2/datasets/arsenic/metatranscriptome/alignments"

output2="featurecounts"

mkdir $output2

cd /data

featureCounts --verbose -T $SLURM_CPUS_PER_TASK -t CDS -g ID -a $input2/EA_WTA_prodigal_genes.gff -o $output2/EA_MT_on_EA_WTA_MG_featurecounts.tsv $input3/EA*sorted*.bam

featureCounts --verbose -T $SLURM_CPUS_PER_TASK -t CDS -g ID -a $input2/TSB_WTA_prodigal_genes.gff -o $output2/TSB_MT_on_TSB_WTA_MG_featurecounts.tsv $input3/TSB*sorted*.bam

```
# Feature Counting Commands

Here, we are using the alignments files (sorted bam files) to count the number of sample reads mapping to our predicted genes. 

Can you guess from what kind of assembly the genes are predicted?

```
featureCounts --verbose -T $SLURM_CPUS_PER_TASK -t CDS -g ID -a $input2/EA_WTA_prodigal_genes.gff -o $output2/EA_MT_on_EA_WTA_MG_featurecounts.tsv EA*sorted*.bam
```
* `-T` = processors (threads, default 1)
* `-t` = type of feature (e.g: exon)
* `-g` = gene ID
* `-a` = annotation file (here in GFF format)
* `-o` = output

# Launch the job

Now let's use a previous script to launch the job. We will use one that can activate an environment within a container and execute the script inside it. 

```
cp /work/eml-course_bioinfo_metaomics/scripts/05_metabolic.sh .
```
<details>
  <summary>Set the script</summary>
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --time 3:00:00
#SBATCH --account=bioinformatics-meta-omics1
#SBATCH --reservation=bioinformatics_meta-omics2
```
```
PATH_TO_THE_IMAGE="/work/eml-course_bioinfo_metaomics/images/featurecounts.sif"
```
```
MICROMAMBA_ENV="featurecounts"
```
Ensure that the main command is set correctly:
```
srun apptainer exec \
--bind $PATH_TO_HOST_WORKING_DIRECTORY:/data,"/work/eml-course_bioinfo_metaomics/":/data2 \
$PATH_TO_THE_IMAGE \
micromamba run -n $MICROMAMBA_ENV bash /data/$PATH_TO_YOUR_SCRIPT
```
</details>

The script should only take a few minutes to run. After this, you can have a look at the outputs.

```
less -S featurecounts/* 
```
You can do ":" and then press N to go to the next file.

The output should look similar to this:

```
# Program:featureCounts v1.5.3; Command:"featureCounts" "--verbose" "-T" "48" "-t" "CDS" "-g" "ID" "-a" "metagenome/assemblies/prodigal/EA_WTA_prodigal_genes.gff" "-o" "metatranscriptome/featurecounts/EA_MT_on_EA_WTA_MG_featurecounts.ts>
Geneid  Chr     Start   End     Strand  Length  metatranscriptome/alignments/EA_MT_WOA_G_1__EA_WTA.sorted.bam   metatranscriptome/alignments/EA_MT_WOA_G_2__EA_WTA.sorted.bam   metatranscriptome/alignments/EA_MT_WOA_G_3__EA_WTA.sorted.ba>
1_1     k119_1  1       204     +       204     2       0       1       0       0       0       0       0       3       0       0       0
2_1     k119_2  1       171     -       171     0       0       0       0       0       0       0       0       0       0       0       0
2_2     k119_2  180     389     -       210     0       0       0       0       0       0       0       0       0       0       0       0
3_1     k119_3  227     418     +       192     0       0       0       0       0       0       0       0       0       0       0       0
4_1     k119_4  1       600     +       600     297     160     252     16      16      20      310     227     275     19      23      42
```

## Differential expression

Now we will use those raw counts to find which genes might be differentially expressed at RNA level using the DeSeq2 package in R.

To do so, you need to copy on your laptop the featureCounts output as well as the Rscript inside: /work/eml-course_bioinfo_metaomics/scripts/deseq2.R
You will also need the following designed matrices present in: /work/eml-course_bioinfo_metaomics/datasets/arsenic/metatranscriptome/matrice/

Finally install the dependencies by doing:
```
install.packages(PACKAGE)
```
Or this if it is a bioconductor package:
```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PACKAGE")
```

The script provided includes steps such as data preprocessing, quality control, statistical testing, and result visualization, producing summary reports and plots for further analysis and interpretation.

Here are some explanations about the different parts: 

```
####--Getting today date--####

date<-Sys.Date()
time<-format(Sys.time(), "%H-%M-%S")

####--Setting working directory--####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####--Libraries--####

library(xlsx)
library(stringr)
library('DESeq2')
library("RColorBrewer")
library(pheatmap)
library(ReportingTools)
library(EnhancedVolcano)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("EnhancedVolcano")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

```

Here are just the libraries, that you have to install if you do not already have them.

```
####--Setting global variable--####

#data2=/work/eml-course_bioinfo_metaomics/datasets/arsenic
#raw=data2

#metatranscriptome/alignments

# MG
#Ends : protein_id bin kegg description MG_DNA_TPM(kallisto or feature counts) MT_TPM(kallisto) 

headers=c("gene_id","KO_id","KO_description","ncbi_genus_id","kegg_gene_id","GHOSTX_score")

#headers:https://www.kegg.jp/blastkoala/help_ghostkoala.html

# Column1	User's gene id
# Column2.0 K number (KO identifier) assinged
# Column2.1 K number (KO identifier) assinged
# Column3	Second level of the KEGG Organisms hierarchy <-Are in wide format /!\
# Column4	Third level of the KEGG Organisms hierarchy <-Are in wide format /!\
# Column5	Genus in the NCBI taxonomy
# Column6	KEGG GENES ID
# Column7	GHOSTX score

#KO Identifier (K number): Represents a functional category or orthologous group.
#KEGG GENES ID: Represents a unique identifier for an individual gene.

####--Loading data--####

#--From the metagenome

#-Ghostkoala

ghostkoala_table<-read.delim2(
  "raw/metagenome/annotation/ghostkoala/EA_WTA_ghostkoala.tsv",
  header=F,
  col.names = headers)

#Here we just clean a bit the data
ghostkoala_table<-ghostkoala_table[ghostkoala_table$gene_id!="",]

  ```
All of this part is just to load the annotation data, generated with ghostkoala. You will need to download the EA_WTA_ghostkoala.tsv present at:
/work/eml-course_bioinfo_metaomics/datasets/arsenic/

And place it at "raw/metagenome/annotation/ghostkoala" or else simply change the path to wherever you put it.

It is not mandatory to run the scripts but it would be good to have it.

The path after "raw" where data is located in the arsenic folder on the scitas.

``` 
#--From the metatranscriptome

#-Featurecounts

featurecounts_table<-read.delim2("raw/metatranscriptome/featurecounts/EA_MT_on_EA_WTA_MG_featurecounts.tsv",header=T,skip=1)
#Geneid	Chr	Start	End	Strand	Length Sample_1 Sample_2 Sample_3 etc.

colnames(featurecounts_table)
```

This part is to load the featurecounts output inside R.

```
#-Checkm

checkm_table<-read.csv("raw/metagenome/binning/checkm/EA_WTA_checkm_profile.csv",
                       header=T,
                       fileEncoding = "UTF-8-BOM")

colnames(checkm_table)[9]<-"bin_size_Mbp"
```

Optionally you can load checkm profile to find which MAs have genes differentially expressed between the conditions. 

Also, you can find more about the samples with the metadata files present inside the arsenic folders.

```
/work/eml-course_bioinfo_metaomics/datasets/arsenic/
├── metagenome
│   ├── annotation
│   │   ├── eggnog
│   │   │   └── metadata_eggnog.csv
│   │   ├── ghostkoala
│   │   │   └── metadata_ghostkoala.tsv
│   │   └── metaxa2
│   │       └── metadata_metaxa2.csv
│   ├── assemblies
│   │   ├── kallisto
│   │   │   └── metadata_kallisto.csv
│   │   ├── megahit
│   │   │   └── metadata_megahit.tsv
│   │   ├── multiqc
│   │   │   ├── EA_WOA
│   │   │   │   └── multiqc_data
│   │   │   ├── EA_WTA
│   │   │   │   └── multiqc_data
│   │   │   ├── TSB_WOA
│   │   │   │   └── multiqc_data
│   │   │   └── TSB_WTA
│   │   │       └── multiqc_data
│   │   └── prodigal
│   └── binning
│       ├── checkm
│       │   └── metadata_checkm.csv
│       ├── contigs2bins
│       │   ├── EA_WOA_metawrap_bins
│       │   ├── EA_WTA_metawrap_bins
│       │   ├── metadata_contigs2bins.csv
│       │   ├── TSB_WOA_metawrap_bins
│       │   └── TSB_WTA_metawrap_bins
│       ├── drep
│       └── ghostkoala2bins
│           └── metadata_ghostkoala2bins.csv
├── metaproteome
│   ├── differential_expression
│   │   ├── EA
│   │   │   └── metadata_differential_expression_EA.csv
│   │   └── TSB
│   │       └── metadata_differential_expression_TSB.csv
│   └── protein_abundance
│       └── metadata_protein_abundance.tsv
└── metatranscriptome
    ├── alignments
    ├── differential_expression
    │   ├── EA
    │   │   └── metadata_differential_expression_EA.csv
    │   └── TSB
    │       └── metadata_differential_expression_TSB.csv
    ├── featurecounts
    ├── featurecounts_ghostkoala
    │   └── metadata_featurescounts.csv
    ├── matrice
    ├── multiqc
    └── reads
        ├── metatranscriptome_G
        └── metatranscriptome_R
            └── metadata_metatranscriptome_R.csv
```


```
#--From the proteome

protein_table<-read.csv("raw/EA_DE_proteins.xlsx",header = F)

```

You can also include one of the proteomic files inside /metaproteome/differential_expression to see if the expression is the same at the protein level.



## Proteomic part

In this dataset protein abundance has been acquired by sending the prodigal output to another laboratory specialized in high pressure liquid chromatography coupled to inductively-coupled plasma mass spectrometry.

Here are the steps that have been applied to the data:
* Data Transformation: Protein abundance values were transformed using a log2 function.
* Normalization: Values were normalized among biological replicates using LOESS normalization, and mean-centered across all conditions.
* Data Filtering: Proteins without abundance values in at least two of the biological triplicates in at least one condition were removed.
* Missing Data Imputation: Any remaining missing data were replaced with random numbers drawn from a normal distribution, using specific parameters.
* Statistical Analysis: Differentially abundant proteins were identified using Student’s t-test method with an adjusted q-value threshold of ≤0.05.
* Additional Filtering: Proteins were further filtered based on an absolute log2 fold change threshold of ≥1.

We will not go into more detail in this course and focus more on visualizing the results.

So similarly to the differential abundance analysis, you can use these results to see which genes are expressed at protein levels.


