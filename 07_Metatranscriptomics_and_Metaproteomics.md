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

output2="featurecounts"

mkdir $output2

cd /data

featureCounts --verbose -T $SLURM_CPUS_PER_TASK -t CDS -g ID -a $input2/EA_WTA_prodigal_genes.gff -o $output2/EA_MT_on_EA_WTA_MG_featurecounts.tsv EA*sorted*.bam

featureCounts --verbose -T $SLURM_CPUS_PER_TASK -t CDS -g ID -a $input2/TSB_WTA_prodigal_genes.gff -o $output2/TSB_MT_on_TSB_WTA_MG_featurecounts.tsv TSB*sorted*.bam

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
You will also need the following designed matrices present in: /work/eml-course_bioinfo_metaomics/matrice

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


