# 16S rRNA gene amplicon analysis

For this workshop you will work with a dataset of amplicon libraries (2x250bp) sequenced on an Illumina MiSeq.

## [Requirements](00_Requirements.md)
* Install miniconda
* Install R Studio with packages DADA2 and Ampvis2

## Input files you'll be provided
* A set of R1 and R2 fastq files
* A list of sample names - **`samples.txt`**
* Sample metadata -  **`metadata.txt`**
* A cutadapt script - **`cutadapt.sh`**

## Trimming primers

For this step you need to know the length of your reads and the primers used for amplification. Our reads are 2x250bp and amplicons were generated with the primers 515F and 806R. 

We'll use cutadapt to remove primers.

### 1. Run cutadapt

Open terminal and create a conda environment for [cutadapt](https://cutadapt.readthedocs.io/en/stable/).
```
conda create -n cutadapt -c bioconda cutadapt
```

Once installed, activate the environment:

```
conda activate cutadapt
```
You can view and edit the cutadapt script in nano:
```
nano cutadapt_515_806.sh
```

To run the cutadapt script:
```
bash cutadapt_515_806.sh
```

We can look at the proportion of reads removed by looking at the file **`cutadapt_primer_trimming_stats.txt`** in the **`trimmed_reads`** directory. This directory also contains your trimmed reads.

```
more trimmed_reads/cutadapt_primer_trimming_stats.txt
```
The second column shows the fraction of reads retained in each sample and the third column shows the fraction of bps retained in each sample.

**Q: How many reads were retained?**

## Processing amplicon data with DADA2
We'll process the trimmed reads in R Studio with [DADA2](https://benjjneb.github.io/dada2/index.html) following the guidance in this [tutorial](https://benjjneb.github.io/dada2/tutorial.html).

DADA2 (Divisive Amplicon Denoising Algorithm 2) generates amplicon sequence variants (ASVs) from amplicon seqeuencing data.

First, open R Studio and load the required packages. If you haven't already, check the [requirements page](00_Requirements.md) for instructions how to install the R packages we'll use.

```
library(dada2)
library(tidyverse)
library(ampvis2)
```

### 2. Load our trimmed reads into R

First we'll get a copy of our sample names:
```
samples <- scan("samples.txt", what = "character")
```
Then we'll load the forward reads we trimmed with cutadapt:
```
forward_reads <- sort(list.files("trimmed_reads", pattern = "R1_trimmed.fastq.gz", full.names = TRUE))
```
And the reverse reads we trimmed with cutadapt:
```
reverse_reads <- sort(list.files("trimmed_reads", pattern = "R2_trimmed.fastq.gz", full.names = TRUE))
```
Now we can inspect read quality profiles for the first two samples:
```
plotQualityProfile(forward_reads[1:2])
plotQualityProfile(reverse_reads[1:2])
```
You can change the numbers in the square brackets to look at the read quality profiles of other samples.

**Q: Do you notice a difference in quality between the forward and reverse reads?**

### 3. Quality control the reads and learn error rates
We'll remove low-quality reads using the thresholds we specify below.
First we assign the filenames for the filtered fastq.gz files we'll generate:
```
filtered_forward_reads <- paste0(samples, "_R1_filt.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filt.fastq.gz")
```

Then perform quality filtering on the reads. It's important to adjust these parameters for your samples. We're using 2x250bp reads amplified with the primers 515f-806r.

```
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                maxN=0, rm.phix=TRUE, minLen=175, truncLen=c(210,190))

head(filtered_out)
```
* `maxEE` is a quality filtering threshold based on [expected errors](https://www.drive5.com/usearch/manual/exp_errs.html)
* `rm.phix` removes reads that match the PhiX bacteriophage genome added to Illumina sequencing runs for quality monitoring
* `minLen` is the minimum length of reads to keep after trimming
* `maxN` removes any sequences contains N
* `trunclen` sets the minimum size to trim the forward and reverse reads to keep the quality scores roughly above 30 overall

Dada2 will now learn the error rates for our sequencing data:
```
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
```
The error rates can be visualised:
```
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
```
The black line shows the estimated error rates and the red line shows the error rates expected. Ideally, the estimated error rates (black line) will be a good fit to the observed rates (points), and the error rates will drop with increased quality.

### 4. Dereplication and inferrings ASVs

Calculate unique reads and determine amplicon sequence variants:

```
dada_forward <- dada(filtered_forward_reads, err=err_forward_reads, pool=TRUE, multithread=TRUE)

dada_reverse <- dada(filtered_reverse_reads, err=err_reverse_reads, pool=TRUE, multithread=TRUE)
```

### 5. Merge paired reads

We can now merge the forward and reverse ASVs we calculated in the previous step and reconstruct the full target amplicon.
```
merged_reads <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse, filtered_reverse_reads, verbose=TRUE)

head(merged_reads[[1]])
```
The `nmatch` column in **`merged_reads`** shows the length of the overlap when we merge the reads. Our reads were 2x250bp of an 291bp region so we should have no problem overlapping.

### 6. Construct sequence table
We're now ready to make an ASV table:

```
seqtab <- makeSequenceTable(merged_reads)
dim(seqtab)
```
**Q: How many ASVs do we have?**

We should also inspect distribution of sequence lengths:
```
table(nchar(getSequences(seqtab)))
```
We're expecting amplicons to be ~252 bp, we can remove reads that are not in the expected size range, e.g.,

```
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:260]
table(nchar(getSequences(seqtab2)))
```

### 7. Remove chimeras
The next step is to remove chimeras - artificial sequences formed from two or more biological sequences joined together. [Here's](https://drive5.com/usearch/manual/chimeras.html) a good explanation of chimeras if you want to know more.

```
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```
**Q: How many chimeras were identifed?**

The total chimeras identified can often sound like a lot, but we can also check what they represented in terms of abundance:
```
sum(seqtab.nochim)/sum(seqtab2)
```
**Q: What percentage of reads do we have after chimera removal?**

### 8. Track reads through the pipeline
As a final progress check, we can check the number of reads that were removed during each step of the pipeline. If lot's of reads are getting lost at a certain step, we should go back and investigate what is happening.

```
getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_reads, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))


summary_tab
```

### 9. Assign taxonomy
Before we can assign taxonomy we have to download a database available [here](https://benjjneb.github.io/dada2/training.html). We're going to use Silva, but as you can see there are a number of options to choose from.

Once you've downloaded the database, we can tell R where to find it and assign taxonomy to our ASVs:
```
taxa <- assignTaxonomy(seqtab.nochim, "~/Databases/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
```
DADA2 also implements a method to make species level assignments, but only ASVs with exacts matches (100% identity) to sequenced reference strains will be assigned a species level classification:
```
taxa <- addSpecies(taxa, "~/Databases/silva_species_assignment_v138.1.fa.gz")
```
Now we can inspect the taxonomic assignments:
```
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

### 10. Write output
We'll follow the steps from [this tutorial](https://astrobiomike.github.io/amplicon/dada2_workflow_ex#extracting-the-standard-goods-from-dada2) to get a fasta file **`ASVs.fa`**, a count table **`ASVs_counts.tsv`**, and a taxonomy table **`ASVs_taxonomy.tsv`**. [This tutorial](https://astrobiomike.github.io/amplicon/dada2_workflow_ex#extracting-the-standard-goods-from-dada2) also discusses each of the steps in detail and is a useful resource if you want to learn more.

Change sequence headers to ASV_1, ASV_2...
```
#Extract column names (ASV sequences) from 'seqtab.nochim'
asv_seqs <- colnames(seqtab.nochim)

#Create a vector 'asv_headers' storing the name of each ASV
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

#Iterate over each column (ASV) in 'seqtab.nochim' and generate a header for each ASV in the format ">ASV_i"
for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}
```

Make and write out a fasta of our ASV sequences by combining 'asv_headers' and 'asv_seqs':
```
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")
```
Make and write out a counts table:
```
#'seqtab.nochim' is a table containing our sequence data
#'asv_headers' is a vector of header information associated with each sequence

#Transpose (t) the sequence table to make rows correspond to ASVs and columns to samples

asv_tab <- t(seqtab.nochim)

#Remove the ">" character from the beginning of each row name
row.names(asv_tab) <- sub(">", "", asv_headers)

#Iterate over each column (sample) in the ASV table and remove the "_R1_filt.fastq.gz" suffix
for (col in 1:ncol(asv_tab)){
    colnames(asv_tab)[col] <-  sub("_R1_filt.fastq.gz", "", colnames(asv_tab)[col])
}

#Write the ASV table to a tab-separated values (tsv) file
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
```
Make and write out a taxonomy table:
```
asv_tax <- taxa

#Remove the ">" character from the beginning of each row name in 'asv_tax'
row.names(asv_tax) <- sub(">", "", asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
```

You can find a copy of the expected output files on the github under **amplicon_output**.

## Visualise the results
There are lots of different ways to visualise amplicon data. Today we'll use [ampvis2](https://kasperskytte.github.io/ampvis2/). [Here](amplicon_output) is a copy of the data we generated above that we'll load into R to process with ampvis2.

First let's load the ASV table and add "ASV" to the first column header:
```
amp_asvtab <- read.delim("ASVs_counts.tsv")
colnames(amp_asvtab)[1]="ASV"
```

Next, we'll do the same for the taxonomy table:
```
amp_taxtab <- read.delim("ASVs_taxonomy.tsv")
colnames(amp_taxtab)[1]="ASV"
```

Ampvis can take a metadata file. This can include as many different paraments as you'd like to visualise but the first column should be "SampleID". Load the metadata file:
```
amp_meta <- read.delim("metadata.txt", header = TRUE)
```
Now we're ready to load the data into ampvis:
```
amp_data <- amp_load(
  otutable = amp_asvtab,
  metadata = amp_meta,
  taxonomy = amp_taxtab,
  fasta = "ASVs.fa"
)
```
We can inspect the ampvis object:
```
data
```
If everything loaded correctly it should look something like this:
```
ampvis2 object with 4 elements. 
Summary of OTU table:
     Samples         OTUs  Total#Reads    Min#Reads    Max#Reads Median#Reads    Avg#Reads 
          26         3986      2660996        65915       138104     104075.5       102346 

Assigned taxonomy:
     Kingdom       Phylum        Class        Order       Family        Genus      Species 
  3982(100%) 3760(94.33%) 3523(88.38%) 2988(74.96%) 1862(46.71%) 1186(29.75%)    60(1.51%) 

Metadata variables: 4 
 SampleID, Drillhole, PoreSize, Depth
 ```

 We can also inspect the fasta sequences we included to see how many ASVs we have and their sequence lengths:
```
amp_data$refseq
```

Let's see what taxa we have in our samples. We can do a quick check with the command `amp_heatmap` to get a first impression. This will show us the top 10 most abundant phyla.
```
amp_heatmap(amp_data)
```

But we can also change many of the settings with ampvis:
```
amp_heatmap(amp_data,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            tax_show = 25,
            facet_by = "Drillhole",
            color_vector = c("white", "navy"),
            plot_colorscale = "sqrt",
            plot_values = FALSE) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")
```

Check out the [ampvis documentation](https://kasperskytte.github.io/ampvis2/reference/index.html) to see more options! Have a go plotting a rarefaction curve with **`amp_rarecurve`**.

**Q: Are our samples well represented by our sequencing?**

Whether or not you should rarefy your data is often debated as it "throws away" data. Reading [this paper](https://doi.org/10.1128/msphere.00355-23) is recommended if you want to learn more about rarefying. Here's how we rarefy with ampvis:
```
amp_rare <- amp_filter_samples(amp_data, minreads = 65915, rarefy = 65915)
```

If you inspect `amp_rare` after running this command you'll see all samples now have 65915 reads. This number was chosen as it is our smallest library size.

There are many more visualisation options, check out the [ampvis documentation](https://kasperskytte.github.io/ampvis2/reference/index.html) and try making a boxplot or a ordination plot!

