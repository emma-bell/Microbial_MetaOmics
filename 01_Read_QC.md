# Read QC
We will start by evaluating the sequence quality of the raw reads which are provided in zipped fastq files (extension .fastq.gz).

## FastQC on raw reads
FastQC performs quality control checks on raw sequence data and can give a quick impression of your data before doing further analysis. We'll also use MultiQC which summarises the FASTQC output from multiple samples.

1. From your working directory, make a folder for read quality control.
```
mkdir -p 01_ReadQC/fastqc_pass1
#the -p flag makes the parent directory if it does not already exist
```

2. Run FastQC

We can run fastqc on all of our samples using a for loop.

```
for sample in $(cat samples.txt)

do

fastqc 00_raw/${sample}_{R1,R2}.fastq.gz -o 01_ReadQC/fastqc_pass1 --threads 20

done
```
You can look at the output from each sample by opening the "your_sample_fastqc.html" file. We can also run MultiQC to get a summary of the FastQC reports from all samples.

```
#move into the fastqc directory and run multiqc
cd 01_ReadQC/fastqc_pass1

#run multiqc in the current directory
multiqc .
```

Now the file 'multiqc_report.html' contains an overall view of your data. You will need to copy the directory "multiqc_data" with the html file to your laptop to view the output interactively.

## Quality filtering with fastp
We'll now perform quality filtering on the reads using [fastp](https://github.com/OpenGene/fastp#adapters). This step will trim low-quality sequences using their phred score, remove adaptors, and _optionally_ remove PCR duplicates. Some of the key parameters fastp uses are outlined below.

* Adapters are auto-detected by default 
* Quality filtering is enabled by default (with phred 15)
* Length filtering is enabled by default (min length 15)
* PolyG tail trimming is enabled by default for NextSeq/NovaSeq data which is auto-detected 
* PCR duplicate removal is disabled by default 

1. From your working directory, make a folder for the QC reads and fastp output

```
#make a directory for QC reads
mkdir 01_ReadQC/fastp_reads

#make a directory for the fastp report
mkdir 01_ReadQC/fastp_report
```
2. Run fastp

We'll use a for loop again to run fastp on all of our samples.

```
for sample in $(cat samples.txt)

do

fastp -i 00_raw/${sample}_R1.fastq.gz -I 00_raw/${sample}_R2.fastq.gz -o 01_ReadQC/fastp_reads/${sample}_R1.fastq.gz -O 01_ReadQC/fastp_reads/${sample}_R2.fastq.gz \
--qualified_quality_phred 20 --trim_tail=1 \
-h 01_ReadQC/fastp_report/${sample}.fastp.html -j 01_ReadQC/fastp_report/${sample}.fastp.json

done
```

**Q: Take note of the flags we've chosen. What default settings have we changed and why?**

## FastQC on quality filtered reads
Now we've done quality control on our reads, we'll do a second pass with FastQC.

1. From your working directory, make a folder for the FastQC results
```
mkdir 01_ReadQC/fastqc_pass2
```

2. Run FastQC and MultiQC on the QC reads
We can use the same commands we used earlier, now on the quality-controlled reads.  

```
for sample in $(cat samples.txt)

do

fastqc 01_ReadQC/fastp_reads/${sample}_{R1,R2}.fastq.gz -o 01_ReadQC/fastqc_pass2 --threads 20

done
```

```
cd 01_ReadQC/fastqc_pass2
multiqc .
```

**Q: What has changed in your report after running fastp?**

## What is a fastq file?
A fastq file contains both the biological (nucleotide) sequence and it's corresponding quality score. Each entry in a fastq files has 4 lines:

* A sequence identifier with information about the sequencing run and the cluster
* The sequence base calls (A, T, C, G)
* A separator (+)
* The base call quality scores

Here is an example:

```
@NB500928:54:H7WWCBGXY:1:11101:23433:1061 1:N:0:GATGTC
GGCGCNCCATTTCCAGGGCCTTTTGCGGTTGACCGGTCTGCATAAATGCCGTGGCGATCATCAGTCGCAGCTTCGGGTCCAGCTCCCGGCCCTGGTACAGGTTTTCGTGATCGGTCCAGTGCTGTGCCGCACCCTCGAAATCCTCTTGTTC
+
AAAAA#EEEEEEEEEEEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEAEEEEEEE<AEAEEEEEEEEEEEE<EEE<EEEEEAEEEE<AEAAAAEAAA<AEA<<AEAAAAEEA<<6/<<AAEEA/E/AAEA<EEE//
```

You can view your fastq file with the following command

```
gzcat <sample_name>.fastq.gz | head -20
#may be zcat depending on your operating system
```