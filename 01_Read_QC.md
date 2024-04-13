# Read QC
We will start by evaluating the sequence quality of the raw reads which are provided in zipped fastq files (extension **`.fastq.gz`**).

## FastQC on raw reads
[FastQC](https://anaconda.org/bioconda/fastqc) performs quality control checks on raw sequence data and can give a quick impression of your data before doing further analysis. We'll also use [MultiQC](https://anaconda.org/bioconda/multiqc) which summarises the FASTQC output from multiple samples.

### 1. Make a directory for read quality control

Check you are in your working directory.
```
pwd
```
Make a folder for read quality control.
```
mkdir -p 01_ReadQC/fastqc_pass1
```
* the **`-p`** flag makes the parent directory if it does not already exist

### 2. Run FastQC
We can run fastqc on each of our samples using a for loop.

```
for sample in $(cat samples.txt)

do

fastqc 00_raw/${sample}_{R1,R2}.fastq.gz -o 01_ReadQC/fastqc_pass1 --threads 20

done
```
In the directory `01_ReadQC/fastqc_pass1` we now have a **`.zip`** and a **`.html`** for each sample. We can summarise each of those files using MultiQC.

### 3. Run MultiQC
We only need to tell MultiQC where to look for the FastQC files and where to output the results.
```
multiqc --outdir 01_ReadQC/fastqc_pass1/ 01_ReadQC/fastqc_pass1/
```
Now we'll look at the output by copying the **`.zip`** and a **`.html`** files to our laptops.

To copy a file from the server, open a new terminal window and navigate to the directory you would like to save the results in locally.
```
cd Directory/on/your/laptop
```
We can now use the following command copy the files. You will be prompted for your password.
```
scp 'username@jed.epfl.ch:/home/username/01_ReadQC/fastqc_pass1/*.html .'
```

You can now view the output interactively. You can look at the output for each sample by opening the **`your_sample_fastqc.html`**. You can also run look at the summary by opening **`multiqc_report.html`**.

## Quality filtering with fastp
We'll now perform quality filtering on the reads using [fastp](https://github.com/OpenGene/fastp#adapters). This step will trim low-quality sequences using their phred score, remove adaptors, and _optionally_ remove PCR duplicates. Some of the key parameters fastp uses are outlined below.

* Adapters are detected for paired end data with `--detect_adapter_for_pe`
* Quality filtering is enabled. Default is phred 15 and can be adjusted with `--qualified_quality_phred`
* Length filtering is enabled. Default is min length 15
* PolyG tail trimming is enabled by default for NextSeq/NovaSeq data which is auto-detected
* PCR duplicate removal is disabled by default 

### 4. From your working directory, make a folder for the QC reads and fastp output

```
#make a directory for QC reads
mkdir 01_ReadQC/fastp_reads

#make a directory for the fastp report
mkdir 01_ReadQC/fastp_report
```
### 5. Run fastp

We'll use a for loop again to run fastp on all of our samples.

```
for sample in $(cat samples.txt)

do

fastp -i 00_raw/${sample}_R1.fastq.gz -I 00_raw/${sample}_R2.fastq.gz -o 01_ReadQC/fastp_reads/${sample}_R1.fastq.gz -O 01_ReadQC/fastp_reads/${sample}_R2.fastq.gz \
--detect_adapter_for_pe --qualified_quality_phred 20 --trim_tail1=1 --trim_tail2=1 \
-h 01_ReadQC/fastp_report/${sample}.fastp.html -j 01_ReadQC/fastp_report/${sample}.fastp.json

done
```

**Q: Take note of the flags we've chosen. What default settings have we changed and why?**

## FastQC on quality filtered reads
Now we've done quality control on our reads, we'll do a second pass with FastQC.

### 6. From your working directory, make a folder for the FastQC results
```
mkdir 01_ReadQC/fastqc_pass2
```

### 7. Run FastQC and MultiQC on the QC reads
We can use the same commands we used earlier, now on the quality-controlled reads.  We just need to change the input and output directories.

```
for sample in $(cat samples.txt)

do

fastqc 01_ReadQC/fastp_reads/${sample}_{R1,R2}.fastq.gz -o 01_ReadQC/fastqc_pass2 --threads 20

done
```
Summarise with MultiQC
```
multiqc --outdir 01_ReadQC/fastqc_pass1/ 01_ReadQC/fastqc_pass1/
```
To view your files, you can use the above `scp` command again with the directory locations changed.
```
scp 'username@jed.epfl.ch:/home/username/01_ReadQC/fastqc_pass2/*.html /location/on_your_laptop/'
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

**Next:** [02_Assembly](02_Assembly.md)