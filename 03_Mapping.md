# Mapping
We're going to use [Strobealign](https://github.com/ksahlin/strobealign) to map your QC reads to your assembled contigs and [samtools](https://github.com/samtools/samtools) to format our alignments.

We'll do cross mapping, which is aligning all reads to all samples.

First, we need to make a directory (`mkdir`) for our alignments and a subdirectory for each of our samples. We can do that with a quick `for loop`:

```
for sample in $(cat samples.txt); do mkdir -p 03_Mapping/${sample}; done
```

### 1. Perform mapping with Strobealign

This time we'll use the `for loop` to iterate through each set of paired QC reads. We'll run this command for each of our samples. Make sure you change the parts of the command that say *`your_sample`* to your sample name.

```
nano 03_strobealign.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

    for reads in $(cat samples.txt)

    do

    strobealign 02_Assembly/filtered_contigs/${sample}_contigs.fa \
    01_ReadQC/fastp_reads/${reads}_R1.fastq.gz 01_ReadQC/fastp_reads/${reads}_R2.fastq.gz \
    -t $SLURM_CPUS_PER_TASK | samtools sort -o 03_Mapping/${sample}/${reads}.sorted.bam

    samtools index 03_Mapping/${sample}/${reads}.sorted.bam

    done

done
```

* `strobealign ... | samtools sort -o ...` Calls the 'strobealign' command to align reads to filtered contigs, and then pipes (`|`) the output to `samtools sort` to create a sorted BAM file.
* `-t` Specifies the number of threads for parallel processing.
* `samtools index` Creates an index for the sorted BAM file using

Again, you'll need to create an sbatch script to submit this job. The image should be `samtools-strobealign.sif` and your script should be `03_strobealign.sh`

We are being liberal with 10 hours to make sure it will finish but it will likely take less time.
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 10:00:00
#SBATCH --mem=140G
```

You can check on the progress of the script by looking to see that files are being produced in **`03_Mapping/sample`**.


When strobealign is finished, you should have a directory in **`03_Mapping/sample`** for each of your samples. In each sample directory there should be the same number of **`.bam`** (binary alignment) and **`.bam.bai`** (indexed alignment) files as you have samples. We'll need each of these files for binning.

### 2. Check the proportion of reads that aligned to our assembly.

We can use samtools to check the proportion of reads that mapped to our contigs.

```
nano 03_samtools.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

samtools flagstat 03_Mapping/${sample}/${sample}.sorted.bam -O tsv > 03_Mapping/${sample}/${sample}_stats.tsv

done
```

This command shouldn't take long and you can use the same image `samtools-strobealign.sif`, just change the name of your script to `03_samtools.sh`:
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --time 2:00:00
#SBATCH --mem=35G
```

You can view the proportion of mapped reads by typing `more 03_Mapping/your_sample/your_sample_stats.tsv`. "Mapped" shows the count as a percentage of the total number of QC-passed or QC-failed reads after the category name e.g.,

```
29279394    0   total (QC-passed reads + QC-failed reads)
24397974    0   mapped
83.33%      N/A primary mapped %
```

If you want a copy of these stats you can download them (make sure you are in a new terminal window not logged in to the server):

```
scp 'username@jed.epfl.ch:/scratch/username/your_dataset/03_Mapping/your_sample/*stats.tsv' .
```

**Q: What proportion of your reads were mapped to your assembly? Do you consider that "good"?**

You do not need to wait for this to finish before you move onto binning.