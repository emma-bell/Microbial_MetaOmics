# Mapping
We're going to use [Strobealign](https://github.com/ksahlin/strobealign) to map your QC reads to your assembled contigs and [samtools](https://github.com/samtools/samtools) to format our alignments.

We'll do cross mapping, which is aligning all reads to all samples.

First, lets make a directory for our alignments and for each of our samples.

`mkdir -p 03_Mapping/your_sample_name`

1. Perform mapping with Strobealign

This time we'll use a for loop to iterate through each set of paired QC reads. We'll run this command for each of our samples. Make sure you change the parts of the command that say *`your_sample`* to your sample name.

```
for reads in $(cat samples.txt)

do

strobealign 02_Assembly/trimmed_contigs/your_sample_contigs.fa \
01_ReadQC/readsqc/${reads}_R1.fq.gz 01_ReadQC/readsqc/${reads}_R2.fq.gz \
-t 20 | samtools sort -o 03_Mapping/your_sample/${reads}.sorted.bam

samtools index 03_Mapping/your_sample/${reads}.sorted.bam

done
```

You should now have a directory in **`03_Mapping`** for each of your samples. In each sample directory there should be the same number of **`.bam`** (binary alignment) and **`.bam.bai`** (indexed alignment) files as you have samples. We'll need each of these files for binning.

2. Check the proportion of reads that aligned to our assembly.

We can use samtools to check the proportion of reads that mapped to our contigs.

```
#From within your 03_Mapping/your_sample/ directory

samtools flagstat your_sample.sorted.bam -O tsv > your_sample_stats.tsv
```

You can view the proportion of mapped reads by typing `more your_sample_stats.tsv`. "Mapped" shows the count as a percentage of the total number of QC-passed or QC-failed reads after the category name e.g.,

```
29279394    0   total (QC-passed reads + QC-failed reads)
24397974    0   mapped
83.33%      N/A primary mapped %
```

**Q: What proportion of your reads were mapped to your assembly? Do you consider that "good"?**