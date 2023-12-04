# Mapping
We're going to use [Strobealign](https://github.com/ksahlin/strobealign) to map your QC reads to your assembled contigs and [samtools](https://github.com/samtools/samtools) to format our alignments.

We're going to do cross mapping, which is aligning all reads to all samples.

First, lets make a directory for our alignments and for each of our samples.

`mkdir -p 03_Mapping/your_sample_name`

This time we'll use a for loop to iterate through each set of paired QC reads. We'll run this command for each of our samples. Make sure your change the parts of the command that say "your_sample".

```
for sample in $(cat samples.txt)

do

strobealign 02_Assembly/trimmed_contigs/your_sample_contigs.fa \
01_ReadQC/readsqc/${sample}_R1.fq.gz 01_ReadQC/readsqc/${sample}_R2.fq.gz \
-t 20 | samtools sort -o 03_Mapping/your_sample/${sample}.sorted.bam

samtools index 03_Mapping/your_sample/${sample}.sorted.bam

done
```

We can use samtools to check the proportion of reads that mapped to our contigs.
```
samtools flagstat your_sample.sorted.bam -O tsv > your_sample_stats.tsv
```

Mapped shows the count as a percentage of the total number of QC-passed or QC-failed reads after the category name e.g.,

```
29279394    0   total (QC-passed reads + QC-failed reads)
24397974    0   mapped
83.33%      N/A primary mapped %
```

**Q: What proportion of your reads were mapped to your assembly? Do you consider that "good"?** 

You should now have a directory in 03_Mapping for each of your samples. In each directory should be a file in the format ".bam" (binary alignment) and ".bam.bai" (indexed alignment). We'll need each of these files for binning.