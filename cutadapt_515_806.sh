#cutadapt --version 4.4

mkdir -p trimmed_reads

for sample in $(cat samples.txt)
do

	echo "Running cutadapt on sample: $sample"

	cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
	-A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
    	-m 200 -M 235 --discard-untrimmed \
    	-o trimmed_reads/${sample}_R1_trimmed.fastq.gz -p trimmed_reads/${sample}_R2_trimmed.fastq.gz \
    	raw_reads/${sample}_R1.fastq raw_reads/${sample}_R2.fastq \
    	>> trimmed_reads/cutadapt_primer_trimming_stats.txt 2>&1

done

paste samples.txt <(grep "passing" trimmed_reads/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" trimmed_reads/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
