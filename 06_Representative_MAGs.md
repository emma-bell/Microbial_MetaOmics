# Representative MAGs
We have bins from each of our samples but it would be useful to have a dereplicated (or non-redundant) set of bins as we've likely binned the same organism from multiple samples.

To do that we'll use a tool called [dRep](https://github.com/MrOlm/drep). DRep has great documentation [here](https://drep.readthedocs.io/en/latest/) if you want to read more about the process.

Let's start by making our new output directory:

`mkdir -p 06_Representative_MAGs/drep`

## Dereplication with dRep
DRep performs a binQC step using CheckM as part of its evaluation. Since we've already ran CheckM2 we'll supply dRep with our CheckM2 results but we need to format it appropriately first:

```
nano 06_prepare_checkm2_for_drep.sh
```

```
#!/bin/bash

# Concatenate all checkm2 quality reports into a single temporary file

for sample in $(cat samples.txt)

do

    cat 04_Binning/checkm2/${sample}/quality_report.tsv | tail -n+2 >> 06_Representative_MAGs/drep/checkm2.tmp

done


# Extract columns 1 (sample name), 2 (completeness), and 3 (contamination) from the concatenated file, add ".fa," to the sample names, and redirect to a temporary file
awk '{print $1, $2, $3}' OFS="," 06_Representative_MAGs/drep/checkm2.tmp | sed 's/,/.fa,/' > 06_Representative_MAGs/drep/checkm2_for_drep.tmp


# Define headers
headers=genome,completeness,contamination
echo $headers > 06_Representative_MAGs/drep/headers.tmp


# Concatenate the header file and the modified checkm2_for_drep file into a CSV file. This is the file we keep for dRep.
cat 06_Representative_MAGs/drep/headers.tmp 06_Representative_MAGs/drep/checkm2_for_drep.tmp > 06_Representative_MAGs/drep/checkm2_for_drep.csv


# Remove temporary files
rm 06_Representative_MAGs/drep/*.tmp

```

You can run the above script directly in your terminal with `bash`.

Afterwards, check that everything looks as expected with:

```
head 06_Representative_MAGs/drep/checkm2_for_drep.csv
```

You should have a comma separated file that looks like this:

```
genome,completeness,contamination
KR46_June_concoct_bin.11.fa,94.21,0.52
KR46_June_concoct_bin.2.fa,88.45,1.91
KR46_June_concoct_bin.26.fa,77.23,3.49
KR46_June_concoct_bin.39.fa,77.09,3.45
```

Look good? Ok, now make a directory for all of your MAGs from all samples, which we will then dereplicate:

`mkdir 06_Representative_MAGs/MAGs`

And make a copy all of your MAGs in directory:

```
for sample in $(cat samples.txt); do cp 04_Binning/dastool/${sample}/${sample}_DASTool_bins/*.fa 06_Representative_MAGs/MAGs/; done
```
* We are using semicolons (`;`) in our `for loop` here where we would usually start a new line. The semicolon acts as a separator in the same way but also makes it easier to type the command directly into the terminal when we are not running from a script.

You should now have a copy of your bins (.fa files) in the MAGs directory. Check with:

```
ls  06_Representative_MAGs/MAGs/
```

Want to know how many MAGs you have? You can quickly count the number of files in a directory like this:
```
ls -1 06_Representative_MAGs/MAGs/ | wc -l
```

Now we can prepare our dRep script:

```
nano 06_drep.sh
```

```
#!/bin/bash
cd /data

dRep dereplicate 06_Representative_MAGs/drep/out -p $SLURM_CPUS_PER_TASK -g 06_Representative_MAGs/MAGs/*.fa -comp 50 -con 10 -sa 0.98 --genomeInfo 06_Representative_MAGs/drep/checkm2_for_drep.csv
```
* We put the command `dRep dereplicate` followed by the output directory (which will be created by dRep)
* `-g` = the location of our bins
* `-p` = processors (threads, default 6)
* `-comp` = minimum completeness
* `-con` = maximum contamination
* `-sa` = the Average Nucleotide Identity (ANI) of the clustering i.e., 98%

And prepare an sbatch script with the image `drep.sif`:
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 04:00:00
#SBATCH --mem=140G
```

This will produce 5 folders in your output directory:

```
06_Representative_MAGs/drep/out/
├── data
├── data_tables
├── dereplicated_genomes
├── figures
└── log
```

You should take some time to look through each of the folders. Some of the key files we're interested in are:
* `06_Representative_MAGs/drep/data_tables/Wdb.csv`: Winning genomes, this file has your dereplicated MAG set
* `06_Representative_MAGs/drep/data_tables/Cdb.csv`: Genomes and cluster designations. This file has all of the input MAGs and tells you which MAG is in which cluster.

You can download the `data_tables` and `figures` directories to your laptop with `scp -r`.

```
scp -r username@jed.epfl.sh:06_Representative_MAGs/drep/data_tables a/directory/on/your/laptop
```

## Calculating the relative abundance of dereplicated MAGs

We now have a dereplicated set of MAGs (`06_Representative_MAGs/drep/dereplicated_genomes/`) that represent all of the microorganisms we've recovered from our samples. It would be useful to know the relative abundance of each of those MAGs in each of our samples. To do that we're going to use [CoverM](https://wwood.github.io/CoverM/coverm-genome.html). CoverM takes BAM files or raw reads as input. We're going to prepare BAM files using [Strobealign](`https://github.com/ksahlin/strobealign`).

### Mapping reads to MAGs

First we need to create a concatenated file of your dereplicated genomes (i.e., put all of our dereplicated MAGs into one fasta file):

```
cat 06_Representative_MAGs/drep/dereplicated_genomes/*.fa > 06_Representative_MAGs/MAGdb.fa
```

Then we'll use Strobelign to map our reads to our representative MAGs. Let's make a directory for the output of Strobelign:

```
mkdir 06_Representative_MAGs/reads_to_mags
```

We have used Strobealign before to map our qc reads to our assembled contigs. Can you edit the [script](03_Mapping.md#1-perform-mapping-with-strobealign) we used previously to now map qc reads to the dereplicated MAGs?

Note: Last time we were mapping all reads to multiple samples. This time we will map all reads to the representative MAGs only (`06_Representative_MAGs/MAGdb.fa`).

```
nano 06_strobealign_repMAGs.sh
```
After you've tried on your own, you can check against the script here.
<details>
  <summary>Check script</summary>

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

strobealign 06_Representative_MAGs/MAGdb.fa 01_ReadQC/fastp_reads/${sample}_R1.fastq.gz 01_ReadQC/fastp_reads/${sample}_R2.fastq.gz -U -t $SLURM_CPUS_PER_TASK \
| samtools sort -o 06_Representative_MAGs/reads_to_mags/${sample}.sorted.bam

samtools index 06_Representative_MAGs/reads_to_mags/${sample}.sorted.bam

done
```

</details>


You'll need to prepare an sbatch script with the image `samtools-strobealign.sif`:
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 10:00:00
#SBATCH --mem=140G
```

### Calculating relative abundance with CoverM

Now we have sorted (`.sorted.bam`) and indexed (`.sorted.bam.bai`) BAM files in **`06_Representative_MAGs/reads_to_mags/`** we're ready to run CoverM.

```
mkdir 06_Representative_MAGs/coverm
```

```
#!/bin/bash
cd /data

coverm genome -d 06_Representative_MAGs/drep/dereplicated_genomes/ -x fa \
-b 06_Representative_MAGs/reads_to_mags/*.bam -o 06_Representative_MAGs/coverm/coverm_drep_mag_rel_abun.tsv \
--threads $SLURM_CPUS_PER_TASK
```
* `-d` : a directory containing FASTA files of each genome
* `-x` : tells CoverM that our genomes end in the extension `.fa`
* `-b` : the directory with our sorted BAM files
* `-o` : Where to output the results
* `--min-covered-fraction` : default is 10. Genomes with less covered bases are reported as 0.
* `--method` : default is `relative_abundance`. You can choose different methods which are explained [here](https://github.com/wwood/CoverM/#calculation-methods)

Prepare an sbatch script with the image `coverm.sif`:
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 2:00:00
#SBATCH --mem=140G
```

Copy the results (`06_Representative_MAGs/coverm/coverm_drep_mag_rel_abun.tsv`) to your laptop with `scp`.

Congratulations! You've just completed a metagenomic workflow from raw Illumina short reads through to taxonomically and functionally annotated metagenome-assembled genomes (MAGs)! You'll often hear this type of workflow referred to as "genome-resolved metagenomics". What happens now? What did you find? Don't worry, we're going to go through the outputs we've generated and think about different ways to present the data on Wednesday when we explore [data visualisation](08_Data_visualisation.md).

**Previous:** [05_Annotation](05_Annotation.md)