# Representative MAGs
We have bins from each of our samples but it would be useful to have a dereplicated (or non-redundant) set of bins as we've likely binned the same organism from multiple samples.

To do that we'll use a tool called [dRep](https://github.com/MrOlm/drep). DRep has great documentation [here](https://drep.readthedocs.io/en/latest/) if you want to read more about the process.

Let's make our new output directory:

`mkdir -p 06_Representative_MAGs/drep`

## Running dRep
DRep performs a binQC step using CheckM as part of its evaluation. Since we've already ran CheckM2 we'll supply dRep with that file, but we need to format it appropriately first:

```
nano 06_prepare_for_drep.sh
```

```
#!/bin/bash

# Concatenate all *filtered.tsv files from checkm2 into a single temporary file
cat 04_Binning/checkm2/*/*filtered.tsv > 06_Representative_MAGs/drep/checkm2.tmp

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

If you run `head 06_Representative_MAGs/drep/checkm2_for_drep.csv `

You should have a file that looks like this:

```
genome,completeness,contamination
KR46_June_concoct_bin.11.fa,94.21,0.52
KR46_June_concoct_bin.2.fa,88.45,1.91
KR46_June_concoct_bin.26.fa,77.23,3.49
KR46_June_concoct_bin.39.fa,77.09,3.45
```

Now make a directory for all of the MAGs you will dereplicate:

`mkdir 06_Representative_MAGs/MAGs`

And copy all of your MAGs to that directory:

```
for sample in $(cat samples.txt); do cp 04_Binning/dastool/${sample}/${sample}_DASTool_bins/*.fa 06_Representative_MAGs/MAGs/; done
```
We are using semicolons (`;`) here in our `for loop` as we are not running as a script, we'll just put the command directly in terminal.

You should now have a copy of your bins (.fa files) in the MAGs directory. Check with:

```
ls  06_Representative_MAGs/MAGs/
```

We can now prepare out dRep script:
```
nano 06_drep.sh
```

```
#!/bin/bash
cd /data

dRep dereplicate 06_Representative_MAGs/drep/out -g 06_Representative_MAGs/MAGs/*.fa -comp 50 -con 10 -sa 0.98 --genomeInfo 06_Representative_MAGs/drep/checkm2_for_drep.csv
```
* `-p` = processors (threads)
* `-comp` = minimum completeness
* `-con` = maximum contamination
* `-sa` = the Average Nucleotide Identity (ANI) of the clustering i.e., 98%

Example: dRep dereplicate output_dir/ -g /path/to/genomes/*.fasta

And prepare an sbatch script to submit the job:
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 02:00:00
#SBATCH --mem=140G
```

The output files we're interested in:
`06_Representative_MAGs/drep/data_tables/Wdb.csv`: Winning genomes, this file has your dereplicated MAG set

`06_Representative_MAGs/drep/data_tables/Cdb.csv`: Genomes and cluster designations. This file has all of the input MAGs and tells you which MAG is in which cluster.

You can download the datatables to your laptop with `scp`.

```
scp -r username@jed.epfl.sh:06_Representative_MAGs/drep/data_tables .
```

Now we'll create a concatenated file your dereplicated genomes:

```
cat 06_Representative_MAGs/drep/dereplicated_genomes/*.fa > 06_Representative_MAGs/MAGdb.fa
```

And we'll use strobealign again to map the reads to our set of representative MAGs.
```
mkdir 06_Representative_MAGs/reads_to_mags
```


```
nano 06_strobealign_repMAGs.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

strobealign 06_Representative_MAGs/MAGdb.fa 01_ReadsQC/fastp_reads/metagenomes/${sample}_R1.fastq.gz 01_ReadsQC/fastp_reads/${sample}_R2.fastq.gz -U -t 18 \
| samtools sort -o 06_Representative_MAGs/reads_to_mags/${sample}.sorted.bam

samtools index ${sample}.sorted.bam

done
```

Next we'll run CoverM to calculate the abundance of our MAGs

```
#!/bin/bash
cd /data

coverm genome -d 06_Representative_MAGs/drep/dereplicated_genomes/ -x fa -b 06_Representative_MAGs/reads_to_mags/*.bam -o 06_Representative_MAGs/coverm/coverm_dereplicated_mags_50.tsv --min-covered-fraction 50 --threads 20
```


