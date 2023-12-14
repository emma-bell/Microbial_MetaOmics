# Binning

Now we have our contigs (.fasta) and our alignments (.sorted.bam) we can move onto binning.
We're going to use two assemblers, [CONCOCT](https://github.com/BinPro/CONCOCT) and [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/).

### 1. Let's start by making the directories we'll need for the output.
```
mkdir -p 04_Binning/concoct
mkdir 04_Binning/metabat
```

## CONCOCT
CONCOCT usually takes longer to run than MetaBAT2, so let's start with that one.
You'll also need a directory for each of your samples inside the **`04_Binning/concoct`** directory.

We'll do that with a quick loop that you can type and run directly in the terminal window:

```
for sample in $(cat samples.txt); do mkdir 04_Binning/concoct/${sample}; done
```

### 2. Running CONCOCT
Make our concoct script:
```
nano 04_concoct.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

#First the contigs are cut up into 10000 bp

cut_up_fasta.py 02_Assembly/filtered_contigs/${sample}_contigs.fa -c 10000 -o 0 --merge_last -b 04_Binning/concoct/${sample}/${sample}_contigs_10K.bed > 04_Binning/concoct/${sample}/${sample}_contigs_10K.fa

#Then a coverage table is made using our alignments

concoct_coverage_table.py 04_Binning/concoct/${sample}/${sample}_contigs_10K.bed 03_Mapping/${sample}/*.sorted.bam --threads > 04_Binning/concoct/${sample}/${sample}_coverage_table.tsv

#Now we're ready to run CONCOCT

concoct --composition_file --threads $SLURM_CPUS_PER_TASK 04_Binning/concoct/${sample}/${sample}_contigs_10K.fa --coverage_file 04_Binning/concoct/${sample}/${sample}_coverage_table.tsv -b 04_Binning/concoct/${sample}/${sample}

#Merge subcontig clustering into original contig clustering

merge_cutup_clustering.py 04_Binning/concoct/${sample}/${sample}_clustering_gt1000.csv > 04_Binning/concoct/${sample}/${sample}_clustering_merged.csv

#Extract the bins as individual fasta

mkdir 04_Binning/concoct/${sample}/bins_fasta
extract_fasta_bins.py 02_Assembly/filtered_contigs/${sample}_contigs.fa 04_Binning/concoct/${sample}/${sample}_clustering_merged.csv --output_path 04_Binning/concoct/${sample}/bins_fasta

done
```

Submit the job with an sbatch script, image: `concoct.sif`, script: `04_concoct.sh`

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 12:00:00
#SBATCH --mem=140G
```

Congratulations, you have bins (or will have soon)! You can find them in the **`04_Binning/concoct/${sample}/bins_fasta`** directory.

## MetaBAT2
We're now ready to do binning with MetaBAT2.

### 3. Running MetaBAT2
Again we'll need a metabat script:

```
nano 04_metabat.sh
```

```
#!/bin/bash
/data

for sample in $(cat samples.txt)

do

#First make a coverage file using our alignments

jgi_summarize_bam_contig_depths --outputDepth 04_Binning/metabat/${sample}/${sample}_depth.txt 03_Mapping/${sample}/*.sorted.bam

#Now we're ready to run MetaBAT

metabat2 -m 1500 -i 02_Assembly/filtered_contigs/${sample}_contigs.fa -a 04_Binning/metabat/${sample}/${sample}_depth.txt -o 04_Binning/metabat/${sample}/bins_fasta/${sample}_metabat.bin

#MetaBAT also has an option to summarise the bin depths

aggregateBinDepths.pl 04_Binning/metabat/${sample}/${sample}_depth.txt 04_Binning/metabat/${sample}/bins_fasta/*.fa > 04_Binning/metabat/${sample}/${sample}_bindepth.txt

done
```
You'll need to submit this with an sbatch script calling the metabat image `metabat.sif` and your metabat script `04_metabat.sh`. You can use most of the same SBATCH parameters you used for concoct - the exception is that metabat is much quicker so you can reduce the time to 2 hours (and it will still likely be quicker than that).

Congratulations, you now have more bins! You can find them in the **`04_binning/metabat/${sample}/bins_fasta`** directory.

## Tidying up before we move on to refinement
The sample name and binning method was added to our METABAT bin file names within the `metabat2` script but our concoct bins are currently called **`1.fa`**, **`2.fa`** etc.

Before we move onto bin refinement, let's add the sample name and binning method our CONCOCT bin file names. 

In your working directory make a script for renaming our concoct bins:

`nano 04_rename_concoct_bins_test.sh`

This for loop will show us what the renamed files will look like, without making a change (so we can check it looks good before proceeding):

```
for sample in $(cat samples.txt)

do

        for bin in 04_Binning/concoct/${sample}/bins_fasta/*.fa

        do

        binid=$(echo ${bin} | cut -f5 -d "/")
        location=$(echo ${bin} | cut -f1-4 -d "/")

        echo ${location}/${sample}_concoct.${binid}

        done

done
```

You can run this script directly in your terminal with bash:

```
bash 04_rename_concoct_bins_test.sh
```

Did everything look right? Your screen should output the new names of your bin files. Check the path is correct and your bin names make sense. It should look something like this:

```
04_Binning/test/concoct/KR46_June/bins_fasta/KR46_June_concoct.1.fa
04_Binning/test/concoct/KR46_June/bins_fasta/KR46_June_concoct.2.fa
04_Binning/test/concoct/KR46_June/bins_fasta/KR46_June_concoct.3.fa
04_Binning/test/concoct/KR46_May/bins_fasta/KR46_May_concoct.1.fa
04_Binning/test/concoct/KR46_May/bins_fasta/KR46_May_concoct.2.fa
04_Binning/test/concoct/KR46_May/bins_fasta/KR46_May_concoct.3.fa
04_Binning/test/concoct/KR46_Sept/bins_fasta/KR46_Sept_concoct.1.fa
04_Binning/test/concoct/KR46_Sept/bins_fasta/KR46_Sept_concoct.2.fa
04_Binning/test/concoct/KR46_Sept/bins_fasta/KR46_Sept_concoct.3.fa
```
If so, we can edit the command to commit to the change.

```
cp 04_rename_concoct_bins_test.sh 04_rename_concoct_bins.sh
```
Then open your copied script and edit the line below the # comment:
```
nano 04_rename_concoct_bins.sh
```

```
for sample in $(cat samples.txt)

do

        for bin in 04_Binning/concoct/${sample}/bins_fasta/*.fa

        do

	    binid=$(echo ${bin} | cut -f5 -d "/")
        location=$(echo ${bin} | cut -f1-4 -d "/")

        #this is the only line we will change
        mv ${bin} ${location}/${sample}_concoct.${binid}

        done

done
```

**What did we change?**

**`04_rename_concoct_bins_test.sh`**: **`echo`** shows the modified file paths but does not actually move or rename any files.

**`04_rename_concoct_bins.sh`**: Moves the bin files with **`mv`** to new locations with the modified file names.

We're now ready to move on to bin refinement.

## Bin refinement with DAS Tool

### 6. Performing bin refinement 

Before we can run [DAS Tool](https://github.com/cmks/DAS_Tool) we need to prepare the input files. Let's start by making a directory for the results.

`mkdir 04_Binning/dastool`

The input files we need are a list of contig headers and their associated bin. We can get this list using a script from the DAS Tool package. Let's make a script to do that:
```
nano 04_prepare_dastool.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

#Get the CONCOCT contig input file

Fasta_to_Contig2Bin.sh -e fa -i 04_Binning/concoct/${sample}/bins_fasta/ > 04_Binning/dastool/${sample}.concoct.contigs2bin.tsv

#Get the MetaBAT2 contig input file

Fasta_to_Contig2Bin.sh -e fa -i 04_Binning/metabat/${sample}/bins_fasta/ > 04_Binning/dastool/${sample}.metabat.contigs2bin.tsv

done
```
This won't take long so let's try running the job interactively -
**Make sure to change the PATH to YOUR_USERNAME and YOUR_DATASET**:

```
srun -A bioinformatics-meta-omics1 --pty apptainer shell --bind /scratch/ebell/Subsurface_KR11:/data /home/nljacque/images/dastool.sif
```

You should now be able to execute the dastool script interactively with:

```
bash 04_prepare_dastool.sh
```
You can type `exit` to exit the interactive apptainer.

You should now have two files for each of your samples:

**`04_Binning/dastool/${sample}.metabat.contigs2bin.tsv`**
**`04_Binning/dastool/${sample}.metabat.contigs2bin.tsv`**

Example of what they look like:
```
KR46_June_NODE_131_length_46115_cov_94.653061	KR46_June_concoct_bin.0
KR46_June_NODE_141_length_44438_cov_58.831016	KR46_June_concoct_bin.0
KR46_June_NODE_376_length_17039_cov_65.949187	KR46_June_concoct_bin.0
```
The first column in your **`contig_ID`**, the second column in the **`bin`** that contig was assigned too.

Now we're ready to make our Das Tool script:
```
nano 04_dastool.sh
```

Make sure to check the name of your contigs file in your working directory as some of you have contigs called **`sample.fa`** and some of you have contigs called **`sample.contigs.fa`** etc.

Then in the `-c` flag below use the appropriate contig file name:

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

DAS_Tool -i 04_Binning/dastool/${sample}.concoct.contigs2bin.tsv,04_Binning/dastool/${sample}.metabat.contigs2bin.tsv -l concoct,metabat -c 02_Assembly/filtered_contigs/${sample}.contigs.fa -o 04_Binning/dastool/${sample}/${sample} --write_bins --threads $SLURM_CPUS_PER_TASK

done
```
We'll submit this as a batch script. Change the image to **`dastool.sif`** and the script to **`04_dastool.sh`**

You can use the parameters below:

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 6:00:00
#SBATCH --mem=140G
```

This shouldn't take 6 hours to run, but now would be a good time to take a coffee break :)

When your Das Tool run has finished, you should have a file **`04_Binning/dastool/${sample}/${sample}_DASTool_summary.tsv`** which summarises the results. If you have any bins ending with *`_sub`* in the **`04_Binning/dastool/${sample}/${sample}_DASTool_bins`** directory, these bins have been refined by DAS Tool.

Now let's do CheckM2 on our refined bins to get some quality measures.

# Bin QC
We want to know how complete our bins are and whether they contain any contamination. To do this we're going to use [CheckM2](https://github.com/chklovski/CheckM2). We'll also run [CheckM](https://github.com/Ecogenomics/CheckM) to compare.

## CheckM2

We'll use the "predict" module of CheckM2 to get the completeness and contamination of our Das Tool bins.

### 4. Running CheckM2

Let's make a directory for our bin QC:

`mkdir 04_Binning/checkm2`

Now lets make our checkm2 script:
```
nano 04_checkm2.sh
```

```
!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

#CheckM2 on DASTool bins

~/checkm2/bin/checkm2 predict --threads $SLURM_CPUS_PER_TASK --input 04_Binning/dastool/${sample}/${sample}_DASTool_bins/ -x fa --output-directory 04_Binning/checkm2/${sample}/

done
```

CheckM2 works on a directory of genome bins in FASTA format. By default, CheckM2 assumes these files end with the extension ‘fna’, so we are changing it to `fa` with the `–x` flag.

Create an sbatch script to submit the job.
```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --time 6:00:00
#SBATCH --mem=70G
```

The file **`04_Binning/checkm2/${sample}/quality_report.tsv`** contains your results. You should have one for each sample and each binning method. We can select the bins that pass our minimum completeness (50%) and maximum contamination (10%) thresholds with an `awk` command.

```
nano 04_filter_checkm2_table.sh
```

```
#!/bin/bash

for sample in $(cat samples.txt)

do

awk -F "\t" '$2 >=50 && $3 <=10' 04_Binning/checkm2/${sample}/quality_report.tsv > 04_Binning/checkm2/${sample}/${sample}_quality_report_filtered.tsv

done
```
* `awk` is widely used for text processing in terminal
* `-F "\t"` is telling the command that our columns in `quality_report.tsv` are separated by a tab (`"\t"`)
* `'$2 >=50 && $3 <=10'` is looking for values >50 in column 2 (i.e., completeness >50) AND values <10 in column 3 (i.e., contamination <10)
* How did I know which columns to select?
`head 04_Binning/checkm2/${sample}/quality_report.tsv `

This will only take a few seconds to run so you can execute it directly in terminal with `bash 04_filter_checkm_table.sh`

Let's download these files to our laptops. In a new terminal window on your laptop:

```
scp 'username@jed.epfl.ch:/scratch/username/your_dataset/04_binning/checkm2/*_quality_report_filtered.tsv' .
```

You can open these tsv files in your preferred programme (excel, R...).

**Q: How many bins good quality bins do you have?**

If you want to look at the bins that did not pass quality filtering (i.e., <50% completeness and >10% contamination) they are still in the `04_Binning/checkm2/${sample}/quality_report.tsv` file.

## CheckM

This step is optional, CheckM and CheckM2 do the same thing, CheckM2 just does it a bit faster and a bit better! I have included it as CheckM is widely used so it's good to be aware of it.

### 5. Running CheckM

```
mkdir 04_Binning/checkm1
```

```
nano 04_checkm1.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

#CheckM on Das Tool bins

checkm lineage_wf -f 04_Binning/checkm1/${sample}/${sample}_checkm.tsv \
04_Binning/dastool/${sample}/${sample}_DASTool_bins/ --tab_table -t $SLURM_CPUS_PER_TASK -x fa 04_Binning/checkm1/${sample}

done
```
CheckM works on a directory of genome bins in FASTA format. By default, CheckM assumes these files end with the extension ‘fna’, so we are changing it to `fa` with the `–x` flag. `--tab_table` specifies we want a tab delimited table with the results and `-t 20` sets the threads. The order of the rest of the command is:
`checkm lineage_wf -f <file_to_put_results_in.tsv> <bin folder> <output folder>`

CheckM takes a bit longer to run than CheckM2 so in your sbatch script you can use:

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --time 10:00:00
#SBATCH --mem=70G
```

The file **`04_Binning/checkm1/${sample}_checkm.tsv`** contains our results. We can select the bins that pass our minimum completeness (50%) and maximum contamination (10%) levels with the same `awk` command adjusted for the appropriate columns (the columns of the output file from CheckM are in a different order to CheckM2).

```
nano 04_filter_checkm1_table.sh
```

```
#!/bin/bash

for sample in $(cat samples.txt)

do

awk -F "\t" '$12 >=50 && $13 <=10' 04_Binning/checkm1/${sample}_checkm.tsv > 04_Binning/checkm1/${sample}/${sample}_checkm_filtered.tsv

done
```
`awk` commands are super quick, so you can execute directly in terminal with `bash 04_filter_checkm1_table.sh`


And download a copy of these files to your laptop with the `scp` command.
```
scp 'username@jed.epfl.ch:/scratch/username/your_dataset/04_binning/checkm1/*_checkm_filtered.tsv' .
```

**Q: Do you see any differences between the CheckM2 and CheckM results?**
