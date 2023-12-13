# Assembly of QC reads

Now we have quality filtered our reads, let's move on to assembly. We are going to use both [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](https://github.com/ablab/spades) and compare the results using [QUAST](https://quast.sourceforge.net).

First, let's make our output directories
```
mkdir -p 02_Assembly/megahit
mkdir -p 02_Assembly/metaspades
```

## Assembly with metaSPAdes

We're going to do assembly on each of our samples individually. This means you will have to run the script below on each of your samples. Make sure you change the parts of the command that say *`your_sample`* to your sample name.

### 1. Run metaSPAdes
```
metaspades.py -1 01_ReadQC/readsqc/your_sample_R1.fastq.gz -2 01_ReadQC/readsqc/your_sample_R2.fastq.gz -o 02_Assembly/metaspades/your_sample \
--threads 20 -m 125
```
Note: Threads and memory should be adjusted for your computing resources. In this example, 1 node has 20 cores and 128 GB RAM so we are running 20 threads using a maximum of 125 GB RAM. Larger samples will need more memory than this.

After assembly, your contig file is **`02_Assembly/metaspades/your_sample/contigs.fasta`**.
Always check the log files after an assembly to make sure no errors occurred. You can do this with `tail 02_Assembly/metaspades/your_sample/spades.log`.

Note: we now have FASTA files instead of FASTQ files. FASTA files only hold the nucleotide sequences and don't have the sequence quality values found in FASTQ. You can use the command `head 02_Assembly/metaspades/your_sample/your_fasta_file.fasta` to look at your fasta file.

## Assembly with MEGAHIT

We'll now assemble each of our samples using megahit. This means you will have to run the script below on each of your samples. Make sure you change the parts of the command that say *`your_sample`* to your sample name.

### 2. Run MEGAHIT

Megahit has two options *`meta-sensitive`* and *`meta-large`*, choose the appropriate one for your samples. 

```
megahit -1 01_ReadsQC/readsqc/your_sample_R1.fastq.gz -2 01_ReadsQC/readsqc/your_sample_R2.fastq.gz \
-o 02_Assembly/megahit/your_sample -t 20 --presets meta-sensitive OR meta-large --min-contig-len 500
```

After assembly, your contig file is **`02_Assembly/megahit/your_sample/final.contigs.fa`**.
Always check the log files after an assembly to make sure there are no errors occurred. You can do this with `tail 02_Assembly/megahit/your_sample/megahit.log.txt`. MEGAHIT also provides some simple stats for the assembly at the end of the log file.

## Compare assemblies
Quast can be used to generate some simple metrics of your assembly (e.g., N50, L50). It can also be used to compare different assemblies. We'll use it to compare the assemblies from metaSPAdes and MEGAHIT.

### 3. Run QUAST

If you had ran both megahit and metaspades yourself you would use the below script to compare the output of each:

```
quast.py 02_Assembly/metaspades/your_sample/contigs.fasta 02_Assembly/megahit/your_sample/final_contigs.fa \
-o 02_Assembly/quast/your_sample
```

As each member of a pair ran megahit or metaspades, we can produce a QUAST output for all of our megahit samples and all of our metaspades samples. You can then compare the output with your partner.

### 4. For everyone:

First make an output directory:
```
mkdir 02_Assembly/quast
```
Then make a script with **`nano`** called **`02_quast.sh`**

### 4a. If you ran metaspades:
```
#!/bin/bash

cd /data

quast.py 02_Assembly/metaspades/*/contigs.fasta -o 02_Assembly/quast/quast_output
```

### 4b. If you ran megahit:
```
#!/bin/bash

cd /data

quast.py 02_Assembly/megahit/*/final.contigs.fa -o 02_Assembly/quast/quast_output
```

### 4c. For everyone:

Make an sbatch script for quast - you'll need to change the image file to **`quast.sif`** and your script to **`02_quast.sh`**. Quast runs quickly so you can set the following parameters:

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:15:00
#SBATCH --mem=7G
#SBATCH --account=bioinformatics-meta-omics1
#SBATCH --reservation=bioinformatics_meta-omics1
#SBATCH --mail-user=your_email_if_you_want_to_be_notified_when_the_job_is_complete
#SBATCH --mail-type=ALL
```

Note: Each node has 500GB mem and 72 CPUs. That means each CPU has ~7GB mem (500/72 = 6.94). So we should try to match the memory to the cpus we are requesting, e.g., if we want to use 20 threads, 20x7 would be 140GB mem. Conversely, if you know you need 250 GB mem, that is half of the total memory on a node (500GB) so you should request half of the cpus (72/2 = 36).

### 4d. Download the quast results

Your results are in **`02_Assembly/quast`**. You can download this folder and view the results interactively by clicking on the html file.

In a new terminal window, change into the directory (`cd`) on your laptop where you want to store the results. Then we can copy the files from the server:

```
scp -r username@jed.epfl.ch:/scratch/username/dataset/02_Assembly/quast .
```
* The `-r` flag is to copy a directory
* The `.` means "current directory" i.e., download the files into the directory on your laptop you are running the command from

**Q: Which assembly was best out of the samples you ran? Which assembler was better?**

## Some housekeeping before we move on to mapping

### 1. Check contig lengths
Let's first check the length of our shortest contigs. 

MEGAHIT contigs are a minimum of 500 bp as we defined this in our MEGAHIT run with the flag `--min-contig-len 500`. We can also check this in the MEGAHIT log which gives us some basic assembly statistics. 

For the metaSPAdes assemblies (or if you want to check the megahit another way) we can check the contig lengths using seqkit.

Make a seqkit script:
```
nano 02_seqkit.sh
```
In your script:
```
#!/bin/bash
cd /data

cat 02_Assembly/metaspades/your_sample/contigs.fa | seqkit seq | seqkit stats > 02_Assembly/metaspades/your_sample/stats.txt
```
* `cat` is reading the contig file. The pipe `|` then passes the output into the seqkit command which will reads the sequences and give you some stats.



If you ran megahit make sure to change the path to `02_Assembly/megahit/final.contigs.fa`

Seqkit also runs very quickly, so you can use the same SBATCH parameters you just used for quast. Make sure to change the image to `seqkit.sif` and the script to `02_seqkit.sh` in your sbatch script.

**Q: What is the minimum contig length should we use?**

### 2. Make a directory for your size filtered contigs

If you are happy with your minimum contig length you don't need to filter your contigs. If you have short contigs (<500 bp) we'll filter them out. A longer minimum contig length is better, but it is also a trade off with how many long contigs you have in your assembly.

Let's make a directory for our size filtered contigs:

```
mkdir 02_Assembly/filtered_contigs
```

### 3. Filter out short contigs

We can use a script provided in the [BBTools suite](https://anaconda.org/bioconda/bbmap) in a for loop to remove contigs that are shorter than our chosen threshold. 

```
nano 02_bbtools.sh
```

```
#!/bin/bash

cd /data

for sample in $(cat samples.txt)

do

reformat.sh in=02_Assembly/your_assembly_method/${sample}/your_contigs_file.fa \
out=02_Assembly/filtered_contigs/${sample}_contigs.fa minlength=your_chosen_threshold

done
```
Make an sbatch script with the image **`bbmap.sif`** and the name of your script **`02_bbtools.sh`**

We can run this for an hour but it will probably run much faster:

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --time 01:00:00
#SBATCH --mem=28G
```

### 4. Add a prefix to contig headers
We'll also add a prefix to the contig headers to make it easier to keep track of our samples during the next steps. It's important that you don't use a space character, or a pipe (|) character, as different tools will treat them differently. Let's use underscores and keep the sequence IDs in your FASTA file as simple as possible.

To append your sample name to each header we'll use a `sed` command from inside your "filtered_contigs" directory.

```
cd 02_Assembly/filtered_contigs
```
Then run:
```
sed 's/>/>samplename_/' sample_name.fa | head
```
* sed is performing a substitution (s/) in 'sample_name.fasta'.
* `>` matches the ">" character at the beginning of each sequence header.
* `>samplename_` is what we are replacing the matched ">" with.
* `|` The pipe redirects the output of the sed command to the input of the next command `head`
* `head` displays the first few lines of the modified content. We do this because 'sample_name.fasta' is a large file and we just want to check the first line

Does the header look how you expected it to turn out? If yes, we can remove the pipe to `head` and add the `-i` flag which means edit "in place" to make the changes inside file.
```
sed -i 's/>/>samplename_/' sample_name.fa
```

Now if you look at your contig file, you should see you sample name in the contig header:
```
head sample_name.fa
```
e.g., ">KR46_Sept_NODE_1_length_1038990_cov_118.723418".

You should do this for each of your samples.

Note: if you're working directory is getting full of slurm-123456789.out files, feel free to delete them
```
rm slurm-*.out
```
* The `*` is a wildcard that will select groups based on a common pattern