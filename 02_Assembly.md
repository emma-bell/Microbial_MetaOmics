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

After assembly, your contig file is **`contigs.fasta`**.
Always check the log files after an assembly to make sure no errors occurred. You can do this with `tail spades.log`.

Note: we now have FASTA files instead of FASTQ files. FASTA files only hold the nucleotide sequences and don't have the sequence quality values found in FASTQ. You can use the command `head your_fasta_file.fasta` to look at your fasta file.

## Assembly with MEGAHIT

We'll now assemble each of our samples using megahit. This means you will have to run the script below on each of your samples. Make sure you change the parts of the command that say *`your_sample`* to your sample name.

### 2. Run MEGAHIT

Megahit has two options *`meta-sensitive`* and *`meta-large`*, choose the appropriate one for your samples. 

```
megahit -1 01_ReadsQC/readsqc/your_sample_R1.fastq.gz -2 01_ReadsQC/readsqc/your_sample_R2.fastq.gz \
-o 02_Assembly/megahit/your_sample -t 20 --presets meta-sensitive OR meta-large --min-contig-len 500
```

After assembly, your contig file is **`final.contigs.fa`**.
Always check the log files after an assembly to make sure there are no errors occurred. You can do this with `tail megahit.log.txt`. MEGAHIT also provides some simple stats for the assembly at the end of the log file.

## Compare assemblies
Quast can be used to generate some simple metrics of your assembly (e.g., N50, L50). It can also be used to compare different assemblies. We'll use it to compare the assemblies from metaSPAdes and MEGAHIT.

### 3. Run QUAST

```
mkdir 02_Assembly/quast

quast.py 02_Assembly/metaspades/your_sample/contigs.fasta 02_Assembly/megahit/your_sample/final_contigs.fa \
-o 02_Assembly/quast/your_sample
```
**Q: Which assembly will you move forward with and why?**

## Some housekeeping before we move on to mapping

### 1. Check contig lengths
Let's first check the length of our shortest contigs. 

MEGAHIT contigs are a minimum of 500 bp as we defined this in our MEGAHIT run with the flag `--min-contig-len 500`. We can also check this in the MEGAHIT log which gives us some basic assembly statistics. 

For the metaSPAdes assemblies we can check the contig lengths using seqkit.

```
cat 02_Assembly/metaspades/your_sample/final.contigs.fa | seqkit seq | seqkit stats
```

**Q: What is the minimum contig length should we use?**

### 2. Make a directory for your size filtered contigs
```
mkdir 02_Assembly/trimmed_contigs
```

### 3. Filter out short contigs

We can use a script provided in the [BBTools suite](https://anaconda.org/bioconda/bbmap) in a for loop to remove contigs that are shorter than our chosen threshold. 

```
for sample in $(cat samples.txt)

do

reformat.sh in=02_Assembly/your_chosen_assembly/${sample}/contigs.fa \
out=02_Assembly/trimmed_contigs/${sample}_contigs.fa minlength=your_threshold

done
```

### 4. Add a prefix to contig headers
We'll also add a prefix to the contig headers to make it easier to keep track of our samples during the next steps. It's important that you don't use a space character, or a pipe (|) character, as different tools will treat them differently. Let's use underscores and keep the sequence IDs in your FASTA file as simple as possible.

To append your sample name to each header we'll use a `sed` command from inside your "trimmed_contigs" directory.

```
sed 's/>/>samplename_/' your_contigs.fasta | head
```

Does the header look how you expected it to turn out? If yes, we can remove the pipe to `head` and add the -i flag (in place) to make the changes inside file.
```
sed -i 's/>/>samplename_/' your_contigs.fasta
```
