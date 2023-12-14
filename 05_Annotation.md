# Taxonomic Annotation
We're going to use [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) to assign taxonomy to our genomes based on the Genome Database Taxonomy Database [(GTDB)](https://gtdb.ecogenomic.org).

Again, we'll make a new directory for the results.

`mkdir 05_Annotation`

## GTDB-Tk

We'll run the *`classify_wf`* to get the taxonomy of our bins.

```
nano 05_gtdb.sh
```

```
#!/bin/bash
cd data/

for sample in $(cat samples.txt)

do

gtdbtk classify_wf --genome_dir 04_Binning/dastool/${sample}/${sample}_DASTool_bins -x fa --out_dir 05_Annotation/gtdb/${sample} --mash_db 05_Annotation/gtdb/${sample} --prefix ${sample} --cpus 20

done
```

Then create a sbatch script:

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 10:00:00
#SBATCH --mem=140G
```

When the run is finished, make sure to check the log file for errors.

`tail 05_Annotation/gtdb/${sample}/gtdbtk.log`.

Taxonomic annotation of bacteria bins are in **`05_Annotation/gtdb/${sample}/${sample}.bac120.summary.tsv`** and archaea bins are in **`05_Annotation/gtdb/${sample}/${sample}.ar53.summary.tsv`**.

Download these files to your laptop:

```
scp 'username@jed.epfl.ch:/scratch/username/your_dataset/05_Annotation/gtdb/*summary.tsv' .
```
You can view this files in excel or R. What taxa do you have in your samples??

# Functional Annotation
We also want to know what the microbes we have might be doing. We're going to perform functional annotation with [METABOLIC](https://github.com/AnantharamanLab/METABOLIC).

We can run METABOLIC starting using **amino acid** sequences of our bins. Amino acid files were generated during our checkm2 run so we'll use those, but if you want to convert from a **nucleotide** file to an **amino acid** file yourself in the future you can use a tool called [Prodigal](https://github.com/hyattpd/Prodigal/tree/GoogleImport).

Our amino acid files are here:
**`04_Binning/checkm2/sample/protein_files/`**. They now end with `.faa` meaning it is a fasta file that contains amino acid sequences. You can use the head command to take a look.

```
>KR46_June_NODE_64_length_85318_cov_103.948078_2 # 497 # 1834 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.389
MIFESLSNKLQGALNKLKGKGKLSEKDIDEAMRDIRLSLLEADVNFKVVKDFVKSVKERS
MGSEVMESLTPGQQVVKIVNEEMIKILGEKESKIQFSSTETTVIMMCGLQGAGKTTTAGK
LALRFKKQNKRPLLVACDIYRPAAIKQLEVVGKSVEVPVFKIEGETNPVIIADKALKEAR
KNGNDVLIIDTAGRLHIDEKLMEELIQIRNKVKPSEVLLVLDAMTGQDAVKIAESFNQNM
EITGLILTKVDGDARGGAAISIRAVTSKPIKFVTTGEKMADLEAFHPDRMASRILGMGDL
LSLIEKAQESFDSKKVKEMEEKLRGQGFTLDDFLDQMEQMKSLGPLDQLLEMIPGANSKQ
```

Make the directory for our annotation:
```
mkdir 05_Annotation/metabolic
```

And make our script:
```
nano 05_metabolic.sh
```

```
#!/bin/bash
cd /data

for sample in $(cat samples.txt)

do

perl METABOLIC-G.pl -in 04_Binning/checkm2/${sample}/protein_files/ -o 05_Annotation/metabolic/${sample}

done
```

```
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --time 10:00:00
#SBATCH --mem=140G
```

METABOLIC will take some time to run, so we can move on to the next steps whilst we are waiting for it to complete.

