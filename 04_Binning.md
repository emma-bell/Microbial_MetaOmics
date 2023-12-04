# Binning

Now we have our contigs (.fasta) and our alignments (.sorted.bam) we can move onto binning.
We're going to use two assemblers, [CONCOCT](https://github.com/BinPro/CONCOCT) and [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/).

1. Let's start by making the directories we'll need for the output.
```
mkdir -p 04_Binning/concoct
mkdir 04_Binning/metabat
```

## CONCOCT
CONCOCT usually takes longer to run than MetaBAT2, so let's start with that one.
You'll also need a directory for each of your samples inside the **`04_Binning/concoct`** directory.

2. Running CONCOCT
```
#First the contigs are cut up into 10000 bp

cut_up_fasta.py 02_Assembly/trimmed_contigs/your_sample_contigs.fa -c 10000 -o 0 --merge_last -b 04_Binning/concoct/your_sample/your_sample_contigs_10K.bed > 04_Binning/concoct/your_sample/your_sample_contigs_10K.fa

#Then a coverage table is made using our alignments

concoct_coverage_table.py 04_Binning/concoct/your_sample/your_sample_contigs_10K.bed 03_Mapping/your_sample/*.sorted.bam > 04_Binning/concoct/your_sample/your_sample_coverage_table.tsv

#Now we're ready to run CONCOCT

concoct --composition_file 04_Binning/concoct/your_sample/your_sample_contigs_10K.fa --coverage_file 04_Binning/concoct/your_sample/your_sample_coverage_table.tsv -b 04_Binning/concoct/your_sample/your_sample

#Merge subcontig clustering into original contig clustering

merge_cutup_clustering.py 04_Binning/concoct/your_sample/your_sample_clustering_gt1000.csv > 04_Binning/concoct/your_sample/your_sample_clustering_merged.csv

#Extract the bins as individual fasta

mkdir 04_Binning/concoct/your_sample/bins_fasta
extract_fasta_bins.py 02_Assembly/trimmed_contigs/your_sample_contigs.fa 04_Binning/concoct/your_sample/your_sample_clustering_merged.csv --output_path 04_Binning/concoct/your_sample/bins_fasta
```

Congratulations, you have bins! You can find them in the **`bins_fasta`** directory. Repeat the above process for each of your samples.

Before we move onto binning with MetaBAT2, let's add the sample name and binning method to the names of our bin files.

```
#Move into the directory of the bin.fasta files you want to rename

cd 04_Binning/concoct/your_sample/bins_fasta/

#This for loop will show us what the renamed files will look like, without making the change

for i in *.fa; do echo your_sample_concoct_bin.$i; done

#Did everything look right? If so, we can change the command to commit to the change

for i in *.fa; do mv $i your_sample_concoct_bin.$i; done
```
Perform the above renaming step on each of your samples.

## MetaBAT2
We're now ready to do binning with MetaBAT2. You should repeat this step for each of your samples.

3. Running MetaBAT2

```
#First make a coverage file using our alignments

jgi_summarize_bam_contig_depths --outputDepth 04_Binning/metabat/your_sample/your_sample_depth.txt 03_Mapping/your_sample/*.sorted.bam

#Now we're ready to run MetaBAT

metabat2 -m 1500 -i 02_Assembly/trimmed_contigs/your_sample_contigs.fa -a 04_Binning/metabat/your_sample/your_sample_depth.txt -o 04_Binning/metabat/your_sample/bins_fasta/your_sample_metabat.bin

#MetaBAT also has an option to summarise the bin depths

aggregateBinDepths.pl 04_Binning/metabat/your_sample/your_sample_depth.txt 04_Binning/metabat/your_sample/bins_fasta/*.fa > 04_Binning/metabat/your_sample/your_sample_bindepth.txt
```

Congratulations, you now have more bins! You can find them in the **`bins_fasta`** directory. We added the sample name and binning method to our file names within the `metabat2` command so we don't need to rename our bin files like we did after running CONCOCT.

We're now ready to move on to Bin QC.

# Bin QC
We want to know how complete our bins are and whether they contain any contamination. To do this we're going to use [CheckM2](https://github.com/chklovski/CheckM2). We can optionally also use [CheckM](https://github.com/Ecogenomics/CheckM) and compare the differences.

## CheckM2

We'll use the "predict" module of CheckM2 to get the completeness and contamination of both our CONCOCT and MetaBAT2 bins.

4. Running CheckM2
```
#CheckM2 on CONCOCT bins

~/checkm2/bin/checkm2 predict --threads 6 --input 04_Binning/concoct/your_sample/bins_fasta/ -x fa --output-directory 04_Binning/concoct/your_sample/checkm2 \
--database_path CheckM2_database/uniref100.KO.1.dmnd

#CheckM2 on MetaBAT2 bins

~/checkm2/bin/checkm2 predict --threads 6 --input 04_Binning/metabat/your_sample/bins_fasta/ -x fa --output-directory 04_Binning/metabat/your_sample/checkm2 \
--database_path CheckM2_database/uniref100.KO.1.dmnd
```

The file **`checkm2/quality_report.tsv`** contains our results. You should have one for each sample and each binning method. We can select the bins that pass our minimum completeness (50%) and maximum contamination (10%) thresholds with an `awk` command.

```
#Run this awk command from within your checkm2 output directory

awk -F "\t" '$2 >=50 && $3 <=10' quality_report.tsv > quality_report_filtered.tsv
```

Repeat the above filtering for each of your CheckM2 outputs. The file **`quality_report_filtered`** now contains a list of our bins that pass quality filtering. There should be two for each of your samples.

**Q: How many bins from CONCOCT compared to MetaBAT2 passed the filtering step?**

## CheckM

We can optionally use CheckM with the following instructions.

5. Running CheckM
```
#CheckM on CONCOCT bins

checkm lineage_wf -f 04_Binning/concoct/your_sample/checkm/your_sample_checkm.tsv \
04_Binning/concoct/your_sample/bins_fasta --tab_table -t 20 -x fa 04_Binning/concoct/your_sample/checkm/

#CheckM on MetaBAT2 bins

checkm lineage_wf -f 04_Binning/metabat/your_sample/checkm/your_sample_checkm.tsv \
04_Binning/metabat/your_sample/bins_fasta --tab_table -t 20 -x fa 04_Binning/metabat/your_sample/checkm/
```

The file **`checkm/your_sample_checkm.tsv`** contains our results. We can select the bins that pass our minimum completeness (50%) and maximum contamination (10%) levels with the same `awk` command adjusted for the appropriate columns.

```
#Run this awk command from within your checkm output directory

awk -F "\t" '$12 >=50 && $13 <=10' your_sample_checkm.tsv > your_sample_checkm_filtered.tsv
```

Repeat the above filtering for each of your CheckM outputs. The file **`your_sample_checkm_filtered.tsv`** now contains a list of our bins that pass quality filtering. There should be two for each of your samples.

**Q: What is the difference between CheckM2 and CheckM?**

## Bin refinement with DAS Tool

6. Performing bin refinement 

Before we can run [DAS Tool](https://github.com/cmks/DAS_Tool) we need to prepare the input files. Let's start by making a directory for the results.

`mkdir 04_Binning/dastool`

The input files we need are a list of contig headers and their associated bin. We can get this list using a script from the DAS Tool package.

```
#To get the CONCOCT contig input file

Fasta_to_Contig2Bin.sh -e fa -i 04_Binning/concoct/your_sample/bins_fasta/ > 04_Binning/dastool/your_sample_concoct.contigs2bin.tsv

#To get the MetaBAT2 contig input file

Fasta_to_Contig2Bin.sh -e fa -i 04_Binning/metabat/your_sample/bins_fasta/ > 04_Binning/dastool/your_sample_metabat.contigs2bin.tsv
```

Now we're ready to run DAS Tool on each of your samples.

```
DAS_Tool -i 04_Binning/dastool/your_sample_concoct.contigs2bin.tsv,04_Binning/dastool/your_sample_metabat.contigs2bin.tsv -l concoct,metabat -c 02_Assembly/trimmed_contigs/your_sample_contigs.fa -o 04_Binning/dastool/your_sample/your_sample --write_bins --threads 6
```

The file **`DASTool_summary.tsv`** summarises the results. If you have any bins ending with *`_sub`* in the **`DASTool_bins`** directory, these bins have been refined by DAS Tool. You can run CheckM2 again to see if the *`_sub`* bins have improved.
