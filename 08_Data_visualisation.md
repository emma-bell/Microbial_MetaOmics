# Data visualisation
We've completed the bioinformatic analyses but we also want to make sense of some of the data we've generated. Below is an outline of the files we're going to use and where you should find them in your working directory:

```
05_Annotation/
├── gtdb
│   └── SAMPLE
│       ├── SAMPLE.ar53.summary.tsv
│       └── SAMPLE.bac120.summary.tsv
└── metabolic
│   └── SAMPLE
│       └── METABOLIC_result_each_spreadsheet
│       └──  METABOLIC_result.xlsx
06_Representative_MAGs/
├── coverm
│   ├── coverm_drep_mag_rel_abun.tsv
└── drep
    ├── checkm2_for_drep.csv
    └── out
        ├── data_tables
        │   ├── Cdb.csv
        │   ├── genomeInfo.csv
        │   ├── Wdb.csv
        └──  dereplicated_genomes
```

We're going to do most of the data visualisation with R Studio. So if you haven't already, create a copy of these files on your laptop. I find it useful to create the same file structure on my laptop to make it easier to keep track of things.

## Generate figures for presentations/papers
In R Studio we're going to be using the packages [Tidyverse](https://tidyverse.tidyverse.org) and [Complex Heatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html). Instructions for installation can be found at the previous links and on the [00_Requirements](00_Requirements.md) page.

## Make a phylogenomic tree with GToTree
We're going to make a phylogenomic tree using [GToTree](https://github.com/AstrobioMike/GToTree/wiki/example-usage). This tool is really easy to use and has excellent documentation so it's worth checking out the link for more information. Below we'll go through an example with some of the data we've been working on.

Let's log into SCITAS and make a directory for data visualisation in your working directory:

`mkdir 08_Visualisation`

In R Studio I'm looking at the taxonomy of my MAGs in the table we created called **`gtdb_tax`**.
Looking at the table I'm interested in the genus "SURF-10". In the data I'm using in this example the genomes came from deep groundwater in Finland. I follow subsurface research so I recognise that "SURF" probably means the genus name came from a genome from the Sanford Underground Research facility. I'd like to know a bit more and see if our SURF-10 MAGs are closely related.

From your dataset you can choose any MAG/taxonomy you want, at any taxonomic level.

We'll need to load an interactive session to use GToTree.

```
srun -A bioinformatics-meta-omics2 --pty apptainer shell --bind /scratch/YOUR_USERNAME/YOUR_DATASET:/data PATH_TO_IMAGE
```

Change into your data folder `cd /data`

With GtoTree, we can look for "SURF-10" genomes in GTDB using the following command:

```
gtt-get-accessions-from-GTDB -t SURF-10 --get-taxon-counts --GTDB-representatives-only
```

GToTree will give you (1) the total number of genomes sequenced and (2) how many representative genome entries there are. In your own research, if you're working with a large breadth of diversity with sequenced entries, you may just want to use the representative genomes. For this class, try to choose a group with relatively few entries. In this example, there are 27 genomes for the genus SURF-10.

E.g., of what you will see on your screen:
```
  Reading in the GTDB info table...
    Using GTDB v214.1: Released Jun 9th, 2023


    The rank 'genus' has 27 SURF-10 entries.

  In considering only GTDB representative genomes:

    The rank 'genus' has 8 SURF-10 representative genome entries.
```

Feel free to try a few different options until you find a group of genomes you'd like to make a phylogenomic tree with.
Once you're happy with your choice remove the `--get-taxon-counts` flag:

```
gtt-get-accessions-from-GTDB -t SURF-10
```

And GToTree will produce the following files which contain the genome accessions for your selected genomes, and their associated metadata:

```
The targeted NCBI accessions were written to:
    GTDB-SURF-10-genus-accs.txt

  A subset GTDB table of these targets was written to:
    GTDB-SURF-10-genus-metadata.tsv
```

Now we need to choose which HMM profile to run. You can look at all of the available HMM sets with:

```
gtt-hmms
```

I'm going to choose "Bacteria" for my SURF10 genomes which have the taxonomy:
`d__Bacteria;p__Desulfobacterota;c__Desulfarculia;o__Desulfarculales;f__Desulfarculaceae;g__SURF-10;s__SURF-10 sp018902755`. You want to choose the most specific HMM set that you can.

Next we will create a file telling GToTree which of our genomes we want to input.
In R studio, I'm filtering the rows of the **`MAG_rep_summary`** by the taxonomy classification I'm interested in to recover the file name of my SURF-10 genome:

```
MAG_rep_summary %>% 
  filter(classification == "d__Bacteria;p__Desulfobacterota;c__Desulfarculia;o__Desulfarculales;f__Desulfarculaceae;g__SURF-10;s__SURF-10 sp018902755") %>%
  %>% select(genome)
```

For me, this returns the genome entry `KR46_May_metabat.bin.4`.

Back in terminal, create a text file:

```
nano MAG_fasta.txt
```

And put the path to your genome:

```
/scratch/USERNAME/06_Representative_MAGs/MAGs/YOURGENOME.fa
```

We also need an outgroup to be able to root the tree. I'm using another genome from my dataset. It doesn't matter which genome you choose but consider **(1)** where possible use a genome that the chosen HMM set will apply to. If too many SCGs are missing from the genome GToTree will exclude it from the analysis **(2)** it should be different enough from your genomes of interest - I've chosen a different phylum, within the domain bacteria.

Add the path of your outgroup to the **`MAG_fasta.txt`** file.

Another cool thing about GToTree is that you can look for genes of interest using an additional HMM search. In this example, let's say I would like to know whether all sequenced members of the SURF-10 genus have the genes _dsrA_ and _dsrB_ that make up the dissimilatory sulfate reductase enzyme. The [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) for those genes are `K11180` and `K11181`. 

You can choose a gene and `KO` ID from your annotation results file **`METABOLIC_result.tsv`**. You can also view different [pathway maps](https://www.genome.jp/kegg/pathway.html) on the KEGG PATHWAY Database website.

When you've chosen a gene (or genes) you'd like to look for within your set of genomes, create a text file and enter the KO ID. If you have more than one ID, enter one per line.

```
nano kofams.txt
```

E.g.,
```
K11180
K11181
```

Now we're ready to run GToTree:

```
GToTree -f MAG_fasta.txt -a GTDB-SURF-10-genus-accs.txt -H Bacteria  -K kofams.txt -o GToTree_output -D -j 6
```

* `-f` : text file containing the location of your MAGs
* `-a` : the accessions of the SURF-10 genomes
* `-H` : the HMM set we want to use
* `-K` : KO ids of genes we're interested in
* `-o` : name of your output folder
* `-D` : flag to use GTDB taxonomy
* `-j` : jobs to run in parallel

With 27 genomes and 6 threads, this took 2 minutes to complete.

The following files will be produced.

```
GToTree_output/
├── Aligned_SCGs_mod_names.faa
├── citations.txt
├── Genes_with_no_hits_after_length_filter.txt
├── Genomes_summary_info.tsv
├── GToTree_output.tre
├── gtotree-runlog.txt
├── KO_search_results
├── run_files
└── SCG_hit_counts.tsv
```

Download **`/scratch/USERNAME/Dataset/08_Visualisation/GToTree_output`** to your laptop with `scp -r`. Then you can use [iTOL](https://itol.embl.de) to view the tree you've created.

On iTOL, you can view a tree without creating an account by going to the `Upload` tab. If you want to be able save your tree, you can create a free account.

Drag **`GToTree_output.tre`** into iTOL tree upload. You should now be able to see your tree. The first thing you want to do is reroot the tree at your outgroup genome. The output genome is often easy to spot as it will have a very long branch. Otherwise you can also search the tree nodes by name by clicking on the magnifer on the left of the screen.

Click on the outgroup branch, select `Tree stucture` followed by `Re-root the tree here`.

GToTree produces some iTOL-ready annotation files in **`KO_search_results/iToL_files/`**. You can drag and drop those files directly on the tree. In this case, the files will colour the branches of genomes with the genes _dsrA_ or _dsrB_ blue. There are many more options and ways to annotate in iTOL either interactively or by uploading different annotation files.

For example, we can highlight our genomes in the tree by clicking on them and changing the `Label`. We can colour different clades by selecting `Coloured ranges`.

Have a play around and check out the iTOL [help pages](https://itol.embl.de/help.cgi).