# Taxonomic Annotation
We're going to use [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) to assign taxonomy to our genomes based on the Genome Database Taxonomy Database [(GTDB)](https://gtdb.ecogenomic.org).

Again, we'll make a new directory for the results.

`mkdir 05_Annotation`

## GTDB-Tk

We'll run the *`classify_wf`* to get the taxonomy of our bins.
```
gtdbtk classify_wf --genome_dir 04_Binning/dastool/your_sample/your_sample_DASTool_bins/ -x fa --out_dir 05_Annotation/gtdb/your_sample --mash_db 05_Annotation/gtdb/your_sample --prefix your_sample --cpus 20
```

When the run is finished, make sure to check the log file for errors.

`tail 05_Annotation/gtdb/KR46_May/gtdbtk.log`.

Taxonomic annotation of bacteria bins are in **`your_sample.bac120.summary.tsv`** and archaea bins are in **`your_sample.ar53.summary.tsv`**.

# Functional Annotation
We're going to perform functional annotation with [METABOLIC](https://github.com/AnantharamanLab/METABOLIC).