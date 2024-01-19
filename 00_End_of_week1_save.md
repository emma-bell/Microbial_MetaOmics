# Making a copy of our data at the end of Week 1
Our scratch directories will be deleted in the next two weeks so let's make a copy of the important folders in our home directories:

```
cd /home/USERNAME
```

### In your home directory make a folder to hold your data
```
mkdir YOURDATASET (e.g., Subsurface_KR46)
cd YOURDATASET
```

### **`01_ReadQC`**: we'll keep the quality filtered reads and the final fastqc report

```
mkdir 01_ReadQC
cp -r /scratch/ebell/Subsurface_KR46/01_ReadQC/fastp_reads 01_ReadQC 
cp -r /scratch/ebell/Subsurface_KR46/01_ReadQC/fastqc_pass2 01_ReadQC
```

### **`02_Assembly`**: we'll keep the final assembled contig files

```
mkdir 02_Assembly
cp -r /scratch/ebell/Subsurface_KR46/02_Assembly/filtered_contigs 02_Assembly
```

We won't keep the bam files from the **`03_Mapping`** as they're large and we won't need them again. But in "real life" don't delete them yet.

# We will keep the whole **`04_Binning`**, **`05_Annotation`** and **`06_Representative_MAGs`** directories
```
cp -r /scratch/ebell/Subsurface_KR46/04_Binning .
cp -r /scratch/ebell/Subsurface_KR46/05_Annotation .
cp -r /scratch/ebell/Subsurface_KR46/06_Representative_MAGs .
```

# Let's not forget our scripts!
```
mkdir scripts
cp /scratch/ebell/Subsurface_KR46/*.sh scripts/
```
