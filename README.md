# ACCESS Data Analysis

---

Scripts for downstream analysis plotting of pipeline output


## Setting up testing directory on JUNO

---


### Softlink and copy bam files, pipeline outputs into one folder

```{bash}
bash softlinking_bams
```

### Processing master file to specify certain file paths

```{bash}
Rscript process_master_file.R 
```


## Steps for analysis

---

### 1. Separating copy number output into individual files

```{bash}
Rscript cna_divide_by_sample.R -i input/maf/file -o output/directory
```
