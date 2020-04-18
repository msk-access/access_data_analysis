# ACCESS Data Analysis

Scripts for downstream analysis plotting of pipeline output

---


## Steps for analysis

### 1. Separating copy number output into individual files

```{bash}
Rscript cna_divide_by_sample.R -i input/maf/file -o output/directory
```
