# ACCESS Data Analysis

Scripts for downstream analysis plotting of pipeline output

---


## Steps for analysis

### 1. Separating copy number output into individual files

```{bash}
Rscript cna_divide_by_sample.R -i input/maf/file -o output/directory
```

### 2. Compiling genotypes for SNVs across all samples within each patients

```{bash}
Rscript R/compile_reads.R -h                                        
usage: R/compile_reads.R [-h] [-m MASTERREF] [-o RESULTSDIR]
                         [-pb POOLEDBAMDIR] [-fa FASTAPATH]
                         [-gt GENOTYPERPATH] [-dmp DMPDIR] [-mb MIRRORBAMDIR]
                         [-dmpk DMPKEYPATH]

optional arguments:
  -h, --help            show this help message and exit
  -m MASTERREF, --masterref MASTERREF
                        File path to master reference file
  -o RESULTSDIR, --resultsdir RESULTSDIR
                        Output directory
  -pb POOLEDBAMDIR, --pooledbamdir POOLEDBAMDIR
                        Directory for all pooled bams [default]
  -fa FASTAPATH, --fastapath FASTAPATH
                        Reference fasta path [default]
  -gt GENOTYPERPATH, --genotyperpath GENOTYPERPATH
                        Genotyper executable path [default]
  -dmp DMPDIR, --dmpdir DMPDIR
                        Directory of clinical DMP IMPACT repository [default]
  -mb MIRRORBAMDIR, --mirrorbamdir MIRRORBAMDIR
                        Mirror BAM file directory [default]
  -dmpk DMPKEYPATH, --dmpkeypath DMPKEYPATH
                        DMP mirror BAM key file [default]

```
Example command: 

```{bash}
Rscript R/compile_reads.R -m /juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv -o /juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020
```


### 3. Generating patient level table from the genotyping information

```{bash}
Rscript R/filter_calls.R -h                                         
usage: R/filter_calls.R [-h] [-m MASTERREF] [-o RESULTSDIR] [-dmpk DMPKEYPATH]
                        [-ch CHLIST] [-c CRITERIA]

optional arguments:
  -h, --help            show this help message and exit
  -m MASTERREF, --masterref MASTERREF
                        File path to master reference file
  -o RESULTSDIR, --resultsdir RESULTSDIR
                        Output directory
  -dmpk DMPKEYPATH, --dmpkeypath DMPKEYPATH
                        DMP mirror BAM key file [default]
  -ch CHLIST, --chlist CHLIST
                        List of signed out CH calls [default]
  -c CRITERIA, --criteria CRITERIA
                        Calling criteria [default]
```
Example commands:

```{bash}
Rscript R/filter_calls.R -m /juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv -o /juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020
```