---
description: Step 3 -- incorporating SVs into patient table
---

# SV Incorporation

The third step takes all the SV variants from all samples within each patient and present them in the same format as SNVs and incorporate SVs in the patient level table. 

## Usage

```text
Rscript R/SV_incorporation.R -h                                     
usage: R/SV_incorporation.R [-h] [-m MASTERREF] [-o RESULTSDIR] [-dmp DMPDIR]
                            [-c CRITERIA]

optional arguments:
  -h, --help            show this help message and exit
  -m MASTERREF, --masterref MASTERREF
                        File path to master reference file
  -o RESULTSDIR, --resultsdir RESULTSDIR
                        Output directory
  -dmp DMPDIR, --dmpdir DMPDIR
                        Directory of clinical DMP IMPACT repository [default]
  -genes GENELIST, --genelist GENELIST
                        File path to genes covered by ACCESS [default]
  -c CRITERIA, --criteria CRITERIA
                        Calling criteria [default]
```

## What `SV_incorporation.R` does

[Gets DMP signed out SV calls](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L20-L24)

### [For each patient](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L30)

1. [Process plasma sample SVs ](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L34-L39)
   1. Only SVs implicating any ACCESS SV calling key genes are retained
2. [Process DMP SVs](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L41-L55) to similar format to ACCESS SV output
3. [Row-bind plasma and DMP SVs](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L56-L77) and make call level info \(similar to SNVs\)
4. [Annotate](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/SV_incorporation.R#L80-L100) call status for each call of each sample
   1. Not Called
   2. Not Covered -- none of the genes in key genes
   3. Called
5. Read in SNV table, row-bind with SV table, write out table

