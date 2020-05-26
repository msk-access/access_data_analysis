---
description: Step 4 -- generating final CNA call set
---

# CNA Processing

This step generates a final CNA call set for plotting. This consists of: 

* Calls passing de novo CNA calling threshold 
  * Significant adjusted p value \( &lt;=  0.05\)
  * Significant fold change \( &gt; 1.5 or &lt; -1.5\)
* Calls that can be rescued based on prior knowledge from IMPACT samples
  * Significant adjusted p value \( &lt;= 0.05\)
  * Lowered threshold for fold change \( &lt; 1.2 or &lt; -1.2\)

## Usage

```text
Rscript R/CNA_processing.R -h                                       
usage: R/CNA_processing.R [-h] [-m MASTERREF] [-o RESULTSDIR] [-dmp DMPDIR]

optional arguments:
  -h, --help            show this help message and exit
  -m MASTERREF, --masterref MASTERREF
                        File path to master reference file
  -o RESULTSDIR, --resultsdir RESULTSDIR
                        Output directory
  -dmp DMPDIR, --dmpdir DMPDIR
                        Directory of clinical DMP IMPACT repository [default]
```

## What `CNA_processing.R` does

1. [Process DMP CNAs](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/CNA_processing.R#L19-L25) 
2. [Process ACCESS CNAs](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/CNA_processing.R#L29-L33) 
3. [CNA Calling](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/CNA_processing.R#L36-L44) \(de novo and rescue\)

Format of the final call set:

| Tumor\_Sample\_Barcode | cmo\_patient\_id | Hugo\_Symbol | p.adj | fc | CNA\_tumor | CNA | dmp\_patient\_id |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  |  |  |  |



