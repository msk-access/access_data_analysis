---
description: Step 1 -- intra-patient genotyping
---

# Compile Reads

The first step of the pipeline is to genotype all the variants of interest in the included samples \(this means plasma, buffy coat, DMP tumor, DMP normal, and donor samples\). Once we obtained the read counts at every loci of every sample, we then generate a table of VAFs and call status for each variant in all samples within a patient in the next step. 

## Usage

```text
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

## Default

Default options can be found [here](../setup/resources.md#compile-reads)

## What `compile_reads.R` does

### [For each patient](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L47)

* [Create a sample sheet](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L54-L69) -- similar to the one for `genotype-variants`

| Sample\_Barcode | duplex\_bams | simplex\_bams | standard\_bam | Sample\_Type | dmp\_patient\_id |
| :--- | :--- | :--- | :--- | :--- | :--- |
| plasma sample id | /duplex/bam | /simplex/bam | NA | duplex | P-xxxxxxx |
| buffy coat id | NA | NA | /unfiltered/bam | unfilterednormal | P-xxxxxxx |
| DMP Tumor ID | NA | NA | /DMP/bam | DMP\_Tumor | P-xxxxxxx |
| DMP Normal ID | NA | NA | /DMP/bam | DMP\_Normal | P-xxxxxxx |

* [Generate all variants of interests](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L70-L84)
  * DMP calls from cbio repo
  * ACCESS calls from SNV pipeline
* [Generate unique variants list](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L82-L83)
* [Tag hotspots on unique variants](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L86-L94)
* Genotype with [genotype-variants](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L95-L105)

### [Afterwards](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L111), for donor bams

* Obtain all variants genotyped in any patient, [generate a all unique list of variants](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L112-L121)
* Genotype with [genotype-variants](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/compile_reads.R#L122-L129)



