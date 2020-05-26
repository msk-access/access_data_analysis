---
description: Step 2 -- filtering
---

# Filter Calls

The second step takes all the genotypes generated from the first step and organized into a patient level variants table with VAFs and call status for each variant of each sample.  

Each call is subjected to:

1. Read depth filter \(hotspot vs non-hotspot\)
2. Systematic artifact filter
3. Germline filters 
   1. If any normal exist -- \(buffy coat and DMP normal\) 2:1 rule
   2. If not -- exac freq &lt; 0.01% and VAF &lt; 30%
4. CH tag

## Usage

```text
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

## What `filter_calls.R` does

[Generate a reference of systematic artifacts](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L35-L39) -- any call with occurrence in more than or equal to 2 donor samples \(occurrence defined as more than or equal to 2 duplex reads\)

### [For each patient](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L43)

1. [Read in sample sheets](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L48-L50) -- reference for downstream analysis
2. Generate a [preliminary patient level variants table](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L52-L77) 
3. Read in and merging in [hotspots, DMP signed out calls and occurrence in donor samples](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L78-L97)
4. Call status annotation
   1. All call passing read depth/genotype filter annotated as[ 'Called' or 'Genotyped'](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L112-L126)
   2. Any call not satisfying germline filters are [overwritten](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L127-L145) with 'Not Called'
      1. Calls with zero coverage in plasma sample also annotated as 'Not Covered'
5. Final processing
   1. [Combining](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L150-L160) duplex and simplex read counts
   2. [CH tags](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/filter_calls.R#L161-L174)
6. Write out table

### Example of the patient level table:

| Hugo\_Symbol | Start\_position | Variant\_Classification | Other variant descriptions | ... | C-xxxxxx-L001-d\_\_\_duplex.called | C-xxxxxx-L001-d\_\_\_duplex.total | C-xxxxxx-L002-d\_\_\_duplex.called | C-xxxxxx-L001-d\_\_\_duplex.total | C-xxxxxx-N001-d\_\_\_unfilterednormal | P-xxxxxxx-T01-IM6\_\_\_DMP\_Tumor | P-xxxxxxx-T01-IM6\_\_\_DMP\_Normal |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| KRAS | xxxxxx | Missense Mutation | ... | ... | Called | 15/1500\(0.01\) | Not Called | 0/1800\(0\) | 0/200\(0\) | 200/800\(0.25\) | 1/700\(0.001\) |

