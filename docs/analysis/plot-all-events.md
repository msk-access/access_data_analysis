---
description: Step 5 -- ex-graph-aganza
---

# Plot All Events

The final step takes the processed data from the previous steps and plots the genomic alterations over all samples of each patient. SNV and SV events are plotted out by VAFs over timepoints on the top panel and CNAs are plotted by fold-change\(fc\) on the bottom panel. The script currently supports the plasma sample only. 

{% hint style="danger" %}
Currently, the script will filter for 'High' call confidence variants \(This might need to be taken out for a fully automated pipeline\)
{% endhint %}

## Usage

```text
Rscript R/plot_all_events.R -h                                      
usage: R/plot_all_events.R [-h] [-m MASTERREF] [-o RESULTSDIR] [-c CRITERIA]

optional arguments:
  -h, --help            show this help message and exit
  -m MASTERREF, --masterref MASTERREF
                        File path to master reference file
  -o RESULTSDIR, --resultsdir RESULTSDIR
                        Output directory
  -c CRITERIA, --criteria CRITERIA
                        Calling criteria [default]
```

## What `plot_all_events.R` does

### [For each patient](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L154)

* [Read in SNV+SV table](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L157-L159)
  * Only reviewed SNV variants will be plotted --  _High_ for `call_confidence` column
  * SV variants with _in-frame protein fusion_ annotation will be included
* [Melting](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L94-L131) patient-level [table](filter-calls.md#example-of-the-patient-level-table) into a maf-like format
  * Melting VAFs and call status
* [Process maf-like format](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L44-L91)
  * [Filter out ACCESS not covered variants](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L65-L69)
  * [Process SV calls](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L72-L87)
    * Change `Hugo_Symbol` \(which for SVs are `gene1__gene2`\) -- [arrange 2 genes by alphabetical order](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L74-L76) \(i.e. EML4\_\_ALK & ALK\_\_EML4 both to ALK\_\_EML4\)
    * Arrange chromosomes in the alphabetical order of two genes
    * [Aggregate allele counts across SV events](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L12-L41) with the same processed`Hugo_Symbol`
* [Converting sample id to timepoints](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L169-L173) in the master reference
* [Graphing](https://github.com/msk-access/access_data_analysis/blob/17a26eea455707c82824493ebc597d9850d47e82/R/plot_all_events.R#L181-L218)
  * If CNA detected
  * If no CNA is detected

#### 

