---
description: Short descriptions on the steps of analysis
---

# Overview of Analysis Workflow

The pipeline aims to generate uniform and useful outputs for analyst in preliminary stage of analysis.

Example command to run through the pipeline:

## [1. Compile reads](compile-reads.md)

```bash
> Rscript R/compile_reads.R -m $PATH/TO/master_file.csv -o $PATH/TO/results_folder
```

## [2. Filter calls](filter-calls.md)

```bash
> Rscript R/filter_calls.R -m $PATH/TO/master_file.csv -o $PATH/TO/results_folder
```

## [3. SV incorporation](sv-incorporation.md)

```bash
> Rscript R/SV_incorporation.R -m $PATH/TO/manifest_file.tsv -o $PATH/TO/results_folder
```

## [4. CNA processing](cna-processing.md)

```bash
> Rscript R/CNA_processing.R -m $PATH/TO/manifest_file.tsv -o $PATH/TO/results_folder
```

## [5. Plot all events](overview-of-analysis-workflow.md)

```bash
> Rscript R/plot_all_events.R -m $PATH/TO/manifest_file.tsv -o $PATH/TO/results_folder
```

## [6. Generate HTML report](overview-of-analysis-workflow.md#6-generate-html-report)

<pre class="language-bash" data-overflow="wrap"><code class="lang-bash"><strong>> Rscript ~/github/access_data_analysis/reports/create_report.R -md -t ~/github/access_data_analysis/reports/template_days.Rmd -p C-L6H8E2 -r ../results_20Jan2023/results_stringent_hc/C-L6H8E2_SNV_table.csv -tt "Melanoma" -m ../manifest_noDate_days.tsv -o C-L6H8E2_days.html -rc ../results_20Jan2023/CNA_final_call_set -d P-0022907 -ds P-0022907-T01-IM6 -dm /juno/work/ccs/shared/resources/impact/facets/all/P-00229/P-0022907-T01-IM6_P-0022907-N01-IM6/default/P-0022907-T01-IM6_P-0022907-N01-IM6.ccf.maf
</strong></code></pre>

##
