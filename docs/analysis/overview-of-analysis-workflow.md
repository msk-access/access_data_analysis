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

## [5. Plot all events]()

```bash
> Rscript R/plot_all_events.R -m $PATH/TO/manifest_file.tsv -o $PATH/TO/results_folder
```

## [6. Generate HTML report](overview-of-analysis-workflow.md#6-generate-html-report)

## 

