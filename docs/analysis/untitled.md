---
description: Short descriptions on the steps of analysis
---

# Overview of Analysis Workflow

The pipeline aims to generate uniform and useful outputs for analyst in preliminary stage of analysis. 

Example command to run through the pipeline:

## [1. Compile reads](compile-reads.md)

```text
Rscript R/compile_reads.R -m /juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv -o /juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020
```

## [2. Filter calls](filter-calls.md)

```text
Rscript R/filter_calls.R -m /juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv -o /juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020
```

## [3. SV incorporation](sv-incorporation.md)

## [4. CNA processing](cna-processing.md)

## [5. Plot all events](plot-all-events.md)

