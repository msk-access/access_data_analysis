---
description: Intermediate files are generated in a internal structure
---

# Intermediate File Organization

There are intermediate files generated with each step in the `/output/directory` , here is a diagram for its organization

```
.s
+-- C-000001
|   +-- C-000001_all_unique_calls.maf
|   +-- C-000001_impact_calls.maf
|   +-- C-000001_sample_sheet.tsv
|   +-- C-000001_genotype_metadata.tsv
        #plasma sample mafs
|   +-- C-000001-L001-d-SIMPLEX_genotyped.maf
|   +-- C-000001-L001-d-DUPLEX_genotyped.maf
|   +-- C-000001-L001-d-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000001-L001-d-ORG-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000001-L002-d-SIMPLEX_genotyped.maf
|   +-- C-000001-L002-d-DUPLEX_genotyped.maf
|   +-- C-000001-L002-d-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000001-L002-d-ORG-SIMPLEX-DUPLEX_genotyped.maf
|   +-- ...
        #buffy coats
|   +-- C-000001-N001-d-STANDARD_genotyped.maf
|   +-- C-000001-N001-d-ORG-STD_genotyped.maf
|   +-- ...
        #DMP samples
|   +-- P-1000000-T01-IM6-STANDARD_genotyped.maf
|   +-- P-1000000-T01-IM6-ORG-STD_genotyped.maf
|   +-- P-1000000-N01-IM6-STANDARD_genotyped.maf
|   +-- P-1000000-N01-IM6-ORG-STD_genotyped.maf
|   +-- ...
+-- C-000002
|   +-- C-000002_all_unique_calls.maf
|   +-- C-000002_impact_calls.maf
|   +-- C-000002_sample_sheet.tsv
|   +-- C-000002_genotype_metadata.tsv
        #plasma sample mafs
|   +-- C-000002-L001-d-SIMPLEX_genotyped.maf
|   +-- C-000002-L001-d-DUPLEX_genotyped.maf
|   +-- C-000002-L001-d-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000002-L001-d-ORG-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000002-L002-d-SIMPLEX_genotyped.maf
|   +-- C-000002-L002-d-DUPLEX_genotyped.maf
|   +-- C-000002-L002-d-SIMPLEX-DUPLEX_genotyped.maf
|   +-- C-000002-L002-d-ORG-SIMPLEX-DUPLEX_genotyped.maf
|   +-- ...
        #buffy coats
|   +-- C-000002-N001-d-STANDARD_genotyped.maf
|   +-- C-000002-N001-d-ORG-STD_genotyped.maf
|   +-- ...
        #DMP samples
|   +-- P-2000000-T01-IM6-STANDARD_genotyped.maf
|   +-- P-2000000-T01-IM6-ORG-STD_genotyped.maf
|   +-- P-2000000-N01-IM6-STANDARD_genotyped.maf
|   +-- P-2000000-N01-IM6-ORG-STD_genotyped.maf
|   +-- ...
+-- ... (other patient directories)        
+-- pooled
|   +-- all_all_unique.maf
|   +-- pooled_metadata.tsv
        #donor samples
|   +-- DONOR1-STANDARD_genotyped.maf
|   +-- DONOR1-ORG-STD_genotyped.maf
|   +-- DONOR2-STANDARD_genotyped.maf
|   +-- DONOR2-ORG-STD_genotyped.maf
|   +-- ...
+-- results_stringent
|   +-- C-000001_SNV_table.csv
|   +-- C-000002_SNV_table.csv
|   +-- ...
+-- results_stringent_combined
|   +-- C-000001_table.csv
|   +-- C-000002_table.csv
|   +-- ...
+-- CNA_final_call_set
|   +-- C-000001_cna_final_call_set.txt
|   +-- C-000002_cna_final_call_set.txt
|   +-- ...
+-- plots
|   +-- C-000001_all_events.pdf
|   +-- C-000002_all_events.pdf
|   +-- ...
```
