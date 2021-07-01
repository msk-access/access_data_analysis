---
description: >-
  Step 5 -- Create a report showing genomic alteration data for all samples of a
  patient.
---

# Create Patient Report

The final step takes the processed data from the previous steps and plots the genomic alterations over all samples of each patient. The report includes several sections with interactive plots:

## 1. Patient information

The first section displays the patient ID, DMP id \(if provided\), tumor type \(if provided\), and each sample. Any provided sample meta-information is also display for each sample.

## 2. Plot of SNV variant allele frequencies

The second section shows SNV/INDEL events are plotted out by VAFs over timepoints. Above the panel it also display sample timepoint annotation, such as treatment information \(if provided\). If you provide IMPACT sample information, it will segregate each mutation by whether it is known to be clonal in IMPACT, subclonal in IMPACT, or is present in ACCESS only. There are additional tabs that display a table of mutation data and methods description.

## 3. Plot of copy number alterations

The third section shows CNAs that are plotted by fold-change\(fc\) for each ACCESS sample and gene. If there are no CNAs, then this section is not displayed.

## 4. Plot of clonal SNV/INDEL VAFs adjusted for copy number

If you provided an IMPACT sample, this last section will show SNV/INDEL events that are plotted out by VAFs over timepoints. However, the VAFs are corrected for IMPACT copy number information. Details of the method are shown under the `Description` tab in this section. Similar to section 2, sample timepoint annotations are shown above the plot.

## Usage

```text
Rscript reports/create_report.R -h                                      
usage: create_report.R [-h] -t TEMPLATE -p PATIENT_ID -r RESULTS -rc
                       CNA_RESULTS_DIR -tt TUMOR_TYPE -m METADATA [-d DMP_ID]
                       [-ds DMP_SAMPLE_ID] [-dm DMP_MAF] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -t TEMPLATE, --template TEMPLATE
                        Path to Rmarkdown template file.
  -p PATIENT_ID, --patient-id PATIENT_ID
                        Patient ID
  -r RESULTS, --results RESULTS
                        Path to CSV file containing mutation and genotype
                        results for the patient.
  -rc CNA_RESULTS_DIR, --cna-results-dir CNA_RESULTS_DIR
                        Path to directory containing CNA results for the
                        patient.
  -tt TUMOR_TYPE, --tumor-type TUMOR_TYPE
                        Tumor type
  -m METADATA, --metadata METADATA
                        Path to file containing meta data. Should contain a
                        'cmo_sample_id_plasma', 'sex', and 'collection_date'
                        columns. Can also optionally include a 'timepoint'
                        column (e.g. for treatment information).
  -d DMP_ID, --dmp-id DMP_ID
                        DMP patient ID (optional).
  -ds DMP_SAMPLE_ID, --dmp-sample-id DMP_SAMPLE_ID
                        DMP sample ID (optional).
  -dm DMP_MAF, --dmp-maf DMP_MAF
                        Path to DMP MAF file (optional).
  -o OUTPUT, --output OUTPUT
                        Output file
```

