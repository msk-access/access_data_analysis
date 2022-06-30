---
description: Master reference file descriptions
---

# Setup for Running Analysis

## Master reference file

An example of [this file](https://github.com/msk-access/access\_data\_analysis/blob/master/data/example\_master\_file.csv) can be found in the `data/` folder

{% hint style="danger" %}
For not required columns, leave the cell blank if you don't have the information
{% endhint %}

| Column Names               | Information Specified               | Specified format (If any)                                                   | Notes                                                                                                        | Required |
| -------------------------- | ----------------------------------- | --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | -------- |
| cmo\_patient\_id           | Patient ID                          | None                                                                        | Results are presented per unique patient ID                                                                  | Y        |
| cmo\_sample\_id\_plasma    | Plasma Sample ID                    | None                                                                        |                                                                                                              | Y        |
| cmo\_sample\_id\_normal    | Buffy Coat Sample ID                | None                                                                        |                                                                                                              | N        |
| bam\_path\_normal          | Unfiltered buffy coat bam           | Absolute file paths                                                         |                                                                                                              | N        |
| paired                     | Whether the plasma has buffy coat   | Paired/Unpaired                                                             |                                                                                                              | Y        |
| sex                        | Sex                                 | M/F                                                                         | Unrequired                                                                                                   | N        |
| collection\_date           | Collection time points for graphing | <p>dates (m/d/y)</p><p>OR</p><p>character strings (i.e. the sample IDs)</p> | the format should be consistent within the file                                                              | Y        |
| dmp\_patient\_id           | DMP patient ID                      | \*Patient IDs\*                                                             | All DMP samples from this patient ID will be pulled                                                          | N        |
| bam\_path\_plasma\_duplex  | Duplex bam                          | Absolute file paths                                                         |                                                                                                              | Y        |
| bam\_path\_plasma\_simplex | Simplex bam                         | Absolute file paths                                                         |                                                                                                              | Y        |
| maf\_path                  | maf file                            | Absolute file paths                                                         | fillout\_filtered.maf (required columns [here](setup-for-running-analysis.md#required-columns-for-maf-file)) | Y        |
| cna\_path                  | cna file                            | Absolute file paths                                                         | sample level cna file ([helper script included](cna-result-processing.md))                                   | N        |
| sv\_path                   | sv file                             | Absolute file paths                                                         | \<code>\</code>                                                                                              | N        |

{% hint style="warning" %}
Creating this file might be a hassle. Helper script could possibly be made to help with this
{% endhint %}

## Required Columns for maf file

```
Hugo_Symbol,Chromosome,Start_Position,End_Position,Tumor_Sample_Barcode,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,D_t_alt_count_fragment
```
