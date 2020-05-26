---
description: Master reference file descriptions
---

# Setup for Running Analysis

## Master reference file 

An example of [this file](https://github.com/msk-access/access_data_analysis/blob/master/data/example_master_file.csv) can be found in the `data/` folder

| Column Names | Information Specified | Specified format \(If any\) | Notes |
| :--- | :--- | :--- | :--- |
| cmo\_patient\_id | Patient ID | None | Results are presented per unique patient ID |
| cmo\_sample\_id\_plasma | Plasma Sample ID | None |  |
| cmo\_sample\_id\_normal | Buffy Coat Sample ID | None |  |
| bam\_path\_normal | Unfiltered buffy coat bam | Absolute file paths |  |
| paired | Whether the plasma has buffy coat | Paired/Unpaired |  |
| sex | Sex | M/F |  |
| collection\_date | Collection time points for graphing | actual dates/character strings |  |
| dmp\_patient\_id | DMP patient ID | \*Patient IDs\* | All DMP samples from this patient ID will be pulled |
| bam\_path\_plasma\_duplex | Duplex bam | Absolute file paths |  |
| bam\_path\_plasma\_simplex | Simplex bam | Absolute file paths |  |
| maf\_path | maf file | Absolute file paths |  |
| cna\_path | cna file | Absolute file paths | sample level cna file \(helper script included\) |
| sv\_path | sv file | Absolute file paths |  |

{% hint style="warning" %}
 Creating this file might be a hassle. Helper script could possibly be made to help with this
{% endhint %}



