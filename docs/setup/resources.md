---
description: Description of resource files and executables
---

# Resources

There are various resource files and executables needed for this pipeline. If you are working for JUNO and in the `bergerm1` user group, you should be fine as default options will work fine for you. For other users, here are a list of resources needed in various steps in the pipeline, and their descriptions

## [Compile Reads](../analysis/compile-reads.md)

* Pooled bam directory
  * Directory containing list of donor bams \(unfiltered\) to be genotyped for systematic artifact filtering
* Fasta
  * Hg19 human reference fasta
* Genotyper 
  * Path to the GBCMS genotyper executable 
* DMP IMPACT Github Repository 
  * Repository of DMP IMPACT data updated daily through the cbio enterprise github 
* DMP IMPACT raw data
  * Mirror bam directory
    * Directory containing list of DMP IMPACT bams 
  * Mirror bam key file
    * File containing DMP ID - BAM ID mapping
  * Need to talk to [Aijaz Syed](mailto:syeda1@mskcc.org) about 12-245 access

## [Filter Calls](../analysis/filter-calls.md)

* Mirror bam key file 
* CH list
  * list of signed out CH calls from DMP

## [SV Incorporation](../analysis/sv-incorporation.md)

* Mirror bam key file 

