---
description: Description of resource files and executables
---

# Resources

There are various resource files and executables needed for this pipeline. If you are working on JUNO, you should be fine as default options will work fine for you. For other users, here are a list of resources needed in various steps in the pipeline, and their descriptions

## [Compile Reads](../analysis/compile-reads.md)

* Pooled bam directory
  * Directory containing list of donor bams (unfiltered) to be genotyped for systematic artifact filtering
  * Default:`/work/access/production/resources/msk-access/current/novaseq_curated_duplex_bams_dmp/current/`
* Fasta
  * Hg19 human reference fasta
  * Default:`/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta`
* Genotyper&#x20;
  * Path to the GBCMS genotyper executable&#x20;
  * Default: `/ifs/work/bergerm1/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample`
* DMP IMPACT Github Repository&#x20;
  * Repository of DMP IMPACT data updated daily through the cbio enterprise github&#x20;
  * Default: `/juno/work/access/production/resources/cbioportal/current/mskimpact`
* DMP IMPACT raw data
  * Mirror bam directory
    * Directory containing list of DMP IMPACT bams&#x20;
    * Default: `/juno/res/dmpcollab/dmpshare/share/irb12_245/`
  * Mirror bam key file -- ONLY 'IM' (SOLID TISSUE) SAMPLES ARE GENOTYPED
    * File containing DMP ID - BAM ID mapping
    * Default: `/juno/res/dmpcollab/dmprequest/12-245/key.txt`
  * Need to talk to [Aijaz Syed](mailto:syeda1@mskcc.org) about 12-245 access

\>

## [Filter Calls](../analysis/filter-calls.md)

* CH list
  * list of signed out CH calls from DMP
  * Default: `/juno/work/access/production/resources/dmp_signedout_CH/current/signedout_CH.txt`

## [SV Incorporation](../analysis/sv-incorporation.md)

* DMP IMPACT Github Repository&#x20;
  * Repository of DMP IMPACT data updated daily through the cbio enterprise github&#x20;
  * Default: `/juno/work/access/production/resources/cbioportal/current/mskimpact`
