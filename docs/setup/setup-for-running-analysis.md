---
description: Master reference file descriptions
---

# Setup for Running Analysis

## Master reference file 

An example of [this file](https://github.com/msk-access/access_data_analysis/blob/master/data/example_master_file.csv) can be found in the `data/` folder

{% hint style="danger" %}
For not required columns, leave the cell blank if you don't have the information
{% endhint %}

<table>
  <thead>
    <tr>
      <th style="text-align:left">Column Names</th>
      <th style="text-align:left">Information Specified</th>
      <th style="text-align:left">Specified format (If any)</th>
      <th style="text-align:left">Notes</th>
      <th style="text-align:left">Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align:left">cmo_patient_id</td>
      <td style="text-align:left">Patient ID</td>
      <td style="text-align:left">None</td>
      <td style="text-align:left">Results are presented per unique patient ID</td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">cmo_sample_id_plasma</td>
      <td style="text-align:left">Plasma Sample ID</td>
      <td style="text-align:left">None</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">cmo_sample_id_normal</td>
      <td style="text-align:left">Buffy Coat Sample ID</td>
      <td style="text-align:left">None</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">N</td>
    </tr>
    <tr>
      <td style="text-align:left">bam_path_normal</td>
      <td style="text-align:left">Unfiltered buffy coat bam</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">N</td>
    </tr>
    <tr>
      <td style="text-align:left">paired</td>
      <td style="text-align:left">Whether the plasma has buffy coat</td>
      <td style="text-align:left">Paired/Unpaired</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">sex</td>
      <td style="text-align:left">Sex</td>
      <td style="text-align:left">M/F</td>
      <td style="text-align:left">Unrequired</td>
      <td style="text-align:left">N</td>
    </tr>
    <tr>
      <td style="text-align:left">collection_date</td>
      <td style="text-align:left">Collection time points for graphing</td>
      <td style="text-align:left">
        <p>dates (m/d/y)</p>
        <p>OR</p>
        <p>character strings (i.e. the sample IDs)</p>
      </td>
      <td style="text-align:left">the format should be consistent within the file</td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">dmp_patient_id</td>
      <td style="text-align:left">DMP patient ID</td>
      <td style="text-align:left">*Patient IDs*</td>
      <td style="text-align:left">All DMP samples from this patient ID will be pulled</td>
      <td style="text-align:left">N</td>
    </tr>
    <tr>
      <td style="text-align:left">bam_path_plasma_duplex</td>
      <td style="text-align:left">Duplex bam</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">bam_path_plasma_simplex</td>
      <td style="text-align:left">Simplex bam</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left"></td>
      <td style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">maf_path</td>
      <td style="text-align:left">maf file</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left">fillout_filtered.maf (required columns <a href="setup-for-running-analysis.md#required-columns-for-maf-file">here</a>)</td>
      <td
      style="text-align:left">Y</td>
    </tr>
    <tr>
      <td style="text-align:left">cna_path</td>
      <td style="text-align:left">cna file</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left">sample level cna file (<a href="cna-result-processing.md">helper script included</a>)</td>
      <td
      style="text-align:left">N</td>
    </tr>
    <tr>
      <td style="text-align:left">sv_path</td>
      <td style="text-align:left">sv file</td>
      <td style="text-align:left">Absolute file paths</td>
      <td style="text-align:left">&lt;code&gt;&lt;/code&gt;</td>
      <td style="text-align:left">N</td>
    </tr>
  </tbody>
</table>

{% hint style="warning" %}
 Creating this file might be a hassle. Helper script could possibly be made to help with this
{% endhint %}

## Required Columns for maf file

```text
Hugo_Symbol,Chromosome,Start_Position,End_Position,Tumor_Sample_Barcode,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,D_t_alt_count_fragment
```



