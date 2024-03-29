# Table of Contents

- [Table of Contents](#table-of-contents)
- [get\_cbioportal\_variants](#get_cbioportal_variants)
      - [subset\_cpt](#subset_cpt)
      - [subset\_cst](#subset_cst)
      - [subset\_cna](#subset_cna)
      - [subset\_sv](#subset_sv)
      - [subset\_maf](#subset_maf)
      - [read\_tsv](#read_tsv)
      - [read\_ids](#read_ids)
      - [filter\_by\_columns](#filter_by_columns)
      - [filter\_by\_rows](#filter_by_rows)
      - [read\_bed](#read_bed)
      - [check\_if\_covered](#check_if_covered)
      - [get\_row](#get_row)

<a id="get_cbioportal_variants"></a>

# get\_cbioportal\_variants

Requirement:

- pandas
- typing
- typer
- bed_lookup(<https://github.com/msk-access/python_bed_lookup>)

```bash
Usage: get_cbioportal_variants.py [OPTIONS] COMMAND [ARGS]...

Options:
  --install-completion  Install completion for the current shell.
  --show-completion     Show completion for the current shell, to copy it or
                        customize the installation.

  --help                Show this message and exit.

Commands:
  subset-cna  Subset data_CNA.txt file for given set of sample ids.
  subset-cpt  Subset data_clinical_patient.txt file for given set of
              patient...

  subset-cst  Subset data_clinical_samples.txt file for given set of sample...
  subset-maf  Subset MAF/TSV file and mark if an alteration is covered by...
  subset-sv   Subset data_sv.txt file for given set of sample ids.
```

<a id="get_cbioportal_variants.subset_cpt"></a>

#### subset\_cpt

```bash
Usage: get_cbioportal_variants.py subset-cpt [OPTIONS]

  Subset data_clinical_patient.txt file for given set of patient ids.

  Tool to do the following operations: A. Get subset of clinical information
  for samples based on PATIENT_ID in data_clinical_patient.txt file

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -p, --cpt FILE    Clinical Patient file generated by cBioportal repo
                    [default: /work/access/production/resources/cbioportal/cur
                    rent/msk_solid_heme/data_clinical_patient.txt]

  -i, --ids PATH    List of ids to search for in the 'PATIENT_ID' column.
                    Header of this file is 'sample_id'  [default: ]

  --sid TEXT        Identifiers to search for in the 'PATIENT_ID' column. Can
                    be given multiple times  [default: ]

  -n, --name TEXT   Name of the output file  [default:
                    output_clinical_patient.txt]

  -c, --cname TEXT  Name of the column header to be used for sub-setting
                    [default: PATIENT_ID]

  --help            Show this message and exit.
```

<a id="get_cbioportal_variants.subset_cst"></a>

#### subset\_cst

```bash
Usage: get_cbioportal_variants.py subset-cst [OPTIONS]

  Subset data_clinical_samples.txt file for given set of sample ids.

  Tool to do the following operations: A. Get subset of clinical information
  for samples based on SAMPLE_ID in data_clinical_sample.txt file

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -s, --cst FILE    Clinical Sample file generated by cBioportal repo
                    [default: /work/access/production/resources/cbioportal/cur
                    rent/msk_solid_heme/data_clinical_sample.txt]

  -i, --ids PATH    List of ids to search for in the 'SAMPLE_ID' column.
                    Header of this file is 'sample_id'  [default: ]

  --sid TEXT        Identifiers to search for in the 'SAMPLE_ID' column. Can
                    be given multiple times  [default: ]

  -n, --name TEXT   Name of the output file  [default:
                    output_clinical_samples.txt]

  -c, --cname TEXT  Name of the column header to be used for sub-setting
                    [default: SAMPLE_ID]

  --help            Show this message and exit.
```

<a id="get_cbioportal_variants.subset_cna"></a>

#### subset\_cna

```bash
Usage: get_cbioportal_variants.py subset-cna [OPTIONS]

  Subset data_CNA.txt file for given set of sample ids.

  Tool to do the following operations: A. Get subset of samples based on
  column header in data_CNA.txt file

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -c, --cna FILE   Copy Number Variant file generated by cBioportal repo
                   [default: /work/access/production/resources/cbioportal/curr
                   ent/msk_solid_heme/data_CNA.txt]

  -i, --ids PATH   List of ids to search for in the 'header' of the file.
                   Header of this file is 'sample_id'  [default: ]

  --sid TEXT       Identifiers to search for in the 'header' of the file. Can
                   be given multiple times  [default: ]

  -n, --name TEXT  Name of the output file  [default: output_CNA.txt]
  --help           Show this message and exit.
```

<a id="get_cbioportal_variants.subset_sv"></a>

#### subset\_sv

```bash
Usage: get_cbioportal_variants.py subset-sv [OPTIONS]

  Subset data_sv.txt file for given set of sample ids.

  Tool to do the following operations: A. Get subset of structural variants
  based on Sample_ID in data_sv.txt file

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -s, --sv FILE     Structural Variant file generated by cBioportal repo
                    [default: /work/access/production/resources/cbioportal/cur
                    rent/msk_solid_heme/data_sv.txt]

  -i, --ids PATH    List of ids to search for in the 'Sample_ID' column.
                    Header of this file is 'sample_id'  [default: ]

  --sid TEXT        Identifiers to search for in the 'Sample_ID' column. Can
                    be given multiple times  [default: ]

  -n, --name TEXT   Name of the output file  [default: output_sv.txt]
  -c, --cname TEXT  Name of the column header to be used for sub-setting
                    [default: Sample_ID]

  --help            Show this message and exit.
```

<a id="get_cbioportal_variants.subset_maf"></a>

#### subset\_maf

```bash
Usage: get_cbioportal_variants.py subset-maf [OPTIONS]

  Subset MAF/TSV file and mark if an alteration is covered by BED file or
  not

  Tool to do the following operations: A. Get subset of variants based on
  Tumor_Sample_Barcode in data_mutations_extended.txt file B. Mark the
  variants as overlapping with BED file as covered [yes/no], by appending
  "covered" column to the subset MAF

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -m, --maf FILE    MAF file generated by cBioportal repo  [default: /work/acc
                    ess/production/resources/cbioportal/current/msk_solid_heme
                    /data_mutations_extended.txt]

  -i, --ids PATH    List of ids to search for in the 'Tumor_Sample_Barcode'
                    column. Header of this file is 'sample_id'  [default: ]

  --sid TEXT        Identifiers to search for in the 'Tumor_Sample_Barcode'
                    column. Can be given multiple times  [default: ]

  -b, --bed FILE    BED file to find overlapping variants  [default:
                    /work/access/production/resources/msk-
                    access/current/regions_of_interest/current/MSK-
                    ACCESS-v1_0-probe-A.sorted.bed]

  -n, --name TEXT   Name of the output file  [default: output.maf]
  -c, --cname TEXT  Name of the column header to be used for sub-setting
                    [default: Tumor_Sample_Barcode]

  --help            Show this message and exit.
```


<a id="get_cbioportal_variants.read_tsv"></a>

#### read\_tsv

```python
def read_tsv(tsv)
```

Read a tsv file

**Arguments**:

- `maf` _File_ - Input MAF/tsv like format file
  
**Returns**:

- `data_frame` - Output a data frame containing the MAF/tsv

<a id="get_cbioportal_variants.read_ids"></a>

#### read\_ids

```python
def read_ids(sid, ids)
```

make a list of ids

**Arguments**:

- `sid` _tuple_ - Multiple ids as tuple
- `ids` _File_ - File containing multiple ids
  
**Returns**:

- `list` - List containing all ids

<a id="get_cbioportal_variants.filter_by_columns"></a>

#### filter\_by\_columns

```python
def filter_by_columns(sid, tsv_df)
```

Filter data by columns

**Arguments**:

- `sid` _list_ - list of columns to subset over
- `tsv_df` _data_frame_ - data_frame to subset from
  
**Returns**:

- `data_frame` - A copy of the subset of the data_frame

<a id="get_cbioportal_variants.filter_by_rows"></a>

#### filter\_by\_rows

```python
def filter_by_rows(sid, tsv_df, col_name)
```

Filter the data by rows

**Arguments**:

- `sid` _list_ - list of row names to subset over
- `tsv_df` _data_frame_ - data_frame to subset from
- `col_name` _string_ - name of the column to filter using names in the sid
  
**Returns**:

- `data_frame` - A copy of the subset of the data_frame

<a id="get_cbioportal_variants.read_bed"></a>

#### read\_bed

```python
def read_bed(bed)
```

Read BED file using bed_lookup

**Arguments**:

- `bed` _file_ - File ins BED format to read
  
**Returns**:

  object : bed file object to use for filtering

<a id="get_cbioportal_variants.check_if_covered"></a>

#### check\_if\_covered

```python
def check_if_covered(bedObj, mafObj)
```

Function to check if a variant is covered in a given bed file

**Arguments**:

- `bedObj` _object_ - BED file object to check coverage
- `mafObj` _data_frame_ - data frame to check coverage against coordinates using column 'Chromosome' and position column is 'Start_Position'
  
**Returns**:

- `data_frame` - _description_

<a id="get_cbioportal_variants.get_row"></a>

#### get\_row

```python
def get_row(tsv_file)
```

Function to skip rows

**Arguments**:

- `tsv_file` _file_ - file to be read
  
**Returns**:

- `list` - lines to be skipped
