## Get cbioportal variants

Tool to do the following operations:
* Get subset of variants based on Tumor_Sample_Barcode in MAF file 
* Mark the variants as overlapping with BED file as covered [yes/no], by appending "covered" column to the subset MAF
    
### Requirements
* pandas
* typing
* typer
* bed_lookup(https://github.com/msk-access/python_bed_lookup)

### Example command

```bash
python get_cbioportal_variants.py  --id "Test1" --id "Test2" --id "Test3" 
```

```bash
python get_cbioportal_variants.py  --ids /path/to/ids.txt
```

### Usage

```bash
> python get_cbioportal_variants.py --help
Usage: get_cbioportal_variants.py [OPTIONS]

  Tool to do the following operations: A. Get subset of variants based on
  Tumor_Sample_Barcode in MAF file  B. Mark the variants as overlapping with
  BED file as covered [yes/no], by appending "covered" column to the subset
  MAF

  Requirement: pandas; typing; typer; bed_lookup(https://github.com/msk-
  access/python_bed_lookup)

Options:
  -m, --maf FILE        MAF file generated by cbioportal repo  [default: /work
                        /access/production/resources/cbioportal/current/msk_so
                        lid_heme/data_mutations_extended.txt]

  -i, --ids PATH        List of ids to search for in the
                        'Tumor_Sample_Barcode' column. Header of this file is
                        'sample_id'  [default: ]

  --id TEXT             Identifiers to search for in the
                        'Tumor_Sample_Barcode' column. Can be given multiple
                        times  [default: ]

  -b, --bed FILE        BED file to find overlapping variants  [default:
                        /work/access/production/resources/msk-
                        access/current/regions_of_interest/current/MSK-
                        ACCESS-v1_0-probe-A.sorted.bed]

  -n, --name TEXT       Name of the output file  [default: output.maf]
  --install-completion  Install completion for the current shell.
  --show-completion     Show completion for the current shell, to copy it or
                        customize the installation.

  --help                Show this message and exit.
```