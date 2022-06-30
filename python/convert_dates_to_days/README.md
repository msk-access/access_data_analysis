# Convert dates to days

Tool to do the following operations:

* Reads meta data file, and based on the timepoint information given, convert them to days for a samples belonging to a given patient\_id
* Supports following date formats:
  * MM/DD/YY
  * M/D/YY
  * MM/D/YY
  * M/DD/YY
  * MM/DD/YYYY
  * YYYY/MM/DD

## Requirements

* pandas
* typing
* arrow

## Example command

```bash
python convert_dates_to_days.py -i ./example_input.txt -t2 "SCREEN"
```

## Usage

```bash
> python convert_dates_to_days.py --help
Usage: convert_dates_to_days.py [OPTIONS]

  Tool to do the following operations: A. Reads meta data file, and based on
  the timepoint information given convert them to days for a samples
  belonging to a given patient_id B. Supports following date formats:
  'MM/DD/YY','M/D/YY','MM/D/YY','M/DD/YY','MM/DD/YYYY','YYYY/MM/DD'

  Requirement: pandas; typer; arrow

Options:
  -i, --input FILE        Input file with the information to convert dates to
                          days  [required]

  -t1, --timepoint1 TEXT  Column name which has timpoint information to use
                          the baseline date, first preference  [default: C1D1]

  -t2, --timepoint2 TEXT  Column name which has timpoint information to use
                          the baseline date, second preference  [default: ]

  -o, --output TEXT       Name of the output file  [default: output.txt]
  --install-completion    Install completion for the current shell.
  --show-completion       Show completion for the current shell, to copy it or
                          customize the installation.

  --help                  Show this message and exit.
```
