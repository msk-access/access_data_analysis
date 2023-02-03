from pathlib import Path

import typer
import pandas as pd
import arrow
from datetime import datetime


def validate_date(date_string):
    date_format = [
        "MM/DD/YY",
        "M/D/YY",
        "MM/D/YY",
        "M/DD/YY",
        "MM/DD/YYYY",
        "YYYY/MM/DD",
        "YYYY-MM-DD",
    ]
    for fmt in date_format:
        try:
            return arrow.get(date_string, fmt).date()
        except ValueError:
            pass
        except:
            print("Something else went wrong")
    raise ValueError("no valid date format found")


def main(
    input: Path = typer.Option(
        ...,
        "--input",
        "-i",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input file with the information to convert dates to days",
    ),
    timepoint_label_for_baseline_first: str = typer.Option(
        "C1D1",
        "--timepoint1",
        "-t1",
        help="timepoint name which in the timepoint column to use for the baseline date, first preference",
    ),
    timepoint_label_for_baseline_second: str = typer.Option(
        "",
        "--timepoint2",
        "-t2",
        help="timepoint name which in the timepoint column to use for the baseline date, second preference",
    ),
    timepoint_label_for_baseline_third: str = typer.Option(
        "",
        "--timepoint3",
        "-t3",
        help="timepoint name which in the timepoint column to use for the baseline date, third preference",
    ),
    output_file: str = typer.Option(
        "output.txt",
        "--output",
        "-o",
        help="Name of the output file",
    ),
):

    """
    Tool to do the following operations:
    A. Reads meta data file, and based on the timepoint information given convert them to days for a samples belonging to a given patient_id
    B. Supports following date formats: 'MM/DD/YY','M/D/YY','MM/D/YY','M/DD/YY','MM/DD/YYYY','YYYY/MM/DD'

    Requirement:
    pandas; typer; arrow

    """

    # Read input file
    i_df = pd.read_csv(input, sep="\t", comment="#", low_memory=False)
    # group by cmo_patient_id
    grouped = i_df.groupby("cmo_patient_id")
    keys = grouped.groups.keys()
    df_list = []
    # tarverse via cmo_patient_id to get associated samples
    for i in keys:
        t_df = pd.DataFrame()
        t_df = grouped.get_group(i)
        baseline_date = None

        # Get the baseline date
        if len(t_df) > 1:
            try:
                baseline_date = t_df.loc[
                    t_df["timepoint"] == timepoint_label_for_baseline_first,
                    "collection_date",
                ].iloc[0]
                baseline_date = validate_date(baseline_date)
            except IndexError:
                print(
                    i,
                    "patient does not have first preference timepoint:",
                    timepoint_label_for_baseline_first,
                )
                print(
                    "We will try to use second timepoint if available to use as baseline\n"
                )
                if timepoint_label_for_baseline_second:
                    try:
                        baseline_date = t_df.loc[
                            t_df["timepoint"] == timepoint_label_for_baseline_second,
                            "collection_date",
                        ].iloc[0]
                        baseline_date = validate_date(baseline_date)
                    except IndexError as e:
                        print(
                            i,
                            "patient does not have second preference timepoint:",
                            timepoint_label_for_baseline_second,
                            "\n",
                        )
                        print(e)
                        if timepoint_label_for_baseline_third:
                            try:
                                baseline_date = t_df.loc[
                                    t_df["timepoint"]
                                    == timepoint_label_for_baseline_third,
                                    "collection_date",
                                ].iloc[0]
                                baseline_date = validate_date(baseline_date)
                            except IndexError as e:
                                print(
                                    i,
                                    "patient does not have third preference timepoint:",
                                    timepoint_label_for_baseline_third,
                                    "\n",
                                )
                                print(e)
                                exit(1)
        else:
            baseline_date = str(t_df["collection_date"])
            baseline_date = validate_date(baseline_date)
        # convert to days
        days_list = []
        for a, b in zip(t_df["collection_date"], t_df["timepoint"]):
            fmt_date = validate_date(a)
            delta = fmt_date - baseline_date
            days_list.append(delta.days)
        # make list of modified dataframes
        t_df_copy = t_df.copy(deep=True)
        t_df_copy["collection_in_days"] = days_list
        df_list.append(t_df_copy)
    # merge and write the dataframe
    results = pd.concat(df_list, axis=0, join="outer")
    results.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    typer.run(main)
