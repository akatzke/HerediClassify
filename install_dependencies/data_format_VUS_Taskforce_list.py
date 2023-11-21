#!/usr/bin/env python3

import argparse
import pathlib
import pandas as pd

parser = argparse = argparse.ArgumentParser()

parser.add_argument(
    "-i", "--input", default="", help="path to input json file", type=str
)
parser.add_argument(
    "-o",
    "--output",
    default="",
    help="output file path. If not given will default to standard output.",
)
args = parser.parse_args()


def format_VUS_Taskforce_list(path: pathlib.Path, out_path: pathlib.Path) -> None:
    """
    Convert .xslx spreadsheet into bed file
    """
    vus_xlsx = pd.read_excel(path)
    vus_xlsx.rename(
        columns={
            "Chr.": "#chr",
            "Position Start (hg38)": "start",
            "Position Ende (hg38)": "end",
            "Protein Anfang": "prot_start",
            "Protein Ende": "prot_end",
            "Quelle": "source",
            "Gen": "gene",
            "Name Proteindomäne": "domain_name",
        },
        inplace=True,
    )
    vus_xlsx["strand"] = ""
    vus_checked_pos = vus_xlsx.apply(lambda row: check_start_end_position(row), axis=1)
    vus_xlsx_bed = vus_checked_pos.loc[
        :,
        [
            "#chr",
            "start",
            "end",
            "gene",
            "domain_name",
            "strand",
            "prot_start",
            "prot_end",
            "source",
        ],
    ]
    vus_rename_chr = vus_xlsx_bed.apply(lambda row: rename_chr(row), axis=1)
    vus_rename_chr.to_csv(out_path, header=True, index=False, sep="\t")


def check_start_end_position(row: pd.Series) -> pd.Series:
    """
    Check that start position is smaller than end position
    If that's not the case exchange the two
    """
    if row.start > row.end:
        start = row.start
        end = row.end
        row.start = end
        row.end = start
        row["strand"] = "-"
    else:
        row["strand"] = "+"
    if row.start == row.end:
        raise ValueError(
            f"Position of {row.domain_name} has the same start and end location: {row.start}"
        )
    return row


def rename_chr(row: pd.Series) -> pd.Series:
    """
    Change chromosome name from "11" to "chr11"
    """
    row["#chr"] = "chr" + str(row["#chr"])
    return row


def main():
    """
    Execute VUS Taksforce list formatting
    """
    if args.input == "":
        raise ValueError("No input file provided.")
    input_path = pathlib.Path(args.input)
    if args.output == "":
        file_name = f"{input_path.stem}.bed"
        out_path = input_path.parent / file_name
    else:
        out_path = pathlib.Path(args.output)
    format_VUS_Taskforce_list(input_path, out_path)


if __name__ == "__main__":
    main()
