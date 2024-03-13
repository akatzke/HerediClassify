#!/usr/bin/env python3

import pathlib
import pandas as pd
import argparse

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


def format_uniprot_rep(path: pathlib.Path, out_path: pathlib.Path) -> None:
    """
    Convert the loaded json into a dataframe with a bed-like strucutre
    Reduce number of columns to 16 in order for the file to be readable by bedtools
    """
    json = pd.read_json(path)
    data = pd.DataFrame(json.unipRepeat.tolist())
    red_data = data.iloc[:, 0:16]
    red_data = red_data.rename(
        columns={"chrom": "#chr", "chromStart": "start", "chromEnd": "end"}
    )
    red_data.to_csv(out_path, sep="\t", index=False)


def main():
    if args.input == "":
        raise ValueError("No input file provided.")
    input_path = pathlib.Path(args.input)
    if args.output == "":
        file_name = f"{input_path.stem}.bed"
        out_path = input_path.parent / file_name
    else:
        out_path = pathlib.Path(args.output)
    format_uniprot_rep(input_path, out_path)


if __name__ == "__main__":
    main()
