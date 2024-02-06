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


def format_MANE_transcript_list(path: pathlib.Path, out_path: pathlib.Path) -> None:
    """
    From full length filter list, filter out all non transcript entries
    """
    gtf = load_gtf(path)
    transcripts = gtf.transcript_id.unique()
    transcritps_df = pd.DataFrame(transcripts)
    transcritps_df.to_csv(out_path, sep="t", index=False, header=["transcript"])


def load_gtf(path: pathlib.Path) -> pd.DataFrame:
    """
    Load gtf file with correct column names
    Format IDs to return gene_id, transcript_id and exon_number as separate columns
    """
    names = [
        "chr",
        "origin",
        "feature",
        "start",
        "end",
        "dot_1",
        "strand",
        "dot_2",
        "IDs",
    ]
    df = pd.read_csv(path, names=names, sep="\t", comment="#")
    df_formatted = format_IDs(df)
    return df_formatted


def format_IDs(data: pd.DataFrame) -> pd.DataFrame:
    """
    Read in the id column and separate it into gene_id, transcript_id and exon_number
    """
    ids = data.IDs
    ids_split = [entry.split(";") for entry in ids]
    processed_ids = [
        [item.replace('"', "").strip().split(" ") for item in entry]
        for entry in ids_split
    ]
    dict_ids = [
        {item[0]: item[1] for item in entry if len(item) > 1} for entry in processed_ids
    ]
    return pd.concat([data, pd.DataFrame(dict_ids)], axis=1)


def main():
    if args.input == "":
        raise ValueError("No input file provided.")
    input_path = pathlib.Path(args.input)
    if args.output == "":
        file_name = f"{input_path.stem}_transcript_list.csv"
        out_path = input_path.parent / file_name
    else:
        out_path = pathlib.Path(args.output)
    format_MANE_transcript_list(input_path, out_path)


if __name__ == "__main__":
    main()
