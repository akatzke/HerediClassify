#!/usr/bin/env python3

import pathlib
import logging
import os
import argparse

from cyvcf2 import VCF, Writer

logger = logging.getLogger("Download_data.filter_clinvar")

parser = argparse = argparse.ArgumentParser()

parser.add_argument(
    "-i", "--input", default="", help="path to input json file", type=str
)
parser.add_argument(
    "-f",
    "--filter_indel",
    help="in addition to filtering for snvs, also filter for indels",
    type=bool,
)
args = parser.parse_args()


def prepare_clinvar_file_snv(clinvar_path: pathlib.Path) -> None:
    """
    Filter Clinvar for SNVs and quality
    """
    quality_three_star = ["reviewed_by_expert_panel"]
    quality_two_star = [
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    clinvar_snv_path = filter_clinvar_snv(
        clinvar_path,
    )
    clinvar_quality_three_star_path = filter_clinvar_quality(
        clinvar_snv_path, quality_three_star, "three_star"
    )
    clinvar_quality_two_star_path = filter_clinvar_quality(
        clinvar_snv_path, quality_two_star, "two_star"
    )


def prepare_clinvar_file_indel(clinvar_path: pathlib.Path) -> None:
    """
    Filter ClinVar for indel and qualtiy
    """
    quality_three_star = ["reviewed_by_expert_panel"]
    quality_two_star = [
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    clinvar_indel_path = filter_clinvar_small_indel(
        clinvar_path,
    )
    clinvar_quality_three_star_path = filter_clinvar_quality(
        clinvar_indel_path, quality_three_star, "three_star"
    )
    clinvar_quality_two_star_path = filter_clinvar_quality(
        clinvar_indel_path, quality_two_star, "two_star"
    )


def filter_clinvar_quality(
    clinvar_path: pathlib.Path, quality_filter: list[str], quality_filter_name: str
) -> pathlib.Path:
    """
    Filter clinvar for quality metric
    """
    clinvar = VCF(clinvar_path)
    out_name = clinvar_path.stem.split(".")[0] + "_" + f"{quality_filter_name}.vcf"
    out_path = clinvar_path.parent / out_name
    w = Writer(out_path, clinvar)
    for entry in clinvar:
        if entry.INFO["CLNREVSTAT"] in quality_filter:
            w.write_record(entry)
    w.close()
    clinvar.close()
    path_comp = compress_vcf(out_path)
    index_vcf(path_comp)
    return path_comp


def filter_clinvar_snv(clinvar_path: pathlib.Path) -> pathlib.Path:
    """
    Filter out all non SNVs from dataset
    """
    clinvar = VCF(clinvar_path)
    out_name = clinvar_path.stem.split(".")[0] + f"_snv.vcf"
    out_path = clinvar_path.parent / out_name
    w = Writer(out_path, clinvar)
    for entry in clinvar:
        if not len(entry.REF) > 1 and len(entry.ALT) == 1:
            if len(entry.ALT[0]) == 1:
                w.write_record(entry)
    w.close()
    clinvar.close()
    path_comp = compress_vcf(out_path)
    index_vcf(path_comp)
    return path_comp


def filter_clinvar_small_indel(clinvar_path: pathlib.Path) -> pathlib.Path:
    """
    Filter out all non SNVs from dataset
    """
    clinvar = VCF(clinvar_path)
    out_name = clinvar_path.stem.split(".")[0] + f"_small_indel.vcf"
    out_path = clinvar_path.parent / out_name
    w = Writer(out_path, clinvar)
    for entry in clinvar:
        if len(entry.ALT) == 1:
            if not (len(entry.REF) == 1 and len(entry.ALT[0]) == 1) and (
                len(entry.REF) <= 15 and len(entry.ALT[0]) <= 15
            ):
                w.write_record(entry)
    w.close()
    clinvar.close()
    path_comp = compress_vcf(out_path)
    index_vcf(path_comp)
    return path_comp


def compress_vcf(path: pathlib.Path) -> pathlib.Path:
    """
    Compress vcf file and output path to compressed file
    Return path of compressed file
    """
    cmd_compress = f"bgzip {path}"
    os.system(cmd_compress)
    path_compressed = path.parent / str(path.stem + ".vcf.gz")
    return path_compressed


def index_vcf(path: pathlib.Path) -> None:
    """
    Index compressed vcf file
    """
    cmd_index = f"tabix {path}"
    os.system(cmd_index)


def main():
    if args.input == "":
        raise ValueError("No input file provided")
    input_path = pathlib.Path(args.input)
    prepare_clinvar_file_snv(input_path)
    if args.filter_indel:
        prepare_clinvar_file_indel(input_path)


if __name__ == "__main__":
    main()
