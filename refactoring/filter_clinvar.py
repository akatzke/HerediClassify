#!/usr/bin/env python3

import datetime
import pathlib
import logging
import os
from cyvcf2 import VCF, Writer
from datetime import date


logger = logging.getLogger("GenOtoScope_Annotate.filter_clinvar")


def prepare_clinvar_file() -> None:
    clinvar_path = download_clinvar()
    quality = ["reviewed_by_expert_panel"]
    clinvar_quality_path = filter_clinvar_quality(clinvar_path, quality)
    clinvar_quality_snv_path = filter_clinvar_snv(clinvar_quality_path)


def download_clinvar() -> pathlib.Path:
    cwd = pathlib.Path.cwd()
    os.system("mkdir clinvar")
    dir_name = create_clinvar_dir_name()
    clinvar_path = cwd / "clinvar" / dir_name / "clinvar.vcf.gz"
    command_download_clinvar = f"wget -O {clinvar_path} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    os.system(command_download_clinvar)
    return clinvar_path


def create_clinvar_dir_name() -> str:
    """
    Create clinvar dir name
    Schema: clinvar_month_year
    """
    today = date.today()
    dir_name = f"clinvar_{today.year}_{today.month}_{today.day}"
    return dir_name


def filter_clinvar_quality(
    clinvar_path: pathlib.Path, quality_filter: list[str]
) -> pathlib.Path:
    """
    Filter clinvar for quality metric
    """
    clinvar = VCF(clinvar_path)
    out_name = clinvar_path.stem.split(".")[0] + f"_quality_filter.vcf"
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
        if (len(entry.REF) > 1 and len(entry.REF) <= 15) and len(entry.ALT) == 1:
            if len(entry.ALT[0]) > 1 and len(entry.ALT[0]) <= 15:
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
    cmd_index = f"bcftools index -t {path}"
    os.system(cmd_index)
