#!/usr/bin/env python3

import logging
import pathlib
import pandas as pd
from typing import Generator
from collections.abc import Iterable

from refactoring.variant import TranscriptInfo, VarType
from refactoring.variant_annotate import ClinVar, ClinVar_Type

logger = logging.getLogger("GenOtoScope_Classify.genotoscope_clinvar")

path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")


def get_affected_transcript(
    transcripts: Iterable[TranscriptInfo], var_types: Iterable[VarType]
) -> TranscriptInfo:
    for transcript in transcripts:
        if any(var_type in transcript.var_type for var_type in var_types):
            return transcript
    raise ValueError(f"No transcript has {var_types}")


def convert_vcf_gen_to_df(vcf_generator: Generator) -> pd.DataFrame:
    """
    Covnerts cyvcf generator into a pd.DataFrame
    """
    names = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
    df = pd.DataFrame(columns=names)
    for entry in vcf_generator:
        clinvar_split_str = str(entry).split("\t")
        clinvar_dict = dict(zip(names, clinvar_split_str))
        df = pd.concat([df, pd.DataFrame([clinvar_dict])], axis=0, ignore_index=True)
    df_format_info = format_info(df)
    return df_format_info


def format_info(data: pd.DataFrame) -> pd.DataFrame:
    """
    Format the info column from ClinVar.vcf file as depicted in cyvcf
    """
    info = data["info"]
    info_split = [entry.split("=") for entry in info]
    info_split = [entry.split("\n")[0].split(";") for entry in info]
    processed_info = [[item.split("=") for item in entry] for entry in info_split]
    dict_info = [
        {item[0]: item[1] for item in entry if len(item) > 1}
        for entry in processed_info
    ]
    return pd.concat([data, pd.DataFrame(dict_info)], axis=1)


def create_ClinVar(clinvar: pd.DataFrame, type: ClinVar_Type) -> ClinVar:
    """
    From clinvar entries, get highest classification and IDs of ClinVar entries with that classification
    """
    is_pathogenic = False
    highest_classification = None
    clinvar_ids = None
    if not clinvar.empty:
        if any(clinvar.CLNSIG == "Pathogenic"):
            is_pathogenic = True
            highest_classification = "Pathogenic"
            clinvar_ids = list(clinvar[clinvar.CLNSIG == "Pathogenic"].id)
        elif any(clinvar.CLNSIG == "Likely_pathogenic") or any(
            clinvar.CLNSIG == "Pathogenic/Likely_pathogenic"
        ):
            is_pathogenic = True
            highest_classification = "Likely_pathogenic"
            clinvar_ids = list(
                clinvar[clinvar.CLNSIG.str.contains("Likely_pathogenic")].id
            )
    return ClinVar(is_pathogenic, type, highest_classification, clinvar_ids)


def filter_gene(clinvar: pd.DataFrame, gene: str) -> pd.DataFrame:
    """
    Filter out ClinVar entries that don't contain gene in GENEINFO
    """
    clinvar_filtered = clinvar[clinvar.GENEINFO.str.contains(gene)]
    return clinvar_filtered
