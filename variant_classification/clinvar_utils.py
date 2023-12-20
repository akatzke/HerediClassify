#!/usr/bin/env python3

import logging
import pandas as pd
from typing import Generator, Optional
from dataclasses import dataclass, field
from enum import Enum
from collections.abc import Iterable

import pyensembl

from variant import TranscriptInfo
from var_type import VARTYPE_GROUPS
from ensembl import ensembl

logger = logging.getLogger("GenOtoScope_Classify.genotoscope_clinvar")

# path_clinvar = pathlib.Path("/home/katzkean/clinvar/clinvar_20230730.vcf.gz")


class ClinVar_Type(Enum):
    SAME_AA_CHANGE = "same_aa_change"
    DIFF_AA_CHANGE = "diff_aa_change"
    SAME_NUCLEOTIDE = "same_nucleotide"
    SAME_SPLICE_SITE = "same_splice_site"
    REGION = "region"


class ClinVar_Status(Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"


@dataclass
class ClinVar:
    pathogenic: bool
    type: ClinVar_Type
    highest_classification: Optional[ClinVar_Status]
    ids: list[str] = field(default_factory=list)
    associated_ids: list[str] = field(default_factory=list)


def get_affected_transcript(
    transcripts: Iterable[TranscriptInfo], var_types: VARTYPE_GROUPS
) -> tuple[TranscriptInfo, pyensembl.transcript.Transcript]:
    for transcript in transcripts:
        if any(var_type in transcript.var_type for var_type in var_types.value):
            try:
                ref_transcript = ensembl.transcript_by_id(transcript.transcript_id)
                ref_transcript.coding_sequence
            except ValueError or AttributeError:
                continue
            return transcript, ref_transcript
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
    clinvar_ids = []
    if not clinvar.empty:
        if any(clinvar.CLNSIG == "Pathogenic"):
            is_pathogenic = True
            highest_classification = ClinVar_Status.PATHOGENIC
            clinvar_ids = list(clinvar[clinvar.CLNSIG == "Pathogenic"].id)
        elif any(clinvar.CLNSIG == "Likely_pathogenic") or any(
            clinvar.CLNSIG == "Pathogenic/Likely_pathogenic"
        ):
            is_pathogenic = True
            highest_classification = ClinVar_Status.LIKELY_PATHOGENIC
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


def summarise_ClinVars(clinvars: list[ClinVar], type: ClinVar_Type) -> ClinVar:
    """
    Summarise a list of ClinVars into one ClinVar object
    """
    if len(clinvars) == 0:
        return ClinVar(pathogenic=False, highest_classification=None, ids=[], type=type)
    elif len(clinvars) == 1:
        return clinvars[0]
    else:
        pathogenic = False
        highest_classification = None
        ids = []
        for entry in clinvars:
            if entry.pathogenic:
                pathogenic = True
                ids.append(entry.ids)
                if not highest_classification == "Pathogenic":
                    highest_classification = entry.highest_classification
        return ClinVar(pathogenic, type, highest_classification, ids)
