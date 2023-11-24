#!/usr/bin/env python3

import pathlib
import logging
from collections.abc import Callable

from information import Classification_Info, Info
from var_type import VARTYPE_GROUPS
from clinvar_missense import check_clinvar_missense
from clinvar_splicing import check_clinvar_splicing
from variant import VariantInfo, TranscriptInfo
from clinvar_utils import ClinVar_Type, ClinVar

logger = logging.getLogger("annotate_clinvar")


def annotate_clinvar(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
) -> dict[ClinVar_Type, ClinVar]:
    """
    Manage ClinVar annotation
    For both splicing and missense variants
    """
    if any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant.var_type):
        clinvar_same_aa, clinvar_diff_aa = check_clinvar_missense(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_aa = ClinVar(False, ClinVar_Type.SAME_AA_CHANGE, None, [])
        clinvar_diff_aa = ClinVar(False, ClinVar_Type.DIFF_AA_CHANGE, None, [])
    if any(var_type in VARTYPE_GROUPS.INTRONIC.value for var_type in variant.var_type):
        clinvar_same_pos, clinvar_splice_site = check_clinvar_splicing(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_pos = ClinVar(False, ClinVar_Type.SAME_NUCLEOTIDE, None, [])
        clinvar_splice_site = ClinVar(False, ClinVar_Type.SAME_SPLICE_SITE, None, [])
    clinvar_dict = {
        ClinVar_Type.SAME_AA_CHANGE: clinvar_same_aa,
        ClinVar_Type.DIFF_AA_CHANGE: clinvar_diff_aa,
        ClinVar_Type.SAME_NUCLEOTIDE: clinvar_same_pos,
        ClinVar_Type.SAME_SPLICE_SITE: clinvar_splice_site,
    }
    return clinvar_dict


def get_annotate_clinvar(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for clinvar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH,
        ),
    )
