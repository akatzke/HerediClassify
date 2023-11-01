#!/usr/bin/env python3

import pathlib
from collections.abc import Callable

from refactoring.information import classification_information
from refactoring.var_type import VARTYPE_GROUPS
from refactoring.clinvar_missense import check_clinvar_missense
from refactoring.clinvar_splicing import check_clinvar_splicing
from refactoring.variant import VariantInfo, TranscriptInfo
from refactoring.clinvar_utils import CLINVAR_TYPE, ClinVar


def annotate_clinvar(
    variant: VariantInfo, transcripts: list[TranscriptInfo], path_clinvar: pathlib.Path
) -> dict[CLINVAR_TYPE, ClinVar]:
    """
    Manage ClinVar annotation
    For both splicing and missense variants
    """
    if any(var_type in VARTYPE_GROUPS.MISSENSE.value for var_type in variant.var_type):
        clinvar_same_aa, clinvar_diff_aa = check_clinvar_missense(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_aa = ClinVar(False, CLINVAR_TYPE.SAME_AA_CHANGE, None, [])
        clinvar_diff_aa = ClinVar(False, CLINVAR_TYPE.DIFF_AA_CHANGE, None, [])
    if any(var_type in VARTYPE_GROUPS.INTRONIC.value for var_type in variant.var_type):
        clinvar_same_pos, clinvar_splice_site = check_clinvar_splicing(
            variant, transcripts, path_clinvar
        )
    else:
        clinvar_same_pos = ClinVar(False, CLINVAR_TYPE.SAME_NUCLEOTIDE, None, [])
        clinvar_splice_site = ClinVar(False, CLINVAR_TYPE.SAME_SPLICE_SITE, None, [])
    clinvar_dict = {
        CLINVAR_TYPE.SAME_AA_CHANGE: clinvar_same_aa,
        CLINVAR_TYPE.DIFF_AA_CHANGE: clinvar_diff_aa,
        CLINVAR_TYPE.SAME_NUCLEOTIDE: clinvar_same_pos,
        CLINVAR_TYPE.SAME_SPLICE_SITE: clinvar_splice_site,
    }
    return clinvar_dict


def get_annotate_clinvar() -> tuple[Callable, tuple[classification_information, ...]]:
    """
    Get function for clinvar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar,
        (
            classification_information.VARIANT,
            classification_information.TRANSCRIPT,
            classification_information.CLINVAR_PATH,
        ),
    )
