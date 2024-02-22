#!/usr/bin/env python3

import pathlib
import logging
from collections.abc import Callable

import pandas as pd

from information import Classification_Info, Info
from variant import VariantInfo, TranscriptInfo
from var_type import VARTYPE_GROUPS
from clinvar_utils import (
    ClinVar,
    ClinVar_Type,
    filter_gene,
    create_ClinVar,
    get_affected_transcript,
)
from custom_exceptions import No_transcript_with_var_type_found
from clinvar_missense import (
    extract_var_codon_info,
    extract_clinvar_entries_missense,
    construct_clinvar_prot_change,
)
from similarity_score import (
    get_similarity_score,
    get_similarity_score_clinvar,
)
from acmg_rules.computation_evidence_utils import Threshold
from format_spliceai import format_spliceai

logger = logging.getLogger("GenOtoScope_Classify.clinvar.missense_similarity_score")


def check_clinvar_missense_similarity(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
    path_similarity_score: pathlib.Path,
    direction_similarity_score: str,
    threshold: Threshold,
) -> ClinVar:
    """
    Check ClinVar for entries supporting pathogenicity of missense variant
    Only looks for different amino acids
    Also ensures that splicing is not predicted for any of the already classified variants
    """
    threshold_spliceAI = max(threshold.thresholds)
    # In case the variant is not an SNV, return empty ClinVar result
    if len(variant.var_obs) != 1 or len(variant.var_ref) != 1:
        logger.warning(
            "Variant is not a SNV. PS1/PM5 currently not implemented for delins."
        )
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return ClinVar_diff_aa
    try:
        affected_transcript, ref_transcript = get_affected_transcript(
            transcripts, VARTYPE_GROUPS.MISSENSE
        )
    except No_transcript_with_var_type_found:
        logger.warning("No transcript with variant type missense found.")
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return ClinVar_diff_aa

    var_codon_info = extract_var_codon_info(
        variant, ref_transcript, affected_transcript
    )
    clinvar_same_codon = extract_clinvar_entries_missense(
        path_clinvar,
        variant.chr,
        var_codon_info["genomic_pos"],
        var_codon_info["intersects_intron_at"],
    )
    clinvar_same_codon_aa: pd.DataFrame = clinvar_same_codon.apply(
        construct_clinvar_prot_change,
        var_codon_info=var_codon_info,
        axis=1,
    )
    if not clinvar_same_codon_aa.empty:
        clinvar_same_codon_aa_splicing = format_spliceai(clinvar_same_codon_aa)
        clinvar_same_codon_aa_filtered = filter_gene(
            clinvar_same_codon_aa_splicing, variant.gene_name
        )
        clinvar_same_codon_aa_similarity: pd.DataFrame = filter_similarity(
            clinvar_same_codon_aa_filtered,
            var_codon_info,
            path_similarity_score,
            direction_similarity_score,
        )
        if clinvar_same_codon_aa_similarity.empty:
            clinvar_diff_aa = pd.DataFrame()
        else:
            clinvar_same_codon_aa_spliceAI = clinvar_same_codon_aa_similarity[
                clinvar_same_codon_aa_similarity.SpliceAI_max < threshold_spliceAI
            ]
            clinvar_diff_aa = clinvar_same_codon_aa_spliceAI[
                clinvar_same_codon_aa_spliceAI.prot_alt != var_codon_info["prot_alt"]
            ]
    else:
        clinvar_diff_aa = pd.DataFrame()

    ClinVar_diff_aa = create_ClinVar(clinvar_diff_aa, ClinVar_Type.DIFF_AA_CHANGE)

    return ClinVar_diff_aa


def filter_similarity(
    clinvar_same_codon_aa: pd.DataFrame,
    var_codon_info: dict,
    path_similarity_score: pathlib.Path,
    direction: str,
) -> pd.DataFrame:
    """
    Filter clinvar dataframe for similarity score
    Direction of score is determined by the name of the tool
    """
    clinvar_similarity = get_similarity_score_clinvar(
        clinvar_same_codon_aa, path_similarity_score
    )
    var_score = get_similarity_score(var_codon_info, path_similarity_score)
    if direction == "less":
        clinvar_filtered = clinvar_similarity[
            clinvar_similarity.similarity_score <= var_score
        ]
    elif direction == "greater":
        clinvar_filtered = clinvar_similarity[
            clinvar_similarity.similarity_score >= var_score
        ]
    else:
        raise ValueError(
            f"No direction for filtering is defiend for the similarity score {path_similarity_score.stem}."
        )
    return clinvar_filtered


def get_check_clinvar_missense_similarity(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for clinvar annotation with similarity score and needed classification_information objects
    """
    return (
        check_clinvar_missense_similarity,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH_SPLICEAI,
            class_info.SIMILARITY_SCORE_PATH,
            class_info.SIMILARITY_SOCRE_DIRECTION,
            class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
        ),
    )
