#!/usr/bin/env python3

import pathlib
import logging
from collections.abc import Callable

import pandas as pd

from cyvcf2 import VCF

from information import Classification_Info, Info
from var_type import VARTYPE_GROUPS
from variant import VariantInfo, TranscriptInfo
from clinvar_utils import (
    ClinVar,
    ClinVar_Type,
    filter_gene,
    create_ClinVar,
    get_affected_transcript,
    convert_vcf_gen_to_df,
)
from custom_exceptions import No_transcript_with_var_type_found
from clinvar_missense import (
    extract_var_codon_info,
    extract_clinvar_entries_missense,
    construct_clinvar_prot_change,
)
from clinvar_splicing import find_corresponding_splice_site
from acmg_rules.computation_evidence_utils import Threshold
from format_spliceai import format_spliceai

logger = logging.getLogger("GenOtoScope_Classify.clinvar_annot_spliceai")


def annotate_clinvar_spliceai_protein(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
    threshold: Threshold,
) -> dict[ClinVar_Type, ClinVar]:
    """
    Manage ClinVar annotation
    For both splicing and missense variants
    """
    threshold_spliceAI = max(threshold.thresholds)
    if len(variant.var_obs) != 1 or len(variant.var_ref) != 1:
        logger.warning(
            "Variant is not a SNV. PS1/PM5 currently not implemented for delins."
        )
        ClinVar_same_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.SAME_AA_CHANGE)
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return {
            ClinVar_Type.SAME_AA_CHANGE: ClinVar_same_aa,
            ClinVar_Type.DIFF_AA_CHANGE: ClinVar_diff_aa,
        }
    try:
        affected_transcript, ref_transcript = get_affected_transcript(
            transcripts, VARTYPE_GROUPS.MISSENSE
        )
    except No_transcript_with_var_type_found:
        logger.warning("No transcript with variant type missense found.")
        ClinVar_same_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.SAME_AA_CHANGE)
        ClinVar_diff_aa = create_ClinVar(pd.DataFrame(), ClinVar_Type.DIFF_AA_CHANGE)
        return {
            ClinVar_Type.SAME_AA_CHANGE: ClinVar_same_aa,
            ClinVar_Type.DIFF_AA_CHANGE: ClinVar_diff_aa,
        }

    var_codon_info = extract_var_codon_info(
        variant, ref_transcript, affected_transcript
    )
    clinvar_same_codon = extract_clinvar_entries_missense(
        path_clinvar,
        variant.chr,
        var_codon_info["genomic_pos"],
        var_codon_info["intersects_intron_at"],
    )
    if "SpliceAI" not in clinvar_same_codon.columns and not clinvar_same_codon.empty:
        logger.warning(f"No SpliceAI column availabel for ClinVar entries.")
        clinvar_diff_aa = pd.DataFrame()
        clinvar_same_aa = pd.DataFrame()
    elif not clinvar_same_codon.empty:
        clinvar_same_codon_splicing = format_spliceai(clinvar_same_codon)
        clinvar_same_codon_aa: pd.DataFrame = clinvar_same_codon_splicing.apply(
            construct_clinvar_prot_change,
            var_codon_info=var_codon_info,
            axis=1,
        )
        # Filter out all ClinVar entries with termination codon
        if not clinvar_same_codon_aa.empty:
            clinvar_same_codon_no_ter = clinvar_same_codon_aa[
                clinvar_same_codon_aa.prot_alt != "Ter"
            ]
        else:
            clinvar_same_codon_no_ter = clinvar_same_codon_aa
        clinvar_same_codon_aa_filtered = filter_gene(
            clinvar_same_codon_no_ter, variant.gene_name
        )
        clinvar_same_codon_aa_spliceAI = clinvar_same_codon_aa_filtered[
            clinvar_same_codon_aa_filtered.SpliceAI_max <= threshold_spliceAI
        ]
        clinvar_same_aa = clinvar_same_codon_aa_spliceAI[
            clinvar_same_codon_aa_spliceAI.prot_alt == var_codon_info["prot_alt"]
        ]
        clinvar_diff_aa = clinvar_same_codon_aa_spliceAI[
            clinvar_same_codon_aa_spliceAI.prot_alt != var_codon_info["prot_alt"]
        ]
    else:
        clinvar_diff_aa = pd.DataFrame()
        clinvar_same_aa = pd.DataFrame()

    ClinVar_diff_aa = create_ClinVar(clinvar_diff_aa, ClinVar_Type.DIFF_AA_CHANGE)
    ClinVar_same_aa = create_ClinVar(clinvar_same_aa, ClinVar_Type.SAME_AA_CHANGE)
    return {
        ClinVar_Type.SAME_AA_CHANGE: ClinVar_same_aa,
        ClinVar_Type.DIFF_AA_CHANGE: ClinVar_diff_aa,
    }


def get_annotate_clinvar_spliceai_protein(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for clinvar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar_spliceai_protein,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH_SPLICEAI,
            class_info.THRESHOLD_SPLICING_PREDICTION_BENIGN,
        ),
    )


def annotate_clinvar_spliceai_splicing(
    variant: VariantInfo,
    transcripts: list[TranscriptInfo],
    path_clinvar: pathlib.Path,
    prediction_dict: dict[str, float],
) -> dict[ClinVar_Type, ClinVar]:
    """
    Manage ClinVar annotation
    For both splicing and missense variants
    """
    if len(variant.var_obs) != 1 or len(variant.var_ref) != 1:
        logger.warning(
            "Variant is not a SNV. PS1/PM5 currently not implemented for delins."
        )
        ClinVar_same_nucleotide = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE
        )
        ClinVar_same_splice_site = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
        )
        return {
            ClinVar_Type.SAME_NUCLEOTIDE: ClinVar_same_nucleotide,
            ClinVar_Type.SAME_SPLICE_SITE: ClinVar_same_splice_site,
        }
    try:
        affected_transcript, ref_transcript = get_affected_transcript(
            transcripts, VARTYPE_GROUPS.INTRONIC
        )
    except No_transcript_with_var_type_found:
        logger.warning(
            "No transcript with variant type splice donor or splice acceptor."
        )
        ClinVar_same_nucleotide = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE
        )
        ClinVar_same_splice_site = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
        )
        return {
            ClinVar_Type.SAME_NUCLEOTIDE: ClinVar_same_nucleotide,
            ClinVar_Type.SAME_SPLICE_SITE: ClinVar_same_splice_site,
        }
    try:
        var_spliceai = prediction_dict["SpliceAI"]
    except KeyError:
        logger.logging(
            "No SpliceAI prediction is given for the variant. ClinVar for splicing variants with SpliceAI assessment can not be assessed."
        )
        ClinVar_same_nucleotide = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE
        )
        ClinVar_same_splice_site = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
        )
        return {
            ClinVar_Type.SAME_NUCLEOTIDE: ClinVar_same_nucleotide,
            ClinVar_Type.SAME_SPLICE_SITE: ClinVar_same_splice_site,
        }

    clinvar = VCF(path_clinvar)
    ### Check ClinVar for pathogenic variants with same nucleotide change
    ### The predicted of the variant under assessment must be similar or higher than the prediction of the known variant
    clinvar_same_pos = clinvar(
        f"{variant.chr}:{variant.genomic_start}-{variant.genomic_end}"
    )
    clinvar_same_pos_df = convert_vcf_gen_to_df(clinvar_same_pos)
    if not clinvar_same_pos_df.empty:
        clinvar_same_pos_formatted = format_spliceai(clinvar_same_pos_df)
        clinvar_same_pos_spliceai = clinvar_same_pos_formatted[
            clinvar_same_pos_formatted.SpliceAI_max >= var_spliceai
        ]
        ClinVar_same_pos = create_ClinVar(
            clinvar_same_pos_spliceai, ClinVar_Type.SAME_NUCLEOTIDE
        )
    else:
        ClinVar_same_pos = create_ClinVar(pd.DataFrame(), ClinVar_Type.SAME_NUCLEOTIDE)

    ### Check ClinVar for pathogenic variant in same / closest splice site
    (start_splice_site, end_splice_site) = find_corresponding_splice_site(
        affected_transcript, ref_transcript, variant
    )
    clinvar_splice_site = clinvar(
        f"{variant.chr}:{start_splice_site}-{end_splice_site}"
    )
    clinvar_splice_site_df = convert_vcf_gen_to_df(clinvar_splice_site)
    if not clinvar_splice_site_df.empty:
        clinvar_splice_site_formatted = format_spliceai(clinvar_splice_site_df)
        clinvar_splice_site_spliceai = clinvar_splice_site_formatted[
            clinvar_splice_site_formatted.SpliceAI_max >= var_spliceai
        ]
        ClinVar_splice_site = create_ClinVar(
            clinvar_splice_site_spliceai, ClinVar_Type.SAME_SPLICE_SITE
        )
    else:
        ClinVar_splice_site = create_ClinVar(
            pd.DataFrame(), ClinVar_Type.SAME_SPLICE_SITE
        )
    return {
        ClinVar_Type.SAME_NUCLEOTIDE: ClinVar_same_pos,
        ClinVar_Type.SAME_SPLICE_SITE: ClinVar_splice_site,
    }


def get_annotate_clinvar_spliceai_splicing(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for clinvar annotation and needed classification_information objects
    """
    return (
        annotate_clinvar_spliceai_splicing,
        (
            class_info.VARIANT,
            class_info.TRANSCRIPT,
            class_info.CLINVAR_PATH_SPLICEAI,
            class_info.VARIANT_PREDICTION,
        ),
    )
