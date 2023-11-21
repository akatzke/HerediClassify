#!/usr/bin/env python3

import pathlib
import logging

import pyensembl

from variant_classification.variant import VariantInfo
from variant_classification.utils import check_intersection_with_bed


logger = logging.getLogger("GenOtoScope_Classify.protein_len_diff_repetitive_region")


def check_prot_len_change_in_repetitive_region_start_loss(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    alt_start_codon: list[int],
    path_rep_uniprot: pathlib.Path,
):
    """
    Assess if prot_len_change caused by alternative start codon is in repetitive region
    """
    ref_start_codon = ref_transcript.start_codon_positions
    if ref_transcript.strand == "-":
        gen_start = min(alt_start_codon)
        gen_end = max(ref_start_codon)
    else:
        gen_start = min(ref_start_codon)
        gen_end = max(alt_start_codon)
    prot_len_in_repetitive_region = check_intersection_with_bed(
        variant, gen_start, gen_end, ref_transcript, path_rep_uniprot
    )
    return prot_len_in_repetitive_region


def check_prot_len_change_in_repetitive_region_exon(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    NMD_affected_exons: list[dict],
    path_rep_uniprot: pathlib.Path,
):
    if NMD_affected_exons:
        are_exons_in_repetitive_region = []
        for exon in NMD_affected_exons:
            gen_start = exon["exon_start"]
            gen_end = exon["exon_end"]
            is_exon_in_repetitive_region = check_intersection_with_bed(
                variant, gen_start, gen_end, ref_transcript, path_rep_uniprot
            )
            are_exons_in_repetitive_region.append(is_exon_in_repetitive_region)
        if all(are_exons_in_repetitive_region):
            return True
        else:
            return False
    else:
        return False
