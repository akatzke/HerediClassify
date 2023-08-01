#!/usr/bin/env python3

import logging
import pyensembl

from refactoring.genotoscope_assess_NMD import search_termination_codon, extract_codons
from refactoring.variant import TranscriptInfo, VariantInfo

logger = logging.getLogger("GenOtoScope_Classify.PVS1.prot_len_diff")


def calculate_prot_len_diff_exonic(
    variant: VariantInfo,
    transcript: TranscriptInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    is_reading_frame_preserved: bool,
    var_coding_seq: str,
    diff_len: int,
) -> float:
    """
    Calculate difference in portein length caused by variant in percent
    """
    ref_prot_len = len(ref_transcript.protein_sequence)
    if is_reading_frame_preserved:
        if "inframe" in transcript.var_type:
            prot_len_diff = diff_len / 3
            return abs(ref_prot_len - prot_len_diff) / ref_prot_len
        elif transcript.var_type == "stop_gain":
            codon_position_ptc = get_position_ptc(ref_transcript, var_coding_seq)
            return abs(ref_prot_len - codon_position_ptc) / ref_prot_len
        elif transcript.var_type == "stop_lost":
            closest_stop_codon_idx, closest_stop_codon = find_inframe_stop(
                ref_transcript.three_prime_utr_sequence
            )
            closest_stop_codon_dist = (closest_stop_codon_idx + 1) * 3
            prot_len_diff = calculate_observed_protein_length(
                variant, closest_stop_codon_dist
            )
            return abs(ref_prot_len - prot_len_diff) / ref_prot_len
        else:
            return 0
    elif not is_reading_frame_preserved:
        codon_position_ptc = get_position_ptc(ref_transcript, var_coding_seq)
        return abs(ref_prot_len - codon_position_ptc) / ref_prot_len


def calculate_prot_len_diff_exonic(
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
) -> float:
    """
    Calculate difference in portein length caused by variant in percent
    """
    ref_prot_len = len(ref_transcript.protein_sequence)
    codon_position_ptc = get_position_ptc(ref_transcript, var_coding_seq)
    return abs(ref_prot_len - codon_position_ptc) / ref_prot_len


def get_position_ptc(
    ref_transcript: pyensembl.transcript.Transcript, var_coding_seq: str
) -> int:
    """
    Calculate the protein length of the observed coding sequence based on the (premature) termination codon
    Called calculate_prot_len_ptc in GenOtoScope
    """
    # after variant edit, the termination codon can be found even in the 3' UTR region
    # thus, search for the very first termination codon on the constructed observed coding sequence pluts the 3' UTR
    var_coding_3_utr = var_coding_seq + ref_transcript.three_prime_utr_sequence
    _, premature_term_codon_index = search_termination_codon(
        extract_codons(var_coding_3_utr), False
    )
    if premature_term_codon_index == -1:
        return 0
    else:
        return premature_term_codon_index


def calculate_observed_protein_length(
    variant_info: VariantInfo, new_stop_distance: int = 0
) -> float:
    if variant_info.var_type == "frameshift_deletion":
        return len(variant_info.var_ref) / 3
    elif variant_info.var_type == "frameshift_insertion":
        return len(variant_info.var_obs) / 3
    else:
        # stop lost
        return new_stop_distance / 3


def find_inframe_stop(codons: list[str]) -> tuple:
    """
    Find in-frame stop codon in the 3 prime UTR region

    Returns
    -------
    closest_stop_codon : str
        closest in-frame stop codon (TAA or TAG or TGA)
    closet_stop_codon_idx : int
        0-based index of closest in-frame stop codon
    """

    logger.debug("Find in-frame stop codon in 3' UTR region")
    stop_codon_indices = {"TAA": -1, "TAG": -1, "TGA": -1}
    closest_stop_codon_idx, closest_stop_codon = len(codons), None
    for stop_codon in stop_codon_indices:
        if stop_codon in codons:
            stop_codon_indices[stop_codon] = codons.index(stop_codon)
            if stop_codon_indices[stop_codon] < closest_stop_codon_idx:
                closest_stop_codon_idx = stop_codon_indices[stop_codon]
                closest_stop_codon = stop_codon
    return closest_stop_codon_idx, closest_stop_codon
