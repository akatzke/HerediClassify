#!/usr/bin/env python3

import logging
import math
import pyensembl

from genotoscope_assess_NMD import (
    search_termination_codon,
    extract_codons,
)

logger = logging.getLogger("GenOtoScope_Classify.PVS1.prot_len_diff")


def calculate_prot_len_diff(
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
    diff_len: int,
) -> tuple[float, int]:
    """
    Calculate difference in portein length caused by variant in percent
    """
    ref_prot_len = len(ref_transcript.protein_sequence)
    codon_position_ptc = get_position_ptc(ref_transcript, var_coding_seq.upper())
    # Substract one from postion of ptc, as ptc does not code for amino acid
    rel_len_protein = abs(codon_position_ptc - 1) / ref_prot_len
    diff_len_protein_percent = 1 - rel_len_protein
    corrected_codon_position_ptc = correct_position_ptc_for_indels(
        diff_len, codon_position_ptc
    )
    return diff_len_protein_percent, corrected_codon_position_ptc


def calculate_prot_len_diff_start_loss(
    ref_transcript: pyensembl.transcript.Transcript,
    alt_start_codon: list[int],
) -> float:
    """
    Differences in protein length caused by start loss variants
    """
    ref_prot_len = len(ref_transcript.protein_sequence)
    codon_position = alt_start_codon[0] / 3
    if alt_start_codon[0] % 3 != 0:
        raise ValueError(
            f"The position of the first codon of the alternative start codon is not divisible by three. Start codon {alt_start_codon}."
        )
    var_prot_len = len(ref_transcript.protein_sequence[int(codon_position) :])
    rel_len_protein = var_prot_len / ref_prot_len
    diff_len_protein_percent = 1 - rel_len_protein
    return diff_len_protein_percent


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


def correct_position_ptc_for_indels(
    diff_len: int,
    codon_positon_ptc: int,
) -> int:
    """
    Correct position of PTC in case of insertion, deletion, indels
    Adjust the PTC codon by codons added/substracted by the variant
    Non-complete codons are rounded up, as the start/end of the codon is still affected
    """
    if diff_len == 0:
        return codon_positon_ptc
    elif diff_len > 0:
        # In case of an insertion the number of codons inserted needs to be substraced to get the position of the PTC in the original transcript
        codon_shift = math.floor(diff_len / 3) * -1
        return codon_positon_ptc + codon_shift
    elif diff_len < 0:
        # In case of a deletion the number of codons deleted needs to be added to get the position of the PTC in the original transcript
        codon_shift = math.floor(abs(diff_len) / 3)
        return codon_positon_ptc + codon_shift
    else:
        raise ValueError("Error whilst correcting PTC position.")
