#!/usr/bin/env python3

import logging
import pyensembl

from variant_classification.genotoscope_assess_NMD import search_termination_codon, extract_codons

logger = logging.getLogger("GenOtoScope_Classify.PVS1.prot_len_diff")


def calculate_prot_len_diff(
    ref_transcript: pyensembl.transcript.Transcript,
    var_coding_seq: str,
) -> float:
    """
    Calculate difference in portein length caused by variant in percent
    """
    ref_prot_len = len(ref_transcript.protein_sequence)
    codon_position_ptc = get_position_ptc(ref_transcript, var_coding_seq)
    # Substract one from postion of ptc, as ptc does not code for amino acid
    rel_len_protein = abs(codon_position_ptc - 1) / ref_prot_len
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
