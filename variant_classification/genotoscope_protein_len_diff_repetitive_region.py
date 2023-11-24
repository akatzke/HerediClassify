#!/usr/bin/env python3

import pathlib
import logging

import pyensembl

from variant import VariantInfo
from utils import check_intersection_with_bed


logger = logging.getLogger("GenOtoScope_Classify.protein_len_diff_repetitive_region")


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
