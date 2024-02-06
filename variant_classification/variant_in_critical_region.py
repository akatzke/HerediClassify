#!/usr/bin/env python3

import pathlib
import logging

import pyensembl

from variant import VariantInfo
from utils import check_intersection_with_bed


logger = logging.getLogger("GenOtoScope_Classify.protein_len_diff_repetitive_region")


def check_variant_in_critical_region_exon(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    NMD_affected_exons: list[dict],
    path_bed: pathlib.Path,
) -> tuple[bool, str]:
    """
    Check if any of the NMD affected are in critical region
    """
    if NMD_affected_exons:
        are_exons_in_critical_region = []
        comment = ""
        for exon in NMD_affected_exons:
            gen_start = exon["exon_start"]
            gen_end = exon["exon_end"]
            is_exon_in_repetitive_region, comment = check_intersection_with_bed(
                variant, gen_start, gen_end, ref_transcript, path_bed
            )
            are_exons_in_critical_region.append(is_exon_in_repetitive_region)
        if any(are_exons_in_critical_region):
            comment = f"Truncated exon overlaps the following clinically significant domains: {comment}."
            return True, comment
        else:
            comment = (
                "Truncated exon does not overlap clinically significant protein domain."
            )
            return False, comment
    else:
        comment = "No truncated region, there truncated region not disease relevant."
        return False, comment
