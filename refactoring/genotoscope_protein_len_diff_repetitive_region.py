#!/usr/bin/env python3

import logging
import pyensembl
import pathlib

from pybedtools import BedTool

from refactoring.variant import TranscriptInfo, VariantInfo

logger = logging.getLogger("GenOtoScope_Classify.protein_len_diff_repetitive_region")


def check_prot_len_change_in_repetitive_region(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    path_rep_uniprot: pathlib.Path,
) -> bool:
    """
    Check if variant overlaps UniProt annotated repetitive region
    """
    variant_interval = BedTool(
        create_bed_line(variant, ref_transcript.strand), from_string=True
    )[0]
    repeats = BedTool(path_rep_uniprot).sort()
    annotation_hits = repeats.all_hits(variant_interval, same_strand=True)
    if len(annotation_hits) > 0:
        return True
    else:
        return False


def create_bed_line(variant: VariantInfo, transcript_strand: str) -> str:
    """
    Create bed line to represent variant info
    Credits: http://daler.github.io/pybedtools/intervals.html
    """
    bed_line = " ".join(
        [
            "chr" + str(variant.chr),
            str(variant.genomic_start),
            str(variant.genomic_end),
            variant.gene_name,
            ".",
            transcript_strand,
        ]
    )
    return bed_line
