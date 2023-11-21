#!/usr/bin/env python3

import pathlib

import pyensembl

from pybedtools import BedTool
from variant_classification.variant import VariantInfo


def check_intersection_with_bed(
    variant: VariantInfo,
    gen_start: int,
    gen_end: int,
    ref_transcript: pyensembl.transcript.Transcript,
    path_uniprot_rep: pathlib.Path,
) -> bool:
    """
    Check if variant overlaps UniProt annotated repetitive region
    """
    variant_interval = BedTool(
        create_bed_line(variant, gen_start, gen_end, ref_transcript.strand),
        from_string=True,
    )[0]
    repeats = BedTool(path_uniprot_rep).sort()
    annotation_hits = repeats.all_hits(variant_interval, same_strand=True)
    if len(annotation_hits) > 0:
        return True
    return False


def create_bed_line(
    variant: VariantInfo, gen_start: int, gen_end: int, transcript_strand: str
) -> str:
    """
    Create bed line to represent variant info
    Credits: http://daler.github.io/pybedtools/intervals.html
    """
    bed_line = " ".join(
        [
            "chr" + str(variant.chr),
            str(gen_start),
            str(gen_end),
            variant.gene_name,
            ".",
            transcript_strand,
        ]
    )
    return bed_line
