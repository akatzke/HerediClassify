#!/usr/bin/env python3

import pathlib

import pyensembl

from pybedtools import BedTool
from variant import VariantInfo


def check_bed_intersect_start_loss(
    variant: VariantInfo,
    ref_transcript: pyensembl.transcript.Transcript,
    alt_start_codon: list[int],
    path_rep_uniprot: pathlib.Path,
) -> bool:
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


def check_intersection_with_bed(
    variant: VariantInfo,
    gen_start: int,
    gen_end: int,
    ref_transcript: pyensembl.transcript.Transcript,
    path_bed: pathlib.Path,
) -> bool:
    """
    Check if variant overlaps UniProt annotated repetitive region
    """
    variant_interval = BedTool(
        create_bed_line(variant, gen_start, gen_end, ref_transcript.strand),
        from_string=True,
    )[0]
    bed = BedTool(path_bed).sort()
    annotation_hits = bed.all_hits(variant_interval, same_strand=True)
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
