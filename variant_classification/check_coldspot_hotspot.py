#!/usr/bin/env python3

import pathlib

from collections.abc import Callable

import pyensembl
import logging

from pybedtools import BedTool

from ensembl import ensembl
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from utils import create_bed_line

logger = logging.getLogger("GenOtoScope_Classify.check_coldspot_hotspot")


def check_variant_intersection_with_bed(
    variant_hotspot_annotation_path: pathlib.Path,
    variant: VariantInfo,
    transcript: list[TranscriptInfo],
) -> bool:
    """
    Check hotspot annotation file for location of variant
    """
    # The reference transcript is only needed for the strand, therefore the specific transcript does not matter
    strand = get_variant_strand(transcript, variant)
    var_start = variant.genomic_start
    var_end = variant.genomic_end
    is_in_hotspot = check_intersection_with_bed_no_strand(
        variant, var_start, var_end, strand, variant_hotspot_annotation_path
    )
    return is_in_hotspot


def get_check_hotspot(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for hotspot annotation and needed classification_information objects
    """
    return (
        check_variant_intersection_with_bed,
        (
            class_info.VARIANT_HOTSPOT_ANNOTATION_PATH,
            class_info.VARIANT,
            class_info.TRANSCRIPT,
        ),
    )


def get_check_coldspot(
    class_info: Classification_Info,
) -> tuple[Callable, tuple[Info, ...]]:
    """
    Get function for hotspot annotation and needed classification_information objects
    """
    return (
        check_variant_intersection_with_bed,
        (
            class_info.VARIANT_COLDSPOT_ANNOTATION_PATH,
            class_info.VARIANT,
            class_info.TRANSCRIPT,
        ),
    )


def check_intersection_with_bed_no_strand(
    variant: VariantInfo,
    gen_start: int,
    gen_end: int,
    strand: str,
    path_bed: pathlib.Path,
) -> bool:
    """
    Check if variant overlaps region in given bed file without checking for strand
    """
    variant_interval = BedTool(
        create_bed_line(variant, gen_start, gen_end, strand),
        from_string=True,
    )[0]
    bed = BedTool(path_bed).sort()
    annotation_hits = bed.all_hits(variant_interval)
    if len(annotation_hits) > 0:
        return True
    return False


def get_variant_strand(transcripts: list[TranscriptInfo], variant: VariantInfo) -> str:
    """
    Try to get strand for gene
    """
    try:
        # Try getting strand from transcript
        ref_transcript = ensembl.transcript_by_id(transcripts[0].transcript_id)
        strand = ref_transcript.strand
        return strand
    except Exception:
        # In case getting strand from transcript fails, try getting strand from gene
        try:
            genes = ensembl.genes_by_name(variant.gene_name)
            for gene in genes:
                # Check if variant location matches genen location
                if (
                    gene.contig in gene.chr
                    and (
                        variant.genomic_start >= gene.start
                        and gene.end >= variant.genomic_start
                    )
                    and (
                        gene.end >= variant.genomic_end
                        and variant.genomic_end >= gene.start
                    )
                ):
                    return gene.strand
            raise ValueError(
                f"None of the genes identified for {variant.gene_name} could be matched to the variant location {variant.chr}:{variant.genomic_start}-{variant.genomic_end}."
            )
        except Exception:
            raise ValueError(
                "The reconstruction of the variant strand failed. Please check."
            )
