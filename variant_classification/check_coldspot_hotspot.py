#!/usr/bin/env python3

import pathlib

from collections.abc import Callable

import pyensembl

from pybedtools import BedTool

from ensembl import ensembl
from information import Classification_Info, Info
from variant import TranscriptInfo, VariantInfo
from utils import create_bed_line


def check_variant_intersection_with_bed(
    variant_hotspot_annotation_path: pathlib.Path,
    variant: VariantInfo,
    transcript: list[TranscriptInfo],
) -> bool:
    """
    Check hotspot annotation file for location of variant
    """
    # The reference transcript is only needed for the strand, therefore the specific transcript does not matter
    ref_transcript = ensembl.transcript_by_id(transcript[0].transcript_id)
    var_start = variant.genomic_start
    var_end = variant.genomic_end
    is_in_hotspot = check_intersection_with_bed_no_strand(
        variant, var_start, var_end, ref_transcript, variant_hotspot_annotation_path
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
    annotation_hits = bed.all_hits(variant_interval)
    if len(annotation_hits) > 0:
        return True
    return False
