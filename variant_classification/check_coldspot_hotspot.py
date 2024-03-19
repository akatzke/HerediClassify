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
        genes = ensembl.genes_by_name(variant.gene_name)
        # If there is only one gene matching, assume it's the correct one and return strand
        if len(genes) == 1:
            return genes[0].strand
        # If more than one gene was found, check if any of the genes is on the same chromosome as the variant
        genes_same_chr = []
        for gene in genes:
            if gene.contig in variant.chr:
                genes_same_chr.append(gene)
        if not len(genes_same_chr):
            raise ValueError(
                "The reconstruction of the varaint strand failed. None of the matching genes is located on the same strand as the variant."
            )
        elif len(genes_same_chr) == 1:
            return genes_same_chr[0].strand
        else:
            # If more than one gene is on the same chromosome, check the location, pick the closest gene
            closest_gene = find_gene_closest_to_variant(variant, genes_same_chr)
            return closest_gene.strand


def find_gene_closest_to_variant(
    variant: VariantInfo, genes: list[pyensembl.gene.Gene]
) -> pyensembl.gene.Gene:
    """
    Find the gene located clostest to the variant
    """
    distance_dict = {}
    for gene in genes:
        dist_to_start = min(
            abs(variant.genomic_start - gene.start),
            abs(variant.genomic_end - gene.start),
        )
        dist_to_end = min(
            abs(variant.genomic_start - gene.end),
            abs(variant.genomic_end - gene.end),
        )
        dist = min(dist_to_start, dist_to_end)
        distance_dict[gene.id] = dist
    gene_id_closest = min(distance_dict, key=distance_dict.get)
    for gene in genes:
        if gene.id == gene_id_closest:
            return gene
