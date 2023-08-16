#!/usr/bin/env python3

from typing import Callable
from attr import dataclass
from cyvcf2 import Variant
from refactoring.variant import (
    VariantInfo,
    TranscriptInfo,
    PopulationDatabases,
    AffectedRegion,
)


@dataclass
class ClinVar:
    same_nucleotide_change_pathogenic: bool
    matching_clinvar_entries: list[str]
    matching_clinvar_entries_highest_classification: str


@dataclass
class ClinVar_exonic(ClinVar):
    same_amino_acid_change_pathogenic: bool
    different_amino_acid_change_in_same_position_pathogenic: bool


@dataclass
class ClinVar_splice(ClinVar):
    same_splice_site_pathogenic: bool


@dataclass
class Variant_annotated(VariantInfo):
    """
    Annotation of variant class
    """

    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
    prediction_tools: dict
    gnomad: PopulationDatabases
    flossies: PopulationDatabases
    affected_region: AffectedRegion
    clinvar: ClinVar
    thresholds: dict

    @classmethod
    def annotate(cls, annotations: list[Callable], variant: Variant) -> Variant:
        for annotation in annotations:
            variant = annotation(variant)
        return variant
