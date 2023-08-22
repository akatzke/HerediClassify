#!/usr/bin/env python3

from typing import Callable, Optional
import pathlib
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
    path_clinvar: pathlib.Path
    filter: Optional[list[str]]


@dataclass
class ClinVar_exonic(ClinVar):
    same_aa_change_pathogenic: bool
    same_aa_change_highest_classification: Optional[str]
    same_aa_change_ids: Optional[list[str]]
    diff_aa_change_in_pathogenic: bool
    diff_aa_change_highest_classification: Optional[str]
    diff_aa_change_ids: Optional[list[str]]


@dataclass
class ClinVar_splice(ClinVar):
    same_nucleotide_change_pathogenic: bool
    same_nucleotide_change_highest_classification: Optional[str]
    same_nucleotide_change_ids: Optional[list[str]]
    same_splice_site_pathogenic: bool
    same_splice_site_highest_classification: Optional[str]
    same_splice_site_change_ids: Optional[list[str]]


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
