#!/usr/bin/env python3

from typing import Callable, Optional, Literal
import pathlib
from attr import dataclass
from cyvcf2 import Variant
from refactoring.variant import (
    VariantInfo,
    TranscriptInfo,
    PopulationDatabases,
    AffectedRegion,
)

ClinVar_Type = Literal[
    "same_aa_change", "diff_aa_change", "same_nucleotide", "same_splice_site", "region"
]

ClinVar_Status = Literal["Pathogenic", "Likely_pathogenic"]


@dataclass
class ClinVar:
    pathogenic: bool
    type: ClinVar_Type
    highest_classification: Optional[ClinVar_Status]
    ids: Optional[list[str]]


@dataclass
class Variant_annotated(VariantInfo):
    """
    Annotation of variant class
    """

    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
    prediction_tools: dict
    same_aa_change_clinvar: ClinVar
    diff_aa_change_clinvar: ClinVar
    same_nuleotide_change: ClinVar
    same_splice_site: ClinVar
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
