#!/usr/bin/env python3

from typing import Callable
from dataclasses import dataclass

from refactoring.transcript_annotated import TranscriptInfo_annot
from clinvar_utils import ClinVar
from refactoring.variant import (
    VariantInfo,
    PopulationDatabases,
    AffectedRegion,
)


@dataclass
class Variant_annotated(VariantInfo):
    """
    Annotation of variant class
    """

    variant_info: VariantInfo
    transcript_info_annot: list[TranscriptInfo_annot]
    prediction_tools: dict
    same_aa_change_clinvar: ClinVar
    diff_aa_change_clinvar: ClinVar
    same_nuleotide_change: ClinVar
    same_splice_site: ClinVar
    gnomad: PopulationDatabases
    flossies: PopulationDatabases
    affected_region: AffectedRegion

    @classmethod
    def annotate(cls, annotations: list[Callable], variant: VariantInfo) -> Variant_annotated:
        for annotation in annotations:
            variant = annotation(variant)
        return variant
