#!/usr/bin/env python3

from typing import Callable, Optional
from dataclasses import dataclass, field
from cyvcf2 import Variant
from enum import Enum
from refactoring.transcript_annotated import TranscriptInfo_annot

from refactoring.variant import (
    VariantInfo,
    TranscriptInfo,
    PopulationDatabases,
    AffectedRegion,
)

class CLINVAR_TYPE(Enum):
    SAME_AA_CHANGE = "same_aa_change"
    DIFF_AA_CHANGE = "diff_aa_change"
    SAME_NUCLEOTIDE = "same_nucleotide"
    SAME_SPLICE_SITE = "same_splice_site"
    REGION = "region"

class CLINVAR_STATUS(Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"

@dataclass
class ClinVar:
    pathogenic: bool
    type: CLINVAR_TYPE
    highest_classification: Optional[CLINVAR_STATUS]
    ids: list[str] = field(default_factory=list)
    associated_ids: list[str] = field(default_factory=list)



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
    def annotate(cls, annotations: list[Callable], variant: Variant) -> Variant:
        for annotation in annotations:
            variant = annotation(variant)
        return variant
