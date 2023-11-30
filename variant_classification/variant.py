#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Optional
import hgvs.posedit

from var_type import VARTYPE


@dataclass
class TranscriptInfo:
    """
    Class containing general transcript information
    """

    transcript_id: str
    var_type: list[VARTYPE]
    var_hgvs: hgvs.posedit.PosEdit
    var_start: int
    var_stop: int
    var_protein: Optional[str]
    exon: Optional[int]
    intron: Optional[int]


@dataclass
class PopulationDatabases:
    name: str
    frequency: Optional[float]
    count: Optional[float]

    def __post_init__(self):
        "Ensure that"
        if self.frequency is None and self.count is None:
            raise ValueError(
                f"Both frequency and count of population database {self.name} were not given. Please give either frequency of count."
            )


@dataclass
class PopulationDatabases_gnomAD(PopulationDatabases):
    popmax: str
    popmax_frequency: float
    popmax_allele_count: int


@dataclass
class AffectedRegion:
    critical_region: Optional[bool] = None
    cancer_hotspot: Optional[bool] = None
    cold_spot: Optional[bool] = None


@dataclass
class VariantInfo:
    chr: str
    genomic_start: int
    genomic_end: int
    gene_name: str
    var_type: list[VARTYPE]
    var_ref: str
    var_obs: str

    def to_string(self) -> str:
        return f"{self.chr}:{self.genomic_start}-{self.genomic_end}{self.var_ref}>{self.var_obs}"


@dataclass
class MultifactorialLikelihood:
    prior: Optional[float] = None
    multifactorial_likelihood: Optional[float] = None
    co_occurrence: Optional[float] = None
    co_segregation: Optional[float] = None


@dataclass
class FunctionalData:
    performed: bool
    pathogenic: bool
    benign: bool


@dataclass
class Variant:
    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
    prediction_tools: Optional[dict[str, float]] = None
    gnomad: Optional[PopulationDatabases_gnomAD] = None
    flossies: Optional[PopulationDatabases] = None
    affected_region: Optional[AffectedRegion] = None
    multifactorial_likelihood: Optional[MultifactorialLikelihood] = None
    functional_assay: Optional[FunctionalData] = None
    splicing_assay: Optional[FunctionalData] = None
