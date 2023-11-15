#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Optional
import hgvs.posedit
import hgvs.parser

from refactoring.var_type import VARTYPE

hgvs_parser = hgvs.parser.Parser()


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
    frequency: float


@dataclass
class PopulationDatabases_gnomAD(PopulationDatabases):
    allele_count: int
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
class Variant:
    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
    prediction_tools: Optional[dict[str, float]] = None
    gnomad: Optional[PopulationDatabases_gnomAD] = None
    flossies: Optional[PopulationDatabases] = None
    affected_region: Optional[AffectedRegion] = None
