#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Literal, Optional
import hgvs.posedit
import hgvs.parser

hgvs_parser = hgvs.parser.Parser()


@dataclass
class TranscriptInfo:
    """
    Class containing general transcript information
    """

    transcript_id: str
    var_type: str
    var_hgvs: hgvs.posedit.PosEdit
    var_start: int
    var_stop: int
    var_protein: Optional[str]
    exon: Optional[int]
    intron: Optional[int]


@dataclass
class PredictionTools:
    name: str
    type: str
    value: float


@dataclass
class PopulationDatabases:
    name: str
    frequency: float


@dataclass
class AffectedRegion:
    repetitive_region: bool
    critical_region: bool
    critical_region_type: str


VarType = Literal["splice_donor", "splice_acceptor"]


@dataclass
class VariantInfo:
    gene_name: str
    var_type: list[VarType]
    chr: str
    genomic_start: int
    genomic_end: int
    var_id: str
    var_ref: str
    var_obs: str

    def to_string(self) -> str:
        return f"{self.chr}:{self.genomic_start}-{self.genomic_end}{self.var_ref}>{self.var_obs}"


@dataclass
class Variant_import:
    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
    prediction_tools: dict
    gnomad: PopulationDatabases
    flossies: PopulationDatabases
    affected_region: AffectedRegion


@dataclass
class RuleResult:
    name: str
    status: bool
    strength: str
    comment: str
