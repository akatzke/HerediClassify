#!/usr/bin/env python3

from dataclasses import dataclass
from typing import List, Optional
import hgvs.posedit
import hgvs.parser

hgvs_parser = hgvs.parser.Parser()


@dataclass
class TranscriptInfo:
    transcript_id: str
    var_type: str
    var_hgvs: hgvs.posedit.PosEdit
    var_start: int
    var_stop: int
    var_protein: Optional[str] = None
    exon: Optional[int] = None
    intron: Optional[int] = None
    skipped_exon_start: Optional[int] = None
    skipped_exon_end: Optional[int] = None
    var_seq: Optional[str] = None
    diff_len: Optional[int] = None
    diff_len_percent: Optional[float] = None
    diff_len_protein_percent: Optional[float] = None
    len_change_in_repetitive_region: Optional[bool] = False
    are_exons_skipped: Optional[bool] = False
    exons_skipped: Optional[List[int]] = None
    stop_codon_exon_skipped: Optional[bool] = False
    start_codon_exon_skipped: Optional[bool] = False
    coding_exon_skipped: Optional[bool] = False
    is_NMD: Optional[bool] = False
    NMD_affected_exons: Optional[List[int]] = None
    transcript_disease_relevant: Optional[bool] = False
    truncated_exon_relevant: Optional[bool] = False
    is_reading_frame_preserved: Optional[bool] = False
    exists_alternative_start_codon: Optional[bool] = False
    pathogenic_variant_between_start_and_stop: Optional[bool] = False


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


@dataclass
class VariantInfo:
    gene_name: str
    var_type: list
    chr: str
    genomic_start: int
    genomic_end: int
    var_id: str
    var_ref: str
    var_obs: str

    def to_string(self) -> str:
        return f"{self.chr}:{self.genomic_start}-{self.genomic_end}{self.var_ref}>{self.var_obs}"


@dataclass
class ClinVar:
    same_amino_acid_change_pathogenic: bool
    same_amino_acid_change_pathogenic_list: List[str]
    different_amino_acid_change_in_same_position_pathogenic: bool
    different_amino_acid_change_in_same_position_pathogenic_list: List[str]


@dataclass
class Variant:
    variant_info: VariantInfo
    transcript_info: List[TranscriptInfo]
    prediction_tools: dict
    gnomad: PopulationDatabases
    flossies: PopulationDatabases
    affected_region: AffectedRegion
    clinvar: ClinVar
    thresholds: dict


@dataclass
class RuleResult:
    name: str
    status: bool
    strength: str
    comment: str
