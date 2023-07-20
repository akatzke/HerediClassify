#!/usr/bin/env python3

from dataclasses import dataclass
from typing import List


@dataclass
class TranscriptInfo:
    transcript_id: str
    var_type: str
    exon: int
    var_hgvs: hgvs_parser
    var_start: int
    var_stop: int
    var_seq: str
    var_protein: str
    diff_len_percent: float
    diff_len_protein_percent: float
    len_change_in_repetitive_region: bool
    exon_skipping: bool
    exons_affected: List[int]
    types_exon_skipped: List[str]
    is_NMD: bool
    NMD_affected_exons: List[int]
    transcript_disease_relevant: bool
    truncated_exon_relevant: bool
    is_reading_frame_preserved: bool
    exists_alternative_start_codon: bool
    pathogenic_variant_between_start_and_stop: bool


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
