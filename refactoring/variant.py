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
    exon_skipping: bool
    exons_affected: List[int]
    types_exon_skipped: List[str]
    is_NMD: bool
    NMD_affected_exons: List[int]
    NMD_affected_transcript_disease_relevant: bool
    truncated_exon_relevant: bool
    is_reading_frame_preserved: bool
    exists_alternative_start_codon: bool
    pathogenic_variant_between_start_and_stop: bool

@dataclass
class PredictionTools:
    name: str
    value: float

@dataclass
class PopulationDatabases:
    name: str
    frequency: float

@dataclass
class AffectedRegion:
    repetitive_region: bool
    critical_region: bool

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
class Variant:
    variant_info: VariantInfo
    transcript_info: List[TranscriptInfo]
    prediction_tools: List[PredictionTools]
    population_data: List[PopulationDatabases]
    affected_region: AffectedRegion
    thresholds: dict

@dataclass
class RuleResult:
    name: str
    status: bool
    strength: str
    comment: str
