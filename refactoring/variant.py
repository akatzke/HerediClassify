#!/usr/bin/env python3

from dataclasses import dataclass, field
from typing import Optional
import hgvs.posedit
import hgvs.parser
import pyensembl

from refactoring.genotoscope_assess_NMD import assess_NMD
from refactoring.genotoscope_construct_variant_sequence import (
    construct_variant_coding_seq,
)
from refactoring.genotoscope_exon_skipping import assess_exon_skipping

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
    var_protein: Optional[str] = None
    exon: Optional[int] = None
    intron: Optional[int] = None
    diff_len: Optional[float] = None
    diff_len_percent: Optional[float] = None
    diff_len_protein_percent: Optional[float] = None
    exons_skipped: list[int] = field(default_factory=list)
    skipped_exon_start: Optional[int] = None
    skipped_exon_end: Optional[int] = None
    stop_codon_exon_skipped: bool = False
    start_codon_exon_skipped: bool = False
    are_exons_skipped: bool = False
    coding_exon_skipped: bool = False
    var_seq: Optional[str] = None
    len_change_in_repetitive_region: bool = False
    is_NMD: bool = False
    NMD_affected_exons: list[int] = field(default_factory=list)
    transcript_disease_relevant: bool = False
    truncated_exon_relevant: bool = False
    is_reading_frame_preserved: bool = False
    exists_alternative_start_codon: bool = False
    pathogenic_variant_between_start_and_stop: bool = False


@dataclass
class TranscriptInfo_exonic(TranscriptInfo):
    """
    Class containing exonic variant specific annotation
    """

    ref_transcript: pyensembl.transcript.Transcript
    diff_len: int = 0
    diff_len_percent: float = 0
    diff_len_protein_percent: Optional[float] = None
    len_change_in_repetitive_region: bool = False

    def __post_init__(self, variant: VariantInfo):
        self.ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            self.transcript_id
        )
        (self.var_seq, self.diff_len) = construct_variant_coding_seq(self, variant)
        self.diff_len_percent = 1
        self.diff_len_protein_percent = 1
        (self.is_NMD, self.NMD_affected_exons) = assess_NMD(self, variant)


@dataclass
class TranscriptInfo_intronic(TranscriptInfo):
    """
    Class containing intronic variant specific annotations
    """

    ref_transcript: pyensembl.transcript.Transcript
    diff_len: float
    diff_len_percent: float
    diff_len_protein_percent: float
    exons_skipped: list[int]
    skipped_exon_start: int
    skipped_exon_end: int
    stop_codon_exon_skipped: bool = False
    start_codon_exon_skipped: bool = False
    are_exons_skipped: bool = False
    coding_exon_skipped: bool = False
    var_seq: str = ""
    len_change_in_repetitive_region: bool = False
    is_NMD: bool = False
    NMD_affected_exons: list[int] = field(default_factory=list)
    transcript_disease_relevant: bool = False
    truncated_exon_relevant: bool = False
    is_reading_frame_preserved: bool = False
    exists_alternative_start_codon: bool = False
    pathogenic_variant_between_start_and_stop: bool = False

    def __post_init__(self, variant: VariantInfo):
        self.ref_transcript = pyensembl.EnsemblRelease(75).transcript_by_id(
            self.transcript_id
        )
        (
            self.exons_skipped,
            self.are_exons_skipped,
            self.skipped_exon_start,
            self.skipped_exon_end,
            self.start_codon_exon_skipped,
            self.stop_codon_exon_skipped,
            self.coding_exon_skipped,
        ) = assess_exon_skipping(self, variant)
        (self.var_seq, self.diff_len) = construct_variant_coding_seq(self, variant)
        self.diff_len_percent = 1
        self.diff_len_protein_percent = 1
        (self.is_NMD, self.NMD_affected_exons) = assess_NMD(self, variant)


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
    same_amino_acid_change_pathogenic_list: list[str]
    different_amino_acid_change_in_same_position_pathogenic: bool
    different_amino_acid_change_in_same_position_pathogenic_list: list[str]


@dataclass
class Variant:
    variant_info: VariantInfo
    transcript_info: list[TranscriptInfo]
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
