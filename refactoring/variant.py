#!/usr/bin/env python3

from dataclasses import dataclass
from typing import Literal, Optional
import hgvs.posedit
import hgvs.parser

hgvs_parser = hgvs.parser.Parser()


VarType = Literal[
    "transcript_ablation",
    "splice_donor_variant",
    "splice_donor",
    "splice_acceptor_variant",
    "splice_acceptor",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature truncation",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
]


@dataclass
class TranscriptInfo:
    """
    Class containing general transcript information
    """

    transcript_id: str
    var_type: list[VarType]
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
