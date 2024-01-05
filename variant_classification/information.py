#!/usr/bin/env python3

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import TypeVar, Optional, Generic


ValueType = TypeVar("ValueType")


class Classification_Info_Groups(Enum):
    THRESHOLD_SINGLE = auto()
    THRESHOLD_PREDICTION = auto()
    THRESHOLD_MULT_STRENGTH = auto()
    CONFIG_ENTRY_STR = auto()
    PATH = auto()
    DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD = auto()


@dataclass(frozen=False)
class Info(Generic[ValueType]):
    name: str
    config_location: Optional[tuple[str, ...]] = field(default=None)
    compute_function: Optional[Callable[[], ValueType]] = field(
        init=False, default=None
    )
    value: Optional[ValueType] = field(init=False, default=None)
    group: Optional[Classification_Info_Groups] = field(default=None)
    optional: bool = False


@dataclass
class Classification_Info:
    ANNOTATED_TRANSCRIPT_LIST: Info
    ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec: Info
    VARIANT_CLINVAR: Info
    VARIANT_HOTSPOT: Info
    VARIANT_HOTSPOT_ANNOTATION: Info
    VARIANT_HOTSPOT_ANNOTATION_PATH: Info
    VARIANT_COLDSPOT_ANNOTATION: Info
    VARIANT_COLDSPOT_ANNOTATION_PATH: Info
    VARIANT_GNOMAD: Info
    VARIANT_FLOSSIES: Info
    VARIANT_PREDICTION: Info
    VARIANT_MULTIFACTORIAL_LIKELIHOOD: Info
    THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC: Info
    THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN: Info
    THRESHOLD_SPLICING_PREDICTION_PATHOGENIC: Info
    THRESHOLD_SPLICING_PREDICTION_BENIGN: Info
    THRESHOLD_LIKELIHOOD_BENIGN: Info
    THRESHOLD_LIKELIHOOD_PATHOGENIC: Info
    THRESHOLD_PM2: Info
    THRESHOLD_BA1: Info
    THRESHOLD_BA1_ABSOLUTE: Info
    THRESHOLD_BS1: Info
    THRESHOLD_BS1_ABSOLUTE: Info
    THRESHOLD_BS1_SUPPORTING: Info
    THRESHOLD_BS2: Info
    THRESHOLD_BS2_SUPPORTING: Info
    THRESHOLD_DIFF_LEN_PROT_PERCENT: Info
    THRESHOLD_NMD: Info
    POS_LAST_KNOWN_PATHO_PTC: Info
    VARIANT: Info
    TRANSCRIPT: Info
    CLINVAR_PATH: Info
    SIMILARITY_SCORE_PATH: Info
    SIMILARITY_SOCRE_DIRECTION: Info
    UNIPROT_REP_REGION_PATH: Info
    CRITICAL_REGION_PATH: Info
    DISEASE_IRRELEVANT_EXONS_PATH: Info
    SPLICE_SITE_TABLE_PATH: Info
    SPLICE_RESULT: Info
    SPLICE_SITE_TABLE_PM5_PATH: Info
    SPLICE_RESULT_PM5: Info
    EXON_PM5_PATH: Info
    PM5_RESULTS_PTC: Info
    FUNCTIONAL_ASSAY: Info
    SPLICING_ASSAY: Info

    def __init__(self):
        self.ANNOTATED_TRANSCRIPT_LIST = Info("annotated_transcript_list")
        self.ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec = Info(
            "annotated_transcript_list_acmg"
        )
        self.VARIANT_CLINVAR = Info("variant_clinvar")
        self.VARIANT_HOTSPOT = Info("variant_hotspot")
        self.VARIANT_HOTSPOT_ANNOTATION = Info("variant_hotspot_annotation")
        self.VARIANT_COLDSPOT_ANNOTATION = Info("variant_coldspot_annotation")
        self.VARIANT_COLDSPOT_ANNOTATION_PATH = Info(
            "variant_coldspot_annotation_path",
            config_location=("annotation_files", "critical_regions", "coldspot_region"),
            group=Classification_Info_Groups.PATH,
        )
        self.VARIANT_HOTSPOT_ANNOTATION_PATH = Info(
            "variant_hotspot_annotation_path",
            config_location=("annotation_files", "critical_regions", "hotspot_region"),
            group=Classification_Info_Groups.PATH,
        )
        self.VARIANT_GNOMAD = Info("variant_gnomad")
        self.VARIANT_FLOSSIES = Info("variant_flossies")
        self.VARIANT_PREDICTION = Info("variant_prediction")
        self.VARIANT_MULTIFACTORIAL_LIKELIHOOD = Info(
            "variant_multifactorial_likelihood"
        )
        self.THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN = Info(
            "prediction_pathogenicity_benign",
            config_location=(
                "prediction_tool_thresholds",
                "pathogenicity_prediction",
                "benign",
            ),
            group=Classification_Info_Groups.THRESHOLD_PREDICTION,
        )
        self.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC = Info(
            "prediction_pathogenicity_pathogenic",
            config_location=(
                "prediction_tool_thresholds",
                "pathogenicity_prediction",
                "pathogenic",
            ),
            group=Classification_Info_Groups.THRESHOLD_PREDICTION,
        )
        self.THRESHOLD_SPLICING_PREDICTION_BENIGN = Info(
            "prediction_splicing_benign",
            config_location=(
                "prediction_tool_thresholds",
                "splicing_prediction",
                "benign",
            ),
            group=Classification_Info_Groups.THRESHOLD_PREDICTION,
        )
        self.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC = Info(
            "prediction_splicing_pathogenic",
            config_location=(
                "prediction_tool_thresholds",
                "splicing_prediction",
                "pathogenic",
            ),
            group=Classification_Info_Groups.THRESHOLD_PREDICTION,
        )
        self.THRESHOLD_LIKELIHOOD_PATHOGENIC = Info(
            "threshold_likelihood_pathogenic",
            config_location=("likelihood_thresholds", "pathogenic"),
            group=Classification_Info_Groups.THRESHOLD_MULT_STRENGTH,
        )
        self.THRESHOLD_LIKELIHOOD_PATHOGENIC = Info(
            "threshold_likelihood_pathogenic",
            config_location=("likelihood_thresholds", "benign"),
            group=Classification_Info_Groups.THRESHOLD_MULT_STRENGTH,
        )
        self.THRESHOLD_PM2 = Info(
            "threshold_pm2",
            config_location=("allele_frequency_thresholds", "threshold_pm2"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BA1 = Info(
            "threshold_ba1",
            config_location=("allele_frequency_thresholds", "threshold_ba1"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BA1_ABSOLUTE = Info(
            "threshold_ba1_absolute",
            config_location=("allele_frequency_thresholds", "threshold_ba1_absolute"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BS1 = Info(
            "threshold_bs1",
            config_location=("allele_frequency_thresholds", "threshold_bs1"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BS1_ABSOLUTE = Info(
            "threshold_bs1_absolute",
            config_location=("allele_frequency_thresholds", "threshold_bs1_absolute"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BS1_SUPPORTING = Info(
            "threshold_bs1_supporting",
            config_location=("allele_frequency_thresholds", "threshold_bs1_supporting"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BS2 = Info(
            "threshold_bs2",
            config_location=("allele_frequency_thresholds", "threshold_bs2"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_BS2_SUPPORTING = Info(
            "threshold_bs2_supporting",
            config_location=("allele_frequency_thresholds", "threshold_bs2_supporting"),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_DIFF_LEN_PROT_PERCENT = Info(
            "threshold_diff_len_prot_percent",
            config_location=(
                "functional_thresholds",
                "threshold_diff_len_prot_percent",
            ),
            group=Classification_Info_Groups.THRESHOLD_SINGLE,
        )
        self.THRESHOLD_NMD = Info(
            "threshold_nmd",
            config_location=("disease_relevant_transcripts", "nmd_threshold"),
            group=Classification_Info_Groups.DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD,
            optional=True,
        )
        self.POS_LAST_KNOWN_PATHO_PTC = Info(
            "pos_last_known_patho_ptc",
            config_location=(
                "disease_relevant_transcripts",
                "pos_last_known_patho_ptc",
            ),
            group=Classification_Info_Groups.DISEASE_RELEVANT_TRANSCRIPT_THRESHOLD,
        )
        self.VARIANT = Info("variant")
        self.TRANSCRIPT = Info("transcript")
        self.CLINVAR_PATH = Info(
            "clinvar_path",
            config_location=("annotation_files", "clinvar", "clinvar_snv"),
            group=Classification_Info_Groups.PATH,
        )
        self.SIMILARITY_SCORE_PATH = Info(
            "similarity_score_path",
            config_location=(
                "annotation_files",
                "similarity_score",
                "similarity_score_file",
            ),
            group=Classification_Info_Groups.PATH,
        )
        self.SIMILARITY_SOCRE_DIRECTION = Info(
            ## Direction in which the similarity score has to be in order for a ClinVar entry to be admissable for PM5
            ## "greater" and "less" always are interpreted as "greater or equal" or "less or equal"
            "similarity_score_direction",
            config_location=(
                "annotation_files",
                "similarity_score",
                "similarity_score_direction",
            ),
            group=Classification_Info_Groups.CONFIG_ENTRY_STR,
        )
        self.UNIPROT_REP_REGION_PATH = Info(
            "uniprot_rep_region_path",
            config_location=("annotation_files", "uniprot", "rep"),
            group=Classification_Info_Groups.PATH,
        )
        self.CRITICAL_REGION_PATH = Info(
            "critical_region_path",
            config_location=("annotation_files", "critical_regions", "critical_region"),
            group=Classification_Info_Groups.PATH,
            optional=True,
        )
        self.DISEASE_IRRELEVANT_EXONS_PATH = Info(
            "disease_irrelevant_exons_path",
            config_location=(
                "annotation_files",
                "critical_regions",
                "disease_irrelevant_exons",
            ),
            group=Classification_Info_Groups.PATH,
            optional=True,
        )
        self.SPLICE_RESULT = Info("splice_result", optional=True)
        self.SPLICE_SITE_TABLE_PATH = Info(
            "splice_site_table_path",
            config_location=("annotation_files", "splice_site_table", "file"),
            group=Classification_Info_Groups.PATH,
        )
        self.SPLICE_RESULT_PM5 = Info("splice_result", optional=True)
        self.SPLICE_SITE_TABLE_PM5_PATH = Info(
            "splice_site_table_path",
            config_location=("annotation_files", "splice_site_table_pm5", "file"),
            group=Classification_Info_Groups.PATH,
        )
        self.EXON_PM5_PATH = Info(
            "exon_pm5_path",
            config_location=("annotation_files", "exon_pm5", "file"),
            group=Classification_Info_Groups.PATH,
        )
        self.PM5_RESULTS_PTC = Info("pm5_results_ptc", optional=True)
        self.FUNCTIONAL_ASSAY = Info("functional_assay")
        self.SPLICING_ASSAY = Info("splicing_assay")
