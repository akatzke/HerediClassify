#!/usr/bin/env python3

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import TypeVar, Optional, Generic

ValueType = TypeVar("ValueType")


@dataclass(frozen=False)
class classification(Generic[ValueType]):
    name: str
    config_location: Optional[tuple[str, ...]] = field(default=None)
    compute_function: Optional[Callable[[], ValueType]] = field(
        init=False, default=None
    )
    value: Optional[ValueType] = field(init=False, default=None)


class classification_information(Enum):
    ANNOTATED_TRANSCRIPT_LIST = classification("annotated_transcript_list")
    ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec = classification(
        "annotated_transcript_list_acmg"
    )
    VARIANT_CLINVAR = classification("variant_clinvar")
    VARIANT_HOTSPOT = classification("variant_hotspot")
    VARIANT_GNOMAD = classification("variant_gnomad")
    VARIANT_FLOSSIES = classification("variant_flossies")
    VARIANT_PREDICTION = classification("variant_prediction")
    THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN = classification(
        "prediction_pathogenicity",
        config_location=(
            "prediction_tool_thresholds",
            "pathogenicity_prediction",
            "benign",
        ),
    )
    THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC = classification(
        "prediction_pathogenicity",
        config_location=(
            "prediction_tool_thresholds",
            "pathogenicity_prediction",
            "pathogenic",
        ),
    )
    THRESHOLD_SPLICING_PREDICTION_BENIGN = classification(
        "prediction_splicing",
        config_location=("prediction_tool_thresholds", "splicing_prediction", "benign"),
    )
    THRESHOLD_SPLICING_PREDICTION_PATHOGENIC = classification(
        "prediction_splicing",
        config_location=(
            "prediction_tool_thresholds",
            "splicing_prediction",
            "pathogenic",
        ),
    )
    THRESHOLD_PM2 = classification(
        "threshold_pm2",
        config_location=("allele_frequency_thresholds", "threshold_pm2"),
    )
    THRESHOLD_BA1 = classification(
        "threshold_ba1",
        config_location=("allele_frequency_thresholds", "threshold_ba1"),
    )
    THRESHOLD_BS1 = classification(
        "threshold_bs1",
        config_location=("allele_frequency_thresholds", "threshold_bs1"),
    )
    THRESHOLD_BS2 = classification(
        "threshold_bs2",
        config_location=("allele_frequency_thresholds", "threshold_bs2"),
    )
    THRESHOLD_DIFF_LEN_PROT_PERCENT = classification(
        "threshold_diff_len_prot_percent",
        config_location=("functional_thresholds", "threshold_diff_len_prot_percent"),
    )
    VARIANT = classification("variant")
    TRANSCRIPT = classification("transcript")
    THRESHOLD_NMD = classification(
        "threshold_nmd",
        config_location=("functional_thresholds", "nmd_threshold"),
    )
    CLINVAR_PATH = classification(
        "clinvar_path", config_location=("annotation_files", "clinvar", "clinvar_snv")
    )
    UNIPROT_REP_REGION_PATH = classification(
        "uniprot_rep_region_path",
        config_location=("annotation_files", "uniprot", "rep"),
    )


class classification_information_groups(Enum):
    THRESHOLDS_SINGLE = {
        classification_information.THRESHOLD_PM2,
        classification_information.THRESHOLD_BA1,
        classification_information.THRESHOLD_BS1,
        classification_information.THRESHOLD_BS2,
        classification_information.THRESHOLD_DIFF_LEN_PROT_PERCENT,
        classification_information.THRESHOLD_NMD,
    }
    THRESHOLD_PREDICTION = {
        classification_information.THRESHOLD_PATHOGENICITY_PREDICTION_BENIGN,
        classification_information.THRESHOLD_PATHOGENICITY_PREDICTION_PATHOGENIC,
        classification_information.THRESHOLD_SPLICING_PREDICTION_BENIGN,
        classification_information.THRESHOLD_SPLICING_PREDICTION_PATHOGENIC,
    }
    PATH = {
        classification_information.CLINVAR_PATH,
        classification_information.UNIPROT_REP_REGION_PATH,
    }
