#!/usr/bin/env python3

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional


@dataclass(frozen=False)
class classification:
    config_location: Optional[tuple[str, ...]] = field(default=None)
    compute_function: Optional[Callable[[], Any]] = field(init=False, default=None)
    value: Optional[Any] = field(init=False)


class classification_information(Enum):
    ANNOTATED_TRANSCRIPT_LIST = classification()
    ANNOTATED_TRANSCRIPT_LIST_ACMG_Spec = classification()
    VARIANT_CLINVAR = classification()
    VARIANT_ANNOT = classification()
    THRESHOLD_PATHOGENICITY_PREDICTION = classification(
        config_location=("prediction_tool_threshold", "pathogenicity_prediction"),
    )
    THRESHOLD_SPLICING_PREDICTION = classification(
        config_location=("prediction_tool_threshold", "splicing_prediction"),
    )
    THRESHOLD_PM2 = classification(
        config_location=("allele_frequency_threshold", "threshold_pm2"),
    )
    THRESHOLD_BA1 = classification(
        config_location=("allele_frequency_threshold", "threshold_ba1"),
    )
    THRESHOLD_BS1 = classification(
        config_location=("allele_frequency_threshold", "threshold_bs1"),
    )
    THRESHOLD_BS2 = classification(
        config_location=("allele_frequency_threshold", "threshold_bs2"),
    )
    THRESHOLD_DIFF_LEN_PROT_PERCENT = classification(
        config_location=("allele_frequency_threshold", "threshold_bs2"),
    )
    VARIANT = classification()
    TRANSCRIPT = classification()
    NMD_THRESHOLD = classification(
        config_location=("functional_thresholds", "nmd_threshold"),
    )
    CLINVAR_PATH = classification()
    UNIPROT_REP_REGION_PATH = classification()


def compute_information(info: classification_information, config: dict):
    try:
        helper_function_access_dict(config, info.value.config_location)
    except KeyError:
        info.value.compute_function()
