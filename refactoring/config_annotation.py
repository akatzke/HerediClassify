#!/usr/bin/env python3

from typing import Protocol
from enum import Enum

from refactoring.variant_annotate import Variant_annotated


class classification_information(Enum):
    ANNOTATED_TRANSCRIPT_LIST = "annotated_transcript_list"
    VARIANT_CLINVAR = "variant_clinvar"
    VARIANT_HOTSPOT = "variant_hotspot"
    VARIANT_ANNOT = "variant_annot"
    THRESHOLD_PATHOGENICITY_PREDICTION = "threshold_pathogenicity_prediction"
    THRESHOLD_SPLICING_PREDICTION = "threshold_splicing_prediction"
    THRESHOLD_PM2 = "threshold_pm2"
    THRESHOLD_BA1 = "threshold_ba1"
    THRESHOLD_BS1 = "threshold_ps1"
    THRESHOLD_BS2 = "threshold_bs2"


class Annotation(Protocol):
    """
    Function that annotates a Variant_annotated object
    """

    def __call__(self, variant: Variant_annotated) -> Variant_annotated:
        """
        Return a further annotated Variant_annotated object
        """
        ...
