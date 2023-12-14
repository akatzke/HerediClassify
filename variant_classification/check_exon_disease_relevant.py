#!/usr/bin/env python3

import pathlib
import pandas as pd


def check_exon_disease_relevant(
    path_disease_irrelevant_exons: pathlib.Path, NMD_affected_exons: list[dict]
) -> bool:
    """
    Check if NMD affected exon has been listed as not relevant for disease
    Usually means there is a rescue transcript for the cases where the affected exon is skipped/truncated
    """
    disease_irrelevant_exons = pd.read_csv(path_disease_irrelevant_exons, sep="\t")
    disease_irrelevant_exon_ids = disease_irrelevant_exons["exon_name"]
    is_exon_disease_relevant = True
    for NMD_affected_exon in NMD_affected_exons:
        if NMD_affected_exon["exon_id"] in disease_irrelevant_exon_ids.values:
            is_exon_disease_relevant = False
        else:
            return True
    return is_exon_disease_relevant
