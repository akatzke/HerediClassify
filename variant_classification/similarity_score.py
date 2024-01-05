#!/usr/bin/env python3

import pathlib
import logging
import pandas as pd
from typing import Optional

from Bio.Seq import IUPACData
import hgvs.parser

hgvs_parser = hgvs.parser.Parser()

logger = logging.getLogger("GenOtoScope_Classify.clinvar.missense")


def get_similarity_score_clinvar(
    data: pd.DataFrame, path_similarity_score: pathlib.Path
) -> pd.DataFrame:
    """
    For a set of amino acid exchanges produce the similarity score
    The dataframe expected is created from ClinVar entries
    """
    similarity_score = pd.read_csv(path_similarity_score, sep="\t")
    data_similarity_score = pd.DataFrame = data.apply(
        get_similarity_score_df, similarity_score=similarity_score, axis=1
    )
    return data_similarity_score


def get_similarity_score_df(
    clinvar_entry: pd.Series, similarity_score: pd.DataFrame
) -> pd.Series:
    """
    Annotate dataframe
    """
    ref_aa = clinvar_entry["prot_ref"]
    alt_aa = clinvar_entry["prot_alt"]
    score = look_up_similarity_score(ref_aa, alt_aa, similarity_score)
    clinvar_entry["similarity_score"] = score
    return clinvar_entry


def get_similarity_score(
    var_codon_info: dict, path_similarity_score: pathlib.Path
) -> Optional[int]:
    """
    Get the Grantham score for the current variant
    """
    similarity_score = pd.read_csv(path_similarity_score, sep="\t")
    ref_aa = var_codon_info["prot_ref"]
    alt_aa = var_codon_info["prot_alt"]
    score = look_up_similarity_score(ref_aa, alt_aa, similarity_score)
    return score


def look_up_similarity_score(
    ref_aa: str, alt_aa: str, similarity_score: pd.DataFrame
) -> Optional[int]:
    """
    Check similarity file for score
    """
    if len(ref_aa) != 1:
        try:
            ref_aa = IUPACData.protein_letters_3to1[ref_aa]
        except KeyError:
            raise ValueError(
                f"For the reference amino acid {ref_aa} no matching one letter amino acid code could be found."
            )
    if len(alt_aa) != 1:
        try:
            alt_aa = IUPACData.protein_letters_3to1[alt_aa]
        except KeyError:
            raise ValueError(
                f"For the alternative amino acid {alt_aa} no matching one letter amino acid code could be found."
            )
    try:
        score = similarity_score[
            (similarity_score.ref_aa == ref_aa) & (similarity_score.alt_aa == alt_aa)
        ].score.values[0]
        return score
    except Exception:
        return None
