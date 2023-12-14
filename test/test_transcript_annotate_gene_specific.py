#!/usr/bin/env python3

import pathlib

from variant_classification.load_config import load_config
import test.paths as paths
from variant_classification.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)


def test_NMD_threshold_true():
    """
    Test NMD threshold
    See that it is true
    """
    pass


def test_NMD_threshold_false():
    """
    Test NMD thresold
    See that it is False
    """
    pass


def test_disease_relevant_transcript():
    """
    Test that filtering of disease relevant transcripts works
    """
    pass


def test_exon_not_disease_relevant():
    """
    Test that not disease relevant exon is recognised
    """
    pass


def test_frameshift():
    """
    Test that frameshift is correctly identified
    """
    pass


def test_ptc_position():
    """
    Test that the correct ptc position is being created
    """
    pass
