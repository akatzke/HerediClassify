#!/usr/bin/env python3

import pathlib
from test.example_variant_indel_GRCh38 import create_example_dup

from test.example_variant_splicing_GRCh38 import create_example_splice_acceptor_BRCA1
from variant_classification.transcript_annotated import (
    TranscriptInfo_exonic,
    TranscriptInfo_intronic,
    TranscriptInfo_start_loss,
)


def test_splicing():
    test_trans, test_var = create_example_splice_acceptor_BRCA1()
    path_clinvar = pathlib.Path("/home/katzkean/databases/Clinvar/clinvar_snv.vcf.gz")
    path_uniprot = pathlib.Path(
        "/home/katzkean/databases/Uniprot/repeats_hg38_uniprot.bed"
    )
    path_critical_region = pathlib.Path(
        "/home/katzkean/variant_classification/data/critical_region/VUS-Task-Force_critical_protein_domains.bed"
    )
    annot_trans = TranscriptInfo_intronic.annotate(
        test_var, path_clinvar, path_uniprot, path_critical_region, test_trans
    )
    assert (
        annot_trans.are_exons_skipped == True
        and annot_trans.is_NMD == True
        and annot_trans.is_truncated_region_disease_relevant == True
        and annot_trans.len_change_in_repetitive_region == False
        and annot_trans.is_reading_frame_preserved == False
    )


def test_indel():
    test_trans, test_var = create_example_dup()
    path_clinvar = pathlib.Path("/home/katzkean/databases/Clinvar/clinvar_snv.vcf.gz")
    path_uniprot = pathlib.Path(
        "/home/katzkean/databases/Uniprot/repeats_hg38_uniprot.bed"
    )
    path_critical_region = pathlib.Path(
        "/home/katzkean/variant_classification/data/critical_region/VUS-Task-Force_critical_protein_domains.bed"
    )
    annot_trans = TranscriptInfo_exonic.annotate(
        test_var, path_clinvar, path_uniprot, path_critical_region, test_trans
    )
