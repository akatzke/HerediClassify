#!/usr/bin/env python3

from enum import Enum


class VARTYPE(Enum):
    STOP_GAINED = "stop_gained"
    STOP_LOST = "stop_lost"
    FRAMESHIFT_VARIANT = "frameshift_variant"
    INFRAME_DELETION = "inframe_deletion"
    INFRAME_INSERTION = "inframe_insertion"
    TRANSCRIPT_ABLATION = "transcript_ablation"
    TRANSCRIPT_AMPLIFICATION = "transcript_amplification"
    FEATURE_ELONGATION = "feature_elongation"
    FEATURE_TRUNCATION = "feature truncation"
    MISSENSE_VARIANT = "missense_variant"
    PROTEIN_ALTERING_VARIANT = "protein_altering_variant"
    INCOMPLETE_TERMINAL_CODON_VARIANT = "incomplete_terminal_codon_variant"
    START_RETAINED_VARIANT = "start_retained_variant"
    STOP_RETAINED_VARIANT = "stop_retained_variant"
    SYNONYMOUS_VARIANT = "synonymous_variant"
    CODING_SEQUENCE_VARIANT = "coding_sequence_variant"
    MATURE_MIRNA_VARIANT = "mature_mirna_variant"
    FIVE_PRIME_UTR_VARIANT = "5_prime_utr_variant"
    THREE_PRIME_UTR_VARIANT = "3_prime_utr_variant"
    NON_CODING_TRANSCRIPT_EXON_VARIANT = "non_coding_transcript_exon_variant"
    INTRON_VARIANT = "intron_variant"
    NMD_TRANSCRIPT_VARIANT = "nmd_transcript_variant"
    NON_CODING_TRANSCRIPT_VARIANT = "non_coding_transcript_variant"
    CODING_TRANSCRIPT_VARIANT = "coding_transcript_variant"
    UPSTREAM_GENE_VARIANT = "upstream_gene_variant"
    DOWNSTREAM_GENE_VARIANT = "downstream_gene_variant"
    TFBS_ABLATION = "tfbs_ablation"
    TFBS_AMPLIFICATION = "tfbs_amplification"
    TP_BINDING_SITE_VARIANT = "tf_binding_site_variant"
    REGULATORY_REGION_ABLATION = "regulatory_region_ablation"
    REGULATORY_REGION_AMPLIFICATION = "regulatory_region_amplification"
    REGULATORY_REGION_VARIANT = "regulatory_region_variant"
    INTERGENIC_VARIANT = "intergenic_variant"
    SEQUENCE_VARIANT = "sequence_variant"
    START_LOST = "start_lost"
    SPLICE_DONOR_5TH_BASE_VARIANT = "splice_donor_5th_base_variant"
    SPLICE_REGION_VARIANT = "splice_region_variant"
    SPLICE_DONOR_REGION_VARIANT = "splice_donor_region_variant"
    SPLICE_POLYPYRIMIDINE_TRACT_VARIANT = "splice_polypyrimidine_tract_variant"
    SPLICE_DONOR_VARIANT = "splice_donor_variant"
    SPLICE_DONOR = "splice_donor"
    SPLICE_ACCEPTOR_VARIANT = "splice_acceptor_variant"
    SPLICE_ACCEPTOR = "splice_acceptor"


class VARTYPE_GROUPS(Enum):
    EXONIC = {
        VARTYPE.STOP_GAINED,
        VARTYPE.STOP_LOST,
        VARTYPE.FRAMESHIFT_VARIANT,
    }
    START_LOST = {VARTYPE.START_LOST}
    INTRONIC = {
        VARTYPE.SPLICE_ACCEPTOR,
        VARTYPE.SPLICE_ACCEPTOR_VARIANT,
        VARTYPE.SPLICE_DONOR,
        VARTYPE.SPLICE_DONOR_VARIANT,
    }
    MISSENSE = {VARTYPE.MISSENSE_VARIANT}
    PREDICTION_PROTEIN = {
        VARTYPE.MISSENSE_VARIANT,
        VARTYPE.INFRAME_DELETION,
        VARTYPE.INFRAME_INSERTION,
    }
    PREDICTION_SPLICING = {
        VARTYPE.SYNONYMOUS_VARIANT,
        VARTYPE.MISSENSE_VARIANT,
        VARTYPE.INFRAME_DELETION,
        VARTYPE.INFRAME_INSERTION,
        VARTYPE.INTRON_VARIANT,
    }
    PREDICTION_NO_PROTEIN = {
        VARTYPE.START_LOST,
        VARTYPE.STOP_GAINED,
        VARTYPE.SPLICE_ACCEPTOR,
        VARTYPE.SPLICE_ACCEPTOR_VARIANT,
        VARTYPE.SPLICE_DONOR,
        VARTYPE.SPLICE_DONOR_VARIANT,
    }
    PREDICTION_BP7 = {VARTYPE.SYNONYMOUS_VARIANT, VARTYPE.INTRON_VARIANT}
