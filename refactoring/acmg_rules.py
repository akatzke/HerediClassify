#!/usr/bin/env python3
from variant import *
from rule_utils import summarise_results_per_transcript

def assess_pvs1(variant: Variant):
    result_start_loss = assess_pvs1_start_loss(variant)
    result_frameshift = assess_pvs1_frameshift(variant)
    result_splice = assess_pvs1_splice(variant)
    result = summarise_results_per_transcript([result_start_loss, result_frameshift, result_splice])
    return result

def assess_pvs1_start_loss(variant: Variant):
    results = []
    for transcript in variant.transcript_info:
        if transcript.exists_alternative_start_codon:
            comment = "Something"
            results.append(RuleResult("PVS1_start_loss", False, "very_strong", comment))
        else:
            if transcript.pathogenic_variant_between_start_and_stop:
                comment = "Something"
                results.append(RuleResult("PVS1_start_loss", True, "moderate", comment))
            else:
                comment = "Something"
                results.append(RuleResult("PVS1_start_loss", True, "supporting", comment))
    result = summarise_results_per_transcript(results)
    return result


def assess_pvs1_splice(variant: Variant):
    results = []
    for transcript in variant.transcript_info:
        if transcript.is_NMD:
            if transcript.NMD_affected_transcript_disease_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "very_strong", comment))
            else:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", False, "very_strong", comment))
        elif transcript.exon_skipping and not transcript.is_NMD and not transcript.is_reading_frame_preserved:
            if transcript.truncated_exon_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "strong", comment))
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "strong", comment))
                else:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "moderate", comment))
        elif transcript.exon_skipping or transcript.is_reading_frame_preserved:
            if transcript.truncated_exon_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "strong", comment))
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "strong", comment))
                else:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "moderate", comment))

        else:
            comment = "Something"
            results.append(RuleResult("PVS1_splice", False, "very_strong", comment))
    result = summarise_results_per_transcript(results)
    return result



def assess_pvs1_frameshift(variant: Variant):
    results = []
    for transcript in variant.transcript_info:
        if transcript.is_NMD:
            if transcript.NMD_affected_transcript_disease_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "very_strong", comment))
            else:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", False, "very_strong", comment))
        else:
            if transcript.truncated_exon_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "strong", comment))
            else:
                if transcript.diff_len_protein_percent > 0.1:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "strong", comment))
                else:
                    comment = "Something"
                    results.append(RuleResult("PVS1_splice", True, "moderate", comment))
    result = summarise_results_per_transcript(results)
    return result
