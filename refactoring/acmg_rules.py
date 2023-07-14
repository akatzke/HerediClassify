#!/usr/bin/env python3
from mmap import PAGESIZE
from variant import *
from rule_utils import (
    Prediction_result,
    aggregate_prediction_results,
    assess_prediction_tool,
    summarise_results_per_transcript,
    assess_one_threshold,
    assess_two_thresholds,
)


def assess_pvs1(variant: Variant):
    result_start_loss = assess_pvs1_start_loss(variant)
    result_frameshift = assess_pvs1_frameshift(variant)
    result_splice = assess_pvs1_splice(variant)
    result = summarise_results_per_transcript(
        [result_start_loss, result_frameshift, result_splice]
    )
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
                results.append(
                    RuleResult("PVS1_start_loss", True, "supporting", comment)
                )
    result = summarise_results_per_transcript(results)
    return result


def assess_pvs1_splice(variant: Variant):
    results = []
    for transcript in variant.transcript_info:
        if transcript.is_NMD:
            if transcript.transcript_disease_relevant:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", True, "very_strong", comment))
            else:
                comment = "Something"
                results.append(RuleResult("PVS1_splice", False, "very_strong", comment))
        elif (
            transcript.exon_skipping
            and not transcript.is_NMD
            and not transcript.is_reading_frame_preserved
        ):
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
            if transcript.transcript_disease_relevant:
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


def assess_ps1(variant: Variant):
    if variant.clinvar.same_amino_acid_change_pathogenic:
        comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {variant.clinvar.same_amino_acid_change_pathogenic_list}."
        result = RuleResult("PS1", True, "strong", comment)
    else:
        comment = "No matches found for variant."
        result = RuleResult("PS1", False, "strong", comment)
    return result


def assess_pm1(variant: Variant):
    if variant.affected_region.critical_region:
        comment = f"Variant in mutational hotspot. {variant.affected_region.critical_region_type}"
        result = RuleResult("PM1", True, "moderate", comment)
    else:
        comment = "Variant not in mutational hotspot"
        result = RuleResult("PM1", False, "moderate", comment)
    return result


def assess_pm2(variant: Variant):
    if variant.gnomad.frequency > variant.thresholds["PM2"]:
        comment = (
            f"Variant occures with {variant.gnomad.frequency} in {variant.gnomad.name}"
        )
        result = RuleResult("PM2", False, "moderate", comment)
    else:
        comment = (
            f"Variant occurs with {variant.gnomad.frequency} in {variant.gnomad.name}"
        )
        result = RuleResult("PM2", True, "moderate", comment)
    return result


def assess_pm4(variant: Variant):
    results = []
    for transcript in variant.transcript_info:
        if (
            transcript.diff_len_protein_percent
            > Variant.thresholds["diff_len_protein_percent"]
            and transcript.transcript_disease_relevant
            and not transcript.len_change_in_repetitive_region
        ):
            comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region does not overlap repetitive region."
            result = RuleResult("PM4", True, "moderate", comment)
        elif (
            transcript.diff_len_protein_percent > 0.1
            and transcript.transcript_disease_relevant
            and transcript.len_change_in_repetitive_region
        ):
            comment = f"Length of disease relevant transcript {transcript.transcript_id} is reduced by {transcript.diff_len_protein_percent}. Deleted region overlaps repetitive region."
            result = RuleResult("PM4", False, "moderate", comment)
        else:
            comment = f"Length of transcript {transcript.transcript_id} altered by {transcript.diff_len_protein_percent}"
            result = RuleResult("PM4", False, "moderate", comment)
        results.append(result)
    final_result = summarise_results_per_transcript(results)
    return final_result


def assess_pm5(variant: Variant):
    if variant.clinvar.different_amino_acid_change_in_same_position_pathogenic:
        comment = f"The following ClinVar entries show the same amino acid change as pathogenic: {variant.clinvar.different_amino_acid_change_in_same_position_pathogenic_list}."
        result = RuleResult("PM5", True, "moderate", comment)
    else:
        comment = "No matches found for variant."
        result = RuleResult("PM5", False, "moderate", comment)
    return result


def assess_pp3(variant: Variant):
    splicing_prediction = []
    pathogenicity_prediction = []
    for entry in variant.prediction_tools:
        if entry.type == "pathogenicity":
            result = assess_prediction_tool(entry, variant.thresholds[entry.name])
            pathogenicity_prediction.append(result)
        elif entry.type == "splcing":
            result = assess_prediction_tool(entry, variant.thresholds[entry.name])
            splicing_prediction.append(result)
    splicing_result = aggregate_prediction_results(splicing_prediction)
    pathogenicity_result = aggregate_prediction_results(pathogenicity_prediction)
    if (
        splicing_result == Prediction_result.PATHOGENIC
        and pathogenicity_prediction == Prediction_result.PATHOGENIC
    ):
        comment = "Varinat is predicted to affect splicing and to be pathogenic."
        result = RuleResult("PP3", True, "supporting", comment)
    elif splicing_result == Prediction_result.PATHOGENIC:
        comment = "Variant is predicted to affect splicing."
        result = RuleResult("PP3", True, "supporting", comment)
    elif pathogenicity_result == Prediction_result.PATHOGENIC:
        comment = "Variant is predicted to be pathogenic."
        result = RuleResult("PP3", True, "supporting", comment)
    else:
        comment = (
            "Varinat is predicted to not affect splicing and to not be pathogenic."
        )
        result = RuleResult("PP3", False, "supporting", comment)
    return result


def assess_ba1(variant: Variant):
    return variant


def assess_bs1(variant: Variant):
    return variant


def assess_bp3(variant: Variant):
    return variant


def assess_bp4(variant: Variant):
    return variant


def assess_bp7(variant: Variant):
    return variant
