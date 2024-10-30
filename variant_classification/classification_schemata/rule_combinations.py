#!/usr/bin/env python3

from classification_schemata.utils import generate_count_rule

# Define all combination of rules for a benign classification
benign_1 = generate_count_rule(min_benign_stand_alone=1)
benign_2 = generate_count_rule(min_benign_strong=2)
benign_brca1_1 = generate_count_rule(min_benign_strong=1, min_benign_moderate=2)
benign_brca1_2 = generate_count_rule(
    min_benign_strong=1, min_benign_moderate=1, min_benign_supporting=1
)
benign_brca1_3 = generate_count_rule(min_benign_strong=1, min_benign_supporting=3)


# Define all combination of rules for a likely benign classification
likely_benign_1 = generate_count_rule(min_benign_strong=1, min_benign_supporting=1)
likely_benign_2 = generate_count_rule(min_benign_supporting=2)
likely_benign_atm_1 = generate_count_rule(min_benign_strong=1)
likely_benign_brca1_1 = generate_count_rule(min_benign_strong=1, min_benign_moderate=1)
likely_benign_brca1_2 = generate_count_rule(
    min_benign_moderate=1, min_benign_supporting=1
)


# Define all combination of rules for a pathogenic classification
pathogenic_1 = generate_count_rule(min_pathogenic_very_strong=2)
pathogenic_2 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_strong=1
)
pathogenic_3 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_moderate=2
)
pathogenic_4 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_moderate=1, min_pathogenic_supporting=1
)
pathogenic_5 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_supporting=2
)
pathogenic_6 = generate_count_rule(min_pathogenic_strong=2)
pathogenic_7 = generate_count_rule(min_pathogenic_strong=1, min_pathogenic_moderate=3)
pathogenic_8 = generate_count_rule(
    min_pathogenic_strong=1, min_pathogenic_moderate=2, min_pathogenic_supporting=2
)
pathogenic_9 = generate_count_rule(
    min_pathogenic_strong=1, min_pathogenic_moderate=1, min_pathogenic_supporting=4
)
pathogenic_brca1_1 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_moderate=1
)
pathogenic_brca1_2 = generate_count_rule(min_pathogenic_strong=3)
pathogenic_brca1_3 = generate_count_rule(
    min_pathogenic_strong=2, min_pathogenic_moderate=1
)
pathogenic_brca1_4 = generate_count_rule(
    min_pathogenic_strong=2, min_pathogenic_supporting=2
)


# Define all combination of rules for a likely pathogenic classification
likely_pathogenic_1 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_moderate=1
)
# This is added in line with recommendations
likely_pathogenic_2 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_supporting=1
)
likely_pathogenic_3 = generate_count_rule(
    min_pathogenic_strong=1, min_pathogenic_moderate=1
)
likely_pathogenic_4 = generate_count_rule(
    min_pathogenic_strong=1, min_pathogenic_supporting=2
)
likely_pathogenic_5 = generate_count_rule(min_pathogenic_moderate=3)
likely_pathogenic_6 = generate_count_rule(
    min_pathogenic_moderate=2, min_pathogenic_supporting=2
)
likely_pathogenic_7 = generate_count_rule(
    min_pathogenic_moderate=1, min_pathogenic_supporting=4
)
likely_pathogenic_brca1_1 = generate_count_rule(
    min_pathogenic_strong=1, min_pathogenic_moderate=1
)
likely_pathogenic_brca1_2 = generate_count_rule(
    min_pathogenic_very_strong=1, min_pathogenic_supporting=2
)
