#!/usr/bin/env python3


class Pyensembl_no_coding_sequence(Exception):
    "Raised when pyensembl fails to produce a coding sequence"
    pass

class Pyensembl_transcript_not_found(Exception):
    "Raised when transcript id is not found in pyensembl"
    pass

class No_transcript_with_var_type_found(Exception):
    "Raise when no variant is found to match variant type."
    pass

class Not_disease_relevant_transcript(Exception):
    "Raised when threshold for disease relevant transcript is accessed but transcript is not disease relevant."
    pass
