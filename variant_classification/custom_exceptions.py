#!/usr/bin/env python3


class Pyensembl_no_coding_sequence(Exception):
    "Raised when pyensembl fails to produce a coding sequence"
    pass


class Not_disease_relevant_transcript(Exception):
    "Raised when threshold for disease relevant transcript is accessed but transcript is not disease relevant."
    pass
