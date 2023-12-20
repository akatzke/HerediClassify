#!/usr/bin/env python3


class Pyensembl_no_coding_sequence(Exception):
    "Raised when pyensembl fails to produce a coding sequence"
    pass


class Pyensembl_transcript_not_found(Exception):
    "Raised when transcript id is not found in pyensembl"
    pass
