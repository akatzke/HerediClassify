#!/usr/bin/env python3

from typing import Protocol
from refactoring.variant_annotate import Variant_annotated


class Annotation(Protocol):
    """
    Function that annotates a Variant_annotated object
    """

    def __call__(self, variant: Variant_annotated) -> Variant_annotated:
        """
        Return a further annotated Variant_annotated object
        """
        ...
