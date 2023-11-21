#!/usr/bin/env python3


def assess_reading_frame_preservation(diff_len: int) -> bool:
    """
    Check if reading frame is preserved
    """

    if diff_len % 3 == 0:
        return True
    else:
        # reading frame is disrupted
        return False
