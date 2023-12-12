#!/usr/bin/env python3


def assess_reading_frame_preservation(diff_len: int) -> tuple[bool, int]:
    """
    Check if reading frame is preserved

    Returns
    -------
    bool:
        Is reading frame preserved
    int:
        Shift in reading frame (either 0, +1 or -1)
    """

    if diff_len % 3 == 0:
        return True, 0
    elif diff_len % 3 == 1:
        # Reading frame is shifted by 1
        return False, +1
    elif diff_len % 3 == 2:
        # Reading frame is shifted by -1
        return False, -1
    else:
        raise ValueError("Error whilst assessing reading frame shift.")
