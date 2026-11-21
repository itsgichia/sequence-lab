"""Pairwise alignment: Needleman-Wunsch (global) and Smith-Waterman (local) with linear gaps."""

from __future__ import annotations

import numpy as np

from .sequence import BioError, clean, kind_of

MAX_CELLS = 4_000_000  # len(a) * len(b)


def align(a_raw: str, b_raw: str, mode: str = "global", match: int = 2, mismatch: int = -1, gap: int = -2) -> dict:
    """Needleman-Wunsch (global) or Smith-Waterman (local) with linear gaps."""
    a, b = clean(a_raw), clean(b_raw)
    kind_a, kind_b = kind_of(a), kind_of(b)
    if kind_a != kind_b:
        raise BioError("Both sequences must be the same type (both DNA or both protein).")
    if len(a) * len(b) > MAX_CELLS:
        raise BioError("These sequences are too long to align here (product of lengths over 4 million).")
    local = mode == "local"
    n, m = len(a), len(b)
    A = np.frombuffer(a.encode(), dtype=np.uint8)
    B = np.frombuffer(b.encode(), dtype=np.uint8)
    sub = np.where(A[:, None] == B[None, :], match, mismatch).astype(np.int32)

    H = np.zeros((n + 1, m + 1), dtype=np.int32)
    cols = np.arange(m + 1, dtype=np.int32)
    if not local:
        H[:, 0] = gap * np.arange(n + 1)
        H[0, :] = gap * cols
    # Each row: best[j] = max(diag, up) is known from the previous row; the "left" move
    # (a gap running along the row) is max_k<=j best[k] + gap*(j-k), which is a cumulative
    # max once gap*j is factored out. That keeps the whole row vectorised.
    for i in range(1, n + 1):
        best = np.empty(m + 1, dtype=np.int32)
        best[0] = H[i, 0]
        best[1:] = np.maximum(H[i - 1, :-1] + sub[i - 1], H[i - 1, 1:] + gap)
        if local:
            np.maximum(best, 0, out=best)
        H[i] = np.maximum.accumulate(best - gap * cols) + gap * cols

    if local:
        i, j = np.unravel_index(int(H.argmax()), H.shape)
        score = int(H[i, j])
    else:
        i, j, score = n, m, int(H[n, m])
    end_a, end_b = int(i), int(j)
    out_a, out_b = [], []
    while i > 0 or j > 0:
        h = H[i, j]
        if local and h == 0:
            break
        if i > 0 and j > 0 and h == H[i - 1, j - 1] + sub[i - 1, j - 1]:
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and h == H[i - 1, j] + gap:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    al_a, al_b = "".join(reversed(out_a)), "".join(reversed(out_b))
    matches = sum(1 for x, y in zip(al_a, al_b) if x == y)
    cols = len(al_a)
    midline = "".join("|" if x == y else " " if "-" in (x, y) else "." for x, y in zip(al_a, al_b))
    return {
        "mode": mode,
        "kind": kind_a,
        "score": score,
        "length": cols,
        "identity_percent": round(100 * matches / max(1, cols), 1),
        "gaps": al_a.count("-") + al_b.count("-"),
        "span_a": [int(i) + 1, end_a],
        "span_b": [int(j) + 1, end_b],
        "a": al_a,
        "b": al_b,
        "midline": midline,
    }
