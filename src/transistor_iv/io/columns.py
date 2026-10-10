"""Column-header matching shared by the Clarius and CSV readers.

Headers are compared after lower-casing and removing every non-alphanumeric
character, so ``"Drain I (A)"`` becomes ``"driania"`` and ``"I_d"`` matches ``"Id"``. A column is found by,
in this order: an exact match of the user-given header, an exact match of
one of the default candidates, a prefix match of either.
"""

from __future__ import annotations

import re

DEFAULT_CANDIDATES: dict[str, tuple[str, ...]] = {
    "i_d": ("DrainI", "DrainI_A", "Id", "Drain I", "Drain I (A)", "Id(A)"),
    "v_g": ("GateV", "Vg", "Gate V", "Gate V (V)", "Vg(V)"),
    "v_d": ("DrainV", "Vd", "Drain V", "Drain V (V)"),
    "i_g": ("GateI", "Ig", "Gate I", "Gate I (A)", "GateI_A"),
}


def normalize(header) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(header).strip().lower())


def find_column(headers: list[str], quantity: str, custom: str | None = None) -> int | None:
    """Index of the column holding ``quantity`` ('i_d', 'v_g', 'v_d', 'i_g'), or None."""
    norm = [normalize(h) for h in headers]
    pool = ([custom] if custom else []) + list(DEFAULT_CANDIDATES[quantity])
    for cand in pool:
        c = normalize(cand)
        if c in norm:
            return norm.index(c)
    for cand in pool:
        c = normalize(cand)
        for idx, h in enumerate(norm):
            if h.startswith(c):
                return idx
    return None


def split_numbered(header) -> tuple[str, int | None]:
    """``'DrainI(3)'`` -> ``('draini', 3)``; ``'DrainI'`` -> ``('draini', None)``."""
    m = re.match(r"^(.+?)\((\d+)\)$", str(header).strip())
    if m:
        return normalize(m.group(1)), int(m.group(2))
    return normalize(header), None


def find_numbered_columns(headers: list[str], quantity: str, custom: str | None = None) -> dict[int, int]:
    """``{measurement number: column index}`` for numbered headers ``Name(1)``, ``Name(2)``, ...

    Clarius stores the curves of one output-curve test (one per gate step)
    in numbered columns. Matching order as in :func:`find_column`.
    """
    parsed = [split_numbered(h) for h in headers]
    pool = ([custom] if custom else []) + list(DEFAULT_CANDIDATES[quantity])
    for exact in (True, False):
        for cand in pool:
            c = normalize(cand)
            found = {num: idx for idx, (base, num) in enumerate(parsed)
                     if num is not None and (base == c if exact else base.startswith(c))}
            if found:
                return found
    return {}
