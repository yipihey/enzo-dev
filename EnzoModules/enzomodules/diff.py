"""Tolerance policy for comparing kernel outputs.

Kept small and dependency-free so any downstream rewrite (Python *or* Julia,
via its own reader of the shared fixtures) can reuse the same notion of
"equal".  The relative+absolute scheme mirrors Enzo's existing yt answer
testing, so the unit layer and the simulation layer agree on what equal means.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class Tolerance:
    """A comparison tolerance.

    ``rtol == atol == 0`` means bitwise equality (use for integer/index
    logic).  Use a small positive ``rtol`` for floating-point kernels that
    need not match a Fortran reference bit-for-bit (FMA / reduction order).
    """

    rtol: float = 0.0
    atol: float = 0.0


#: Bitwise-equality tolerance.
BITWISE = Tolerance(0.0, 0.0)


def isclose(a: float, b: float, tol: Tolerance) -> bool:
    """True when ``abs(a - b) <= atol + rtol*abs(b)``.

    NaN-aware: two NaNs compare equal, so identical reference outputs (which
    may legitimately contain NaN) round-trip as equal.
    """
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if tol.rtol == 0.0 and tol.atol == 0.0:
        return a == b
    return abs(a - b) <= tol.atol + tol.rtol * abs(b)


@dataclass(frozen=True)
class CompareResult:
    """Result of an element-wise array comparison."""

    ok: bool
    n: int
    nfail: int
    maxabs: float          # largest absolute difference
    maxrel: float          # largest relative difference (vs expected)
    worst: int             # 1-based index of the worst element (0 if none)

    def __bool__(self) -> bool:  # so `assert compare(...)` reads naturally
        return self.ok


def compare(actual: Sequence[float], expected: Sequence[float],
            tol: Tolerance) -> CompareResult:
    """Element-wise comparison of two equal-length numeric sequences."""
    a = [float(x) for x in actual]
    e = [float(x) for x in expected]
    if len(a) != len(e):
        raise ValueError(f"length mismatch: actual {len(a)} vs expected {len(e)}")

    nfail = 0
    maxabs = 0.0
    maxrel = 0.0
    worst = 0
    for i, (av, bv) in enumerate(zip(a, e), start=1):
        if not isclose(av, bv, tol):
            nfail += 1
        if math.isnan(av) or math.isnan(bv):
            da = 0.0 if (math.isnan(av) and math.isnan(bv)) else math.inf
        else:
            da = abs(av - bv)
        dr = (0.0 if av == 0.0 else math.inf) if bv == 0.0 else da / abs(bv)
        if da > maxabs:
            maxabs = da
        if dr > maxrel:
            maxrel = dr
            worst = i
    return CompareResult(nfail == 0, len(a), nfail, maxabs, maxrel, worst)
