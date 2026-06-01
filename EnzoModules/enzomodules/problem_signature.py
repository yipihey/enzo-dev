"""Compact, deterministic "signature" of an initialized Enzo problem.

Storing full 2D/3D baryon fields as fixtures would be huge, so instead we
summarize each field by statistics computed over the FULL field (including
ghost zones): count, min, max, mean, sum.  Combined with structural metadata
(per-grid rank/dims/num_particles, problem_type, num_grids, field names/types)
this is enough to catch regressions in the initial-condition generators.

The signature is a plain ``dict`` of JSON-serializable values.  Floats are
stored via :func:`_round`, which keeps ~12 significant digits so the values
round-trip deterministically through JSON.

Used by both ``tools/capture_problems.py`` (capture) and
``tests/test_problem_fixtures.py`` (regression replay).
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

# Schema version, bumped if the signature layout ever changes.
SCHEMA = 1

# Significant digits retained for stored statistics.
_SIGFIG = 12


def _round(x: float) -> float:
    """Round a float to ~12 significant digits so JSON round-trips stably."""
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    # Decimal exponent of the leading digit.
    exp = math.floor(math.log10(abs(x)))
    ndigits = _SIGFIG - 1 - exp
    return round(x, ndigits)


def _field_stats(values: List[float]) -> Dict[str, float]:
    """Statistics over the full flat field (pure Python; no numpy)."""
    n = len(values)
    if n == 0:
        return {"count": 0, "min": 0.0, "max": 0.0, "mean": 0.0, "sum": 0.0}
    vmin = vmax = values[0]
    total = 0.0
    for v in values:
        if v < vmin:
            vmin = v
        if v > vmax:
            vmax = v
        total += v
    return {
        "count": n,
        "min": _round(vmin),
        "max": _round(vmax),
        "mean": _round(total / n),
        "sum": _round(total),
    }


def signature(problem) -> Dict:
    """Compute the compact signature dict for an initialized ``Problem``."""
    grids = []
    for gi in range(problem.num_grids):
        g = problem.grid(gi)
        types = g.field_types
        names = g.field_names
        fields = []
        for fi, (t, name) in enumerate(zip(types, names)):
            stats = _field_stats(g.field(fi))
            fields.append({
                "name": name,
                "type": int(t),
                "dims": list(g.dims),
                "stats": stats,
            })
        grids.append({
            "index": gi,
            "rank": int(g.rank),
            "dims": list(g.dims),
            "num_particles": int(g.num_particles),
            "fields": fields,
        })
    return {
        "schema": SCHEMA,
        "problem_type": int(problem.problem_type),
        "num_grids": int(problem.num_grids),
        "grids": grids,
    }


def _close(a: float, b: float, rtol: float, atol: float) -> bool:
    a = float(a)
    b = float(b)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))


def compare_signatures(a: Dict, b: Dict, rtol: float = 1e-10,
                       atol: float = 1e-300) -> Tuple[bool, str]:
    """Compare two signatures.  Returns ``(ok, message)``.

    Structural fields (problem_type, num_grids, dims, field types/names,
    counts) must match exactly; statistics (min/max/mean/sum) must match
    within the given floating tolerance.
    """
    if a.get("problem_type") != b.get("problem_type"):
        return False, (f"problem_type differs: "
                       f"{a.get('problem_type')} != {b.get('problem_type')}")
    if a.get("num_grids") != b.get("num_grids"):
        return False, (f"num_grids differs: "
                       f"{a.get('num_grids')} != {b.get('num_grids')}")

    ga, gb = a.get("grids", []), b.get("grids", [])
    if len(ga) != len(gb):
        return False, f"grid count differs: {len(ga)} != {len(gb)}"

    for gi, (g1, g2) in enumerate(zip(ga, gb)):
        if g1.get("rank") != g2.get("rank"):
            return False, f"grid {gi}: rank {g1.get('rank')} != {g2.get('rank')}"
        if g1.get("dims") != g2.get("dims"):
            return False, f"grid {gi}: dims {g1.get('dims')} != {g2.get('dims')}"
        if g1.get("num_particles") != g2.get("num_particles"):
            return False, (f"grid {gi}: num_particles "
                           f"{g1.get('num_particles')} != {g2.get('num_particles')}")

        f1, f2 = g1.get("fields", []), g2.get("fields", [])
        if len(f1) != len(f2):
            return False, f"grid {gi}: field count {len(f1)} != {len(f2)}"

        for fa, fb in zip(f1, f2):
            tag = f"grid {gi} field '{fa.get('name')}'"
            if fa.get("type") != fb.get("type"):
                return False, (f"{tag}: type {fa.get('type')} != {fb.get('type')}")
            if fa.get("name") != fb.get("name"):
                return False, (f"{tag}: name {fa.get('name')} != {fb.get('name')}")
            if fa.get("dims") != fb.get("dims"):
                return False, (f"{tag}: dims {fa.get('dims')} != {fb.get('dims')}")
            s1, s2 = fa.get("stats", {}), fb.get("stats", {})
            if s1.get("count") != s2.get("count"):
                return False, (f"{tag}: count {s1.get('count')} != {s2.get('count')}")
            for key in ("min", "max", "mean", "sum"):
                if not _close(s1.get(key), s2.get(key), rtol, atol):
                    return False, (f"{tag}: {key} {s1.get(key)} != {s2.get(key)} "
                                   f"(rtol={rtol}, atol={atol})")
    return True, "signatures match"
