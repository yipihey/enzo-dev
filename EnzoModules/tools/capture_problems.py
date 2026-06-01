#!/usr/bin/env python3
"""Capture golden signature fixtures for Enzo initial conditions.

For a fixed list of self-contained (analytic, no external data) parameter
files, drive Enzo's ``InitializeNew`` via ``enzomodules.problems.Problem``,
compute a compact statistical signature of the resulting grid hierarchy, and
write it as JSON to ``fixtures/Problems/<name>.json``.

The full baryon fields are NOT stored (they are huge); see
``enzomodules.problem_signature`` for what the signature contains.

Run from the EnzoModules directory:

    python3 tools/capture_problems.py
"""
import json
import os
import sys

# repo root = three levels up from this file:
#   <repo>/EnzoModules/tools/capture_problems.py
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
ENZOMODULES = os.path.abspath(os.path.join(_HERE, ".."))

sys.path.insert(0, ENZOMODULES)

from enzomodules import problems                       # noqa: E402
from enzomodules.problem_signature import signature    # noqa: E402

OUT_DIR = os.path.join(ENZOMODULES, "fixtures", "Problems")

# Self-contained problems known to initialize from analytic ICs only.
# Paths are relative to the repo root.
PARAM_FILES = [
    "run/Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo",
    "run/Hydro/Hydro-1D/Toro-2-ShockTube/Toro-2-ShockTube.enzo",
    "run/Hydro/Hydro-2D/SedovBlast/SedovBlast.enzo",
    "run/Hydro/Hydro-2D/Implosion/Implosion.enzo",
    "run/Hydro/Hydro-2D/KelvinHelmholtz/KelvinHelmholtz.enzo",
    "run/Hydro/Hydro-2D/NohProblem2D/NohProblem2D.enzo",
]


def _name(rel):
    return os.path.splitext(os.path.basename(rel))[0]


def main():
    if not problems.available():
        print("grid solver library not built; cannot capture.")
        return 1

    os.makedirs(OUT_DIR, exist_ok=True)
    written = 0
    for rel in PARAM_FILES:
        name = _name(rel)
        path = os.path.join(REPO, rel)
        if not os.path.exists(path):
            print(f"SKIP  {name}: parameter file missing ({rel})")
            continue
        try:
            with problems.Problem(path) as p:
                sig = signature(p)
        except Exception as exc:  # noqa: BLE001
            print(f"SKIP  {name}: failed to initialize ({exc})")
            continue

        sig["paramfile"] = rel  # repo-root-relative; the test re-initializes it.
        sig["name"] = name
        out = os.path.join(OUT_DIR, name + ".json")
        with open(out, "w") as fh:
            json.dump(sig, fh, indent=2, sort_keys=True)
            fh.write("\n")
        nfields = sum(len(g["fields"]) for g in sig["grids"])
        print(f"OK    {name}: problem_type={sig['problem_type']} "
              f"num_grids={sig['num_grids']} fields={nfields} -> "
              f"{os.path.relpath(out, ENZOMODULES)}")
        written += 1

    print(f"\nWrote {written} signature fixture(s) to "
          f"{os.path.relpath(OUT_DIR, ENZOMODULES)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
