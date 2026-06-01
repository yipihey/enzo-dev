"""Read/write the language-neutral EnzoModules ``.fixture`` format.

A fixture stores the captured inputs and reference outputs of a legacy
kernel.  The format is deliberately plain text so any rewrite -- in any
language -- can consume it without an HDF5/JLD stack:

    # comment
    key = value          # scalar: int if integral, else float, else string
    @name v1 v2 v3 ...   # array of floats

``tools/capture_twoshock.py`` and this module must agree on the format.
"""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence


@dataclass
class Fixture:
    """Parsed fixture: scalars (int/float/str) and float arrays."""

    scalars: Dict[str, Any] = field(default_factory=dict)
    arrays: Dict[str, List[float]] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        if key in self.arrays:
            return self.arrays[key]
        return self.scalars[key]

    def __contains__(self, key: str) -> bool:
        return key in self.arrays or key in self.scalars

    @property
    def name(self) -> str:
        return str(self.scalars.get("name", "<unnamed>"))


def _parse_scalar(s: str) -> Any:
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def load_fixture(path: str) -> Fixture:
    fx = Fixture()
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("@"):
                parts = line[1:].split()
                if not parts:
                    continue
                fx.arrays[parts[0]] = [float(p) for p in parts[1:]]
            elif "=" in line:
                key, val = line.split("=", 1)
                fx.scalars[key.strip()] = _parse_scalar(val.strip())
    return fx


def load_dir(directory: str) -> List[Fixture]:
    """Load every ``*.fixture`` in *directory*, sorted by filename."""
    files = sorted(glob.glob(os.path.join(directory, "*.fixture")))
    return [load_fixture(f) for f in files]


def _fmt(x: Any) -> str:
    if isinstance(x, bool):
        return str(int(x))
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return repr(x)            # round-trippable
    return str(x)


# Canonical key ordering so Julia/Python writers stay byte-identical.
SCALAR_ORDER = ["kernel", "name", "idim", "jdim", "i1", "i2", "j1", "j2",
                "ipresfree", "gravity", "idual", "dt", "gamma", "pmin", "eta1"]
ARRAY_ORDER = ["dls", "drs", "pls", "prs", "uls", "urs", "grslice",
               "pbar", "ubar"]


def save_fixture(path: str, fx: Fixture,
                 scalar_order: Sequence[str] = SCALAR_ORDER,
                 array_order: Sequence[str] = ARRAY_ORDER) -> str:
    written_scalars = set()
    written_arrays = set()
    with open(path, "w") as fh:
        fh.write("# EnzoModules fixture\n")
        for k in scalar_order:
            if k in fx.scalars:
                fh.write(f"{k} = {_fmt(fx.scalars[k])}\n")
                written_scalars.add(k)
        for k in fx.scalars:                       # any extras, sorted
            if k not in written_scalars:
                fh.write(f"{k} = {_fmt(fx.scalars[k])}\n")
        for k in array_order:
            if k in fx.arrays:
                fh.write("@" + k + " " + " ".join(_fmt(v) for v in fx.arrays[k]) + "\n")
                written_arrays.add(k)
        for k in fx.arrays:
            if k not in written_arrays:
                fh.write("@" + k + " " + " ".join(_fmt(v) for v in fx.arrays[k]) + "\n")
    return path
