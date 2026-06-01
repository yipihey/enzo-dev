"""EnzoModules -- extracted, fixture-certified, callable legacy Enzo kernels.

This package builds selected legacy C/Fortran Enzo kernels as a shared
library and exposes each as a documented Python function, together with
golden fixtures (captured legacy I/O) and a tolerance/diff utility.

It is deliberately rewrite-agnostic.  Its durable, language-neutral
artifacts -- the compiled bridge ``.so`` and the plain-text fixtures -- can
be reused by *any* port of Enzo (Python or Julia) to (a) call the reference
legacy implementation and (b) verify a rewrite against the captured outputs
using a shared tolerance policy.

See ``README.md`` for the wrap -> test -> capture workflow.
"""
from __future__ import annotations

import os

from . import bridge, diff, examples, fixtures, hydro
from .diff import BITWISE, CompareResult, Tolerance, compare, isclose
from .fixtures import Fixture, load_dir, load_fixture, save_fixture

__all__ = [
    "bridge", "diff", "examples", "fixtures", "hydro",
    "Tolerance", "BITWISE", "CompareResult", "compare", "isclose",
    "Fixture", "load_fixture", "load_dir", "save_fixture",
    "fixturedir", "has_library",
]

__version__ = "0.1.0"

_HERE = os.path.dirname(os.path.abspath(__file__))


def fixturedir(*parts: str) -> str:
    """Absolute path to the fixtures shipped with this package."""
    return os.path.join(_HERE, "..", "fixtures", *parts)


def has_library() -> bool:
    """True if the compiled pilot library is available (kernels runnable)."""
    return bridge.available()
