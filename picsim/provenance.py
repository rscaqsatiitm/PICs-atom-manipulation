"""Provenance-tagged parameters and result writing.

Provenance types follow the project's drafting rules:

USER      supplied in the proposal text (ANRF/ARG/2026/009937, draft v97)
SOURCE    retrieved and read in a working session; ``ref`` holds DOI/arXiv and location
DERIVED   calculated from other parameters; ``ref`` states the formula
SIM       produced by a simulation in this repository; ``ref`` names the script
CHOICE    a design or sweep choice made for the simulation, not a claim about the world
RECALLED  value recalled from memory, not verified in session; must be checked before use
"""
from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

KINDS = ("USER", "SOURCE", "DERIVED", "SIM", "CHOICE", "RECALLED")

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


@dataclass(frozen=True)
class P:
    """A physical parameter with unit and provenance."""

    value: float
    unit: str
    kind: str
    ref: str

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"unknown provenance kind {self.kind!r}")

    def __float__(self):
        return float(self.value)


def _jsonable(x):
    if isinstance(x, P):
        return asdict(x)
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, complex):
        return {"re": x.real, "im": x.imag}
    return x


def git_rev() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "uncommitted"


def save_result(name: str, *, inputs: dict, outputs: dict, script: str, fig=None, notes: str = ""):
    """Write results/<name>.json (inputs with provenance, outputs tagged SIM) and optional figure."""
    RESULTS.mkdir(exist_ok=True)
    record = {
        "name": name,
        "script": script,
        "git": git_rev(),
        "notes": notes,
        "inputs": _jsonable(inputs),
        "outputs": _jsonable(outputs),
        "outputs_kind": "SIM (review before submission)",
    }
    (RESULTS / f"{name}.json").write_text(json.dumps(record, indent=1))
    if fig is not None:
        fig.savefig(RESULTS / f"{name}.png", dpi=150, bbox_inches="tight")
    return record
