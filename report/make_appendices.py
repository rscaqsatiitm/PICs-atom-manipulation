"""Generate the appendices of report/report.tex from the repository record (no hand copying of numbers).

A: parameter ledger (docs/parameters.md, all non-SIM rows; SIM rows chain modules and are in Volume II).
B: reproduction map (README module table).
C: open placeholders and recalled values (the "Open placeholders" and "Recalled" paragraphs of the docs, and the
   RECALLED rows of their verification tables).
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "report"))
from build_docs_volume import fit_table  # noqa: E402

SEC = ROOT / "report" / "sections"


def md2tex(md: str) -> str:
    tex = subprocess.run(["pandoc", "-f", "gfm", "-t", "latex"], input=md, check=True, capture_output=True,
                         text=True).stdout
    blocks = re.split(r"(\\begin\{longtable\}.*?\\end\{longtable\})", tex, flags=re.S)
    return "".join(fit_table(b) if b.startswith("\\begin{longtable}") else b for b in blocks)


def table_rows(md: str):
    return [l for l in md.splitlines() if l.startswith("|")]


def params():
    lines = table_rows((ROOT / "docs" / "parameters.md").read_text())
    head, sep, body = lines[0], lines[1], lines[2:]
    keep = [l for l in body if "| SIM |" not in l]
    md = "\n".join([head, sep] + keep)
    out = ("\\section{Parameter ledger}\\label{app:params}\n"
           "Inputs of every simulation, with their provenance, as recorded in the \\texttt{inputs} of "
           "\\texttt{results/*.json} and collected in \\texttt{docs/parameters.md} by \\texttt{scripts/m7\\_budget.py}. "
           f"Rows of type SIM (outputs of one module used as inputs of another) are omitted here ({len(body) - len(keep)} rows); "
           "they are listed in Volume~II. USER: proposal text; SOURCE: read in session; DERIVED: calculation; "
           "CHOICE: stated design or sweep choice; RECALLED: not verified, to be checked before use.\n\n")
    (SEC / "A_params.tex").write_text(out + md2tex(md))


def repro():
    md = "\n".join(table_rows((ROOT / "README.md").read_text()))
    out = ("\\section{Reproduction map}\\label{app:repro}\n"
           "Each module, its script and its main outputs (from \\texttt{README.md}). Every script writes "
           "\\texttt{results/<name>.json} with its inputs (with provenance), outputs and the git revision; "
           "\\texttt{make test} runs the test suite.\n\n")
    (SEC / "B_repro.tex").write_text(out + md2tex(md))


def open_items():
    parts = []
    for p in sorted((ROOT / "docs").glob("*.md")):
        txt = p.read_text()
        items = []
        for para in re.split(r"\n\s*\n", txt):
            if re.match(r"\s*(\*\*)?(Open placeholders|Recalled, not verified)", para):
                items.append(" ".join(para.split()))
        for l in table_rows(txt):
            cells = [c.strip() for c in l.strip("|").split("|")]
            if len(cells) >= 4 and cells[2].upper().startswith("RECALLED"):
                items.append(f"RECALLED: {cells[0]} ({cells[1]}); {cells[3]}")
        if items:
            parts.append(f"**`docs/{p.name}`**\n\n" + "\n".join(f"- {i}" for i in items))
    md = "\n\n".join(parts)
    out = ("\\section{Open placeholders and recalled values}\\label{app:ledger}\n"
           "Values marked [VALUE?] or [REF?] and values recalled but not verified, collected from the "
           "verification notes of each document. The full provenance ledgers (every quoted value with its type and "
           "pointer) are at the end of each chapter of Volume~II.\n\n")
    (SEC / "C_open.tex").write_text(out + md2tex(md))


if __name__ == "__main__":
    params()
    repro()
    open_items()
    print("appendices written")
