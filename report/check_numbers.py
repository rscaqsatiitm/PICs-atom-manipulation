"""List numbers in the report that do not occur in docs/*.md or results/*.json (for review by hand).

A number is a token like 57, 0.973, 1.8e-4 (also written $1.8\\times10^{-4}$), 2,088. Integers below 10, years,
and section/figure numbering are skipped. A token counts as found if its normalised form appears in the corpus
(docs text, README, and all numbers in results JSON, rounded to the token's precision).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def norm(tok: str) -> str:
    t = tok.replace(",", "")
    m = re.fullmatch(r"(\d+(?:\.\d+)?)e([+-]?\d+)", t)
    if m:
        return f"{float(m.group(1)):g}e{int(m.group(2))}"
    return t


def json_numbers(x, out):
    if isinstance(x, dict):
        for v in x.values():
            json_numbers(v, out)
    elif isinstance(x, list):
        for v in x:
            json_numbers(v, out)
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        out.append(float(x))


def corpus():
    text = " ".join(p.read_text() for p in (ROOT / "docs").glob("*.md")) + (ROOT / "README.md").read_text()
    toks = set(norm(t) for t in re.findall(r"\d[\d,]*(?:\.\d+)?(?:e[+-]?\d+)?", text))
    nums = []
    for p in (ROOT / "results").glob("*.json"):
        try:
            json_numbers(json.loads(p.read_text()), nums)
        except Exception:
            pass
    return toks, nums


def report_tokens(tex: str):
    tex = re.sub(r"%.*", "", tex)
    tex = re.sub(r"\\(label|ref|cref|Cref|cite|includegraphics|input)(\[[^\]]*\])?\{[^}]*\}", "", tex)
    tex = re.sub(r"\$?(\d+(?:\.\d+)?)\s*\\times\s*10\^\{?([+-−]?\d+)\}?\$?", lambda m: f"{m.group(1)}e{m.group(2)}", tex)
    tex = re.sub(r"\\(begin|end)\{[^}]*\}(\{[^}]*\})?|\\[a-zA-Z]+", " ", tex)
    return re.findall(r"(?<![\w.])\d[\d,]*(?:\.\d+)?(?:e[+-]?\d+)?", tex)


def found(tok, toks, nums):
    """Found in the docs text, or equal to a results number at the token's own precision."""
    n = norm(tok)
    if n in toks:
        return True
    try:
        v = float(n)
    except ValueError:
        return False
    if "e" in n:
        mant = n.split("e")[0]
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        target = f"{v:.{sig - 1}e}"
        return any(f"{x:.{sig - 1}e}" == target for x in nums if x)
    dec = len(n.split(".")[1]) if "." in n else 0
    return any(round(x, dec) == round(v, dec) for x in nums)


def main(files):
    toks, nums = corpus()
    miss = {}
    for f in files:
        for t in report_tokens(Path(f).read_text()):
            n = norm(t)
            if re.fullmatch(r"\d", n) or re.fullmatch(r"(19|20)\d\d", n):
                continue
            if not found(t, toks, nums):
                miss.setdefault(f, []).append(t)
    for f, ts in miss.items():
        print(f, ":", " ".join(sorted(set(ts), key=ts.index)))
    print("unmatched:", sum(len(set(v)) for v in miss.values()))


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(str(p) for p in (ROOT / "report" / "sections").glob("[0-9]*.tex")))
