"""Gate: every number in the report must occur in docs/*.md, README.md or results/*.json.

A number is a token like 57, 0.973, 1.8e-4 (also written $1.8\\times10^{-4}$), 2,088. Integers below 10, years,
and section/figure numbering are skipped. A token counts as found if its normalised form appears in the corpus
(docs text, README, and all numbers in results JSON, rounded to the token's precision).

Exit status: 0 if every number is found or allowlisted, 1 if any number is unmatched, 2 if the allowlist is
malformed (an entry without a file, a token or a non-empty reason). The allowlist (report/number_allowlist.json)
is a list of {"file": <tex file name>, "token": <number as written>, "reason": <why it is not a result>}; entries
that no longer match anything are reported as stale (warning only).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ALLOWLIST = ROOT / "report" / "number_allowlist.json"


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


def load_allowlist(path):
    """Return (entries, malformed entries)."""
    path = Path(path)
    data = json.loads(path.read_text()) if path.exists() else []
    if not isinstance(data, list):
        return [], [data]
    bad = [e for e in data if not isinstance(e, dict) or not str(e.get("file", "")).strip()
           or not str(e.get("token", "")).strip() or not str(e.get("reason", "")).strip()]
    return data, bad


def main(files, allowlist=ALLOWLIST):
    allow, bad = load_allowlist(allowlist)
    if bad:
        print("malformed allowlist entries (file, token and a non-empty reason are required):")
        for e in bad:
            print("  ", e)
        return 2
    allowed = {(Path(e["file"]).name, norm(str(e["token"]))) for e in allow}
    used = set()
    toks, nums = corpus()
    miss = {}
    for f in files:
        for t in report_tokens(Path(f).read_text()):
            n = norm(t)
            if re.fullmatch(r"\d", n) or re.fullmatch(r"(19|20)\d\d", n):
                continue
            if not found(t, toks, nums):
                key = (Path(f).name, n)
                if key in allowed:
                    used.add(key)
                    continue
                miss.setdefault(f, []).append(t)
    for f, ts in miss.items():
        print(f, ":", " ".join(sorted(set(ts), key=ts.index)))
    stale = sorted(allowed - used)
    if stale:
        print("warning: allowlist entries that match nothing:", stale)
    n_miss = sum(len(set(v)) for v in miss.values())
    print("unmatched:", n_miss, "| allowlisted:", len(used))
    return 1 if n_miss else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="*", help="tex files (default: report/sections/[0-9]*.tex)")
    ap.add_argument("--allowlist", default=str(ALLOWLIST))
    a = ap.parse_args()
    sys.exit(main(a.files or sorted(str(p) for p in (ROOT / "report" / "sections").glob("[0-9]*.tex")), a.allowlist))
