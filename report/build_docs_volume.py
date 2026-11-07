"""Convert docs/*.md into LaTeX chapters for the compiled-docs volume (report/docs_volume.tex).

Each document becomes one chapter (its H1 title), in topic order. Wording is kept; only the conversion and a few
mechanical fixes are applied: pandoc gfm -> latex with chapters as the top level, \\url for bare paths is not
added, and very wide tables are set in footnote size.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "report" / "docs"
ORDER = ["findings", "parameters", "twenty_tweezers", "independent_steering", "fpr_design", "scaling",
         "scaling_beyond", "stacked_relay", "collimator", "row_coherence", "pm_am", "transport3d", "rearrangement",
         "repair", "megascale"]


def convert(name: str) -> str:
    md = ROOT / "docs" / f"{name}.md"
    tex = subprocess.run(["pandoc", "-f", "gfm", "-t", "latex", "--top-level-division=chapter", f"--id-prefix={name}-", str(md)],
                         check=True, capture_output=True, text=True).stdout
    # the docs number their own sections ("## 3. Platform"); LaTeX numbers them again, so drop the doc's number
    tex = re.sub(r"\\(section|subsection|subsubsection)\{\d+(\.\d+)*\.?\s+", r"\\\1{", tex)
    # pandoc 3.1 ignores --id-prefix for gfm headings: prefix every label and internal link with the document name
    tex = re.sub(r"\\label\{(?!doc:)", "\\\\label{" + name + "-", tex)
    tex = re.sub(r"\\hyperref\[", "\\\\hyperref[" + name + "-", tex)
    # label the chapter after the file, keep the doc's own title
    tex = re.sub(r"\\chapter\{(.*?)\}\\label\{[^}]*\}", lambda m: "\\chapter{" + m.group(1) + "}\\label{doc:" + name + "}",
                 tex, count=1)
    if "\\label{doc:" + name + "}" not in tex:
        tex = re.sub(r"\\chapter\{(.*?)\}", lambda m: "\\chapter{" + m.group(1) + "}\\label{doc:" + name + "}", tex, count=1)
    blocks = re.split(r"(\\begin\{longtable\}.*?\\end\{longtable\})", tex, flags=re.S)
    out = []
    for b in blocks:
        out.append(fit_table(b) if b.startswith("\\begin{longtable}") else b)
    return "".join(out)


def fit_table(b: str) -> str:
    """Give a pandoc longtable wrapping p-columns sized by content when its natural width exceeds the line."""
    head, rest = b.split("\n", 1)
    m = re.match(r"\\begin\{longtable\}\[\]\{@\{\}(l+)@\{\}\}", head)
    if not m:
        return "{\\small\n" + b + "\n}"
    n = len(m.group(1))
    body = re.sub(r"\\(toprule|midrule|bottomrule|endhead|endlastfoot|endfirsthead)(\\noalign\{\})?", "", rest)
    body = body.replace("\\end{longtable}", "")
    rows = [r for r in re.split(r"\\\\\s*\n", body) if "&" in r]
    width = [1] * n
    for r in rows:
        cells = r.replace("\n", " ").split(" & ")
        for k, c in enumerate(cells[:n]):
            txt = re.sub(r"\\[a-zA-Z]+\{?|[{}]", "", c).strip()
            width[k] = max(width[k], min(len(txt), 60), int(2.8 * max((len(w) for w in txt.split()), default=1)))
    size = "\\footnotesize" if n > 5 or sum(width) > 140 else "\\small"
    if sum(width) <= 85:
        return "{" + size + "\n" + b + "\n}"
    tot = sum(width)
    frac = [max(w / tot, 0.06) for w in width]
    s_ = sum(frac)
    avail = f"(\\linewidth-{2 * n}\\tabcolsep)"
    spec = "".join(f">{{\\raggedright\\arraybackslash}}p{{{f / s_:.3f}\\dimexpr{avail}\\relax}}" for f in frac)
    return "{" + size + "\n\\begin{longtable}[]{@{}" + spec + "@{}}\n" + rest + "\n}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = []
    for name in ORDER:
        (OUT / f"{name}.tex").write_text(convert(name))
        inputs.append(f"\\input{{docs/{name}.tex}}")
    (OUT / "all.tex").write_text("\n".join(inputs) + "\n")
    print("wrote", len(ORDER), "chapters")


if __name__ == "__main__":
    main()
