"""report/check_numbers.py as a gate: planted numbers fail, allowlisted ones pass, reasons are required."""
import json
import subprocess
import sys

from picsim.provenance import REPO

SCRIPT = REPO / "report" / "check_numbers.py"
SECTION = REPO / "report" / "sections" / "03_platform.tex"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)


def planted(tmp_path):
    f = tmp_path / "03_platform.tex"
    f.write_text(SECTION.read_text() + "\nA planted value of 123.457~mW that occurs nowhere.\n")
    return f


def test_real_report_passes():
    r = run()
    assert r.returncode == 0, r.stdout
    assert "unmatched: 0" in r.stdout


def test_planted_number_fails(tmp_path):
    r = run(planted(tmp_path))
    assert r.returncode == 1
    assert "123.457" in r.stdout


def test_allowlisted_number_passes(tmp_path):
    allow = tmp_path / "allow.json"
    allow.write_text(json.dumps([{"file": "03_platform.tex", "token": "123.457", "reason": "planted for this test"}]))
    r = run("--allowlist", allow, planted(tmp_path))
    assert r.returncode == 0, r.stdout


def test_allowlist_entry_without_reason_is_rejected(tmp_path):
    allow = tmp_path / "allow.json"
    allow.write_text(json.dumps([{"file": "03_platform.tex", "token": "123.457", "reason": "  "}]))
    r = run("--allowlist", allow, planted(tmp_path))
    assert r.returncode == 2
