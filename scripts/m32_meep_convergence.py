"""Module 32: acceptance of the MEEP results under doubled resolution.

Each reported MEEP quantity must be stable when the resolution is doubled, within a tolerance fixed in TOLERANCES
below before any run (the commit that adds this file precedes the runs). The low-resolution values are the stored
production results; the high-resolution runs never write to the production files.

Cases (the ones the report quotes):
  m3 grating (scripts/m3_grating_fdtd.run_case): stored at 40 px/um, rerun at 80
  m14 facet (scripts/m14_facet_meep.main): air and glass, stored at 40, rerun at 80
  m10 star coupler (scripts/m10_fpr_meep.run): port 0, stored at 20, rerun at 40

Parts (run with the MEEP interpreter, PYTHONPATH=.):
  selftest  the comparator must PASS identical values and FAIL a copy with one value pushed beyond its tolerance
  timing    coarse run of every case at half the stored resolution, timed; production time predicted as
            t_coarse x (r_high / r_coarse)^3 (2D: pixels x time steps)
  run       production runs at the doubled resolution (parallel), raw values to results/m32_meep_highres.json
  check     comparison table with differences and verdicts -> results/m32_meep_convergence.json

A FAIL is a finding to report; no stored result or report number is changed by this module.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

# ---------------------------------------------------------------------------------------------- fixed before running
# kind "abs": |high - low| <= tol; "rel": |high - low| / |low| <= tol; "min": value(high vs low) >= tol
TOLERANCES = {
    "m3": {"directionality_c": ("abs", 0.01), "up_c": ("abs", 0.01), "down_c": ("abs", 0.01),
           "trans_c": ("abs", 0.01), "theta_air_deg": ("abs", 0.1), "L_1e_power_um": ("rel", 0.10)},
    "m14": {"R": ("abs", 0.01), "T": ("abs", 0.01), "T_within_s0.9": ("abs", 0.01),
            "fraction_within_s0.9": ("abs", 0.01), "NA_1e2": ("abs", 0.01)},
    "m10": {"fidelity": ("min", 0.99), "capture": ("abs", 0.01)},
}
RES = {"m3": (40, 80), "m14": (40, 80), "m10": (20, 40)}
M3_CASES = [  # (label, etch um, BOX um, reflector, Bragg pairs), stored file
    ("etch 0.1, BOX 3.0, Si (stored regression case)", 0.1, 3.0, "si", 0, "m3_grating_fdtd_raw.json"),
    ("etch 0.2, BOX 2.0, Si (directionality 0.76)", 0.2, 2.0, "si", 0, "m3_grating_fdtd_raw.json"),
    ("etch 0.05, BOX 2.65, Si (m30 best substrate)", 0.05, 2.65, "si", 0, "m30_directionality.json"),
    ("etch 0.05, BOX 1.9, Bragg 2 pairs", 0.05, 1.9, "dbr", 2, "m30_directionality.json"),
    ("etch 0.05, BOX 1.9, Bragg 3 pairs", 0.05, 1.9, "dbr", 3, "m30_directionality.json"),
    ("etch 0.05, BOX 1.9, metal", 0.05, 1.9, "metal", 0, "m30_directionality.json"),
]
M14_CASES = ["air", "glass"]
M10_PORTS = [0]


# ---------------------------------------------------------------------------------------------- comparator
def fidelity(a, b):
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))


def compare(module, low, high):
    """Rows {quantity, low, high, diff, tolerance, verdict} for one case; low/high are output dicts."""
    rows = []
    for q, (kind, tol) in TOLERANCES[module].items():
        if q == "fidelity":
            f = fidelity(np.array(low["T_re"]) + 1j * np.array(low["T_im"]),
                         np.array(high["T_re"]) + 1j * np.array(high["T_im"]))
            rows.append({"quantity": "arm-field fidelity", "low": 1.0, "high": f, "diff": 1 - f,
                         "tolerance": f">= {tol}", "verdict": "PASS" if f >= tol else "FAIL"})
            continue
        lo, hi = float(low[q]), float(high[q])
        d = hi - lo
        ok = abs(d) <= tol if kind == "abs" else abs(d) / max(abs(lo), 1e-300) <= tol
        rows.append({"quantity": q, "low": lo, "high": hi, "diff": d,
                     "tolerance": f"{'abs' if kind == 'abs' else 'rel'} {tol}", "verdict": "PASS" if ok else "FAIL"})
    return rows


def _stored_m3(case):
    label, etch, tb, kind, pairs, fname = case
    d = json.loads((RESULTS / fname).read_text())
    rows = d["outputs"]["rows"] if isinstance(d, dict) else d
    for r in rows:
        if abs(r["etch_um"] - etch) < 1e-9 and abs(r["t_box_um"] - tb) < 1e-9 and r.get("kind", "si") == kind \
                and r.get("dbr_pairs", 0) == pairs and not r.get("stack"):
            return r
        if fname == "m3_grating_fdtd_raw.json" and abs(r["etch_um"] - etch) < 1e-9 and abs(r["t_box_um"] - tb) < 1e-9:
            return r
    raise KeyError(label)


def stored():
    out = {"m3": {c[0]: _stored_m3(c) for c in M3_CASES}}
    out["m14"] = {m: json.loads((RESULTS / f"m14_facet_meep_{m}.json").read_text())["outputs"] for m in M14_CASES}
    out["m10"] = {f"port {p}": json.loads((RESULTS / f"m10_meep_port{p}_res{RES['m10'][0]}.json").read_text())
                  for p in M10_PORTS}
    return out


# ---------------------------------------------------------------------------------------------- MEEP runs
def _m3(case, res):
    import m3_grating_fdtd as G
    from m30_planar_tile import DBR
    G.RES = res
    _, etch, tb, kind, pairs, _ = case
    t0 = time.time()
    r = G.run_case(etch, tb, stack=DBR(pairs) if kind == "dbr" else None, substrate="metal" if kind == "metal" else "si")
    return {k: v for k, v in r.items() if not isinstance(v, list)}, time.time() - t0


def _m14(medium, res):
    import m14_facet_meep as F
    t0 = time.time()
    r = F.main(medium, res, write=False)
    return {k: v for k, v in r.items() if not k.startswith("spectrum")}, time.time() - t0


def _m10(port, res):
    import m10_fpr_meep as B
    t0 = time.time()
    r = B.run(port, res, write=False)
    return r, time.time() - t0


def _job(args):
    module, key, res = args
    if module == "m3":
        case = next(c for c in M3_CASES if c[0] == key)
        out, dt = _m3(case, res)
    elif module == "m14":
        out, dt = _m14(key, res)
    else:
        out, dt = _m10(int(key.split()[1]), res)
    print(json.dumps({"module": module, "case": key, "res": res, "seconds": round(dt, 1)}), flush=True)
    return {"module": module, "case": key, "res": res, "seconds": dt, "outputs": out}


def _jobs(res_index):
    js = [("m3", c[0], RES["m3"][res_index]) for c in M3_CASES]
    js += [("m14", m, RES["m14"][res_index]) for m in M14_CASES]
    js += [("m10", f"port {p}", RES["m10"][res_index]) for p in M10_PORTS]
    return js


def _save(name, payload):
    payload = {"name": name, "script": "scripts/m32_meep_convergence.py", **payload}
    (RESULTS / f"{name}.json").write_text(json.dumps(payload, indent=1, default=float))


def selftest():
    s = stored()
    out = {}
    for module, cases in s.items():
        key, low = next(iter(cases.items()))
        same = compare(module, low, low)
        bad = dict(low)
        if module == "m10":
            T = np.array(low["T_re"]) + 1j * np.array(low["T_im"])
            T = T * np.exp(1j * np.linspace(0, 2.0, len(T)))  # phase ramp, fidelity well below 0.99
            bad["T_re"], bad["T_im"] = T.real.tolist(), T.imag.tolist()
        else:
            q, (kind, tol) = next(iter(TOLERANCES[module].items()))
            bad[q] = float(low[q]) + (2 * tol if kind == "abs" else 2 * tol * abs(float(low[q])))
        broken = compare(module, low, bad)
        out[module] = {"case": key, "identical_all_pass": all(r["verdict"] == "PASS" for r in same),
                       "perturbed_any_fail": any(r["verdict"] == "FAIL" for r in broken), "perturbed_rows": broken}
        print(module, "identical -> all PASS:", out[module]["identical_all_pass"], "| perturbed -> FAIL detected:",
              out[module]["perturbed_any_fail"], flush=True)
    return out


def timing():
    """Coarse runs at half the stored resolution, one at a time (so the times are not inflated by contention)."""
    rows, total = [], 0.0
    for module, key, res_lo in _jobs(0):
        r_coarse = res_lo // 2
        r = _job((module, key, r_coarse))
        r_high = RES[module][1]
        pred = r["seconds"] * (r_high / r_coarse) ** 3
        total += pred
        rows.append({"module": module, "case": key, "res_coarse": r_coarse, "seconds_coarse": r["seconds"],
                     "res_production": r_high, "predicted_seconds": pred})
        print(f"{module} {key}: {r['seconds']:.1f} s at {r_coarse}/um -> predicted {pred / 60:.1f} min at {r_high}/um",
              flush=True)
    ncpu = int(os.environ.get("NPROC", 4))
    print(f"total predicted {total / 3600:.2f} h serial, about {total / 3600 / ncpu:.2f} h on {ncpu} processes",
          flush=True)
    _save("m32_timing", {"rule": "2D FDTD time ~ resolution^3 (pixels x steps)", "rows": rows,
                         "total_predicted_hours_serial": total / 3600, "processes": ncpu})
    return rows


def run():
    from multiprocessing import Pool
    js = sorted(_jobs(1), key=lambda j: -_predicted(j))
    with Pool(int(os.environ.get("NPROC", 4))) as pool:
        rows = pool.map(_job, js, chunksize=1)
    _save("m32_meep_highres", {"rows": rows})
    return rows


def _predicted(job):
    p = RESULTS / "m32_timing.json"
    if not p.exists():
        return 0.0
    for r in json.loads(p.read_text())["rows"]:
        if r["module"] == job[0] and r["case"] == job[1]:
            return r["predicted_seconds"]
    return 0.0


def check():
    s = stored()
    hi = json.loads((RESULTS / "m32_meep_highres.json").read_text())["rows"]
    tim = {(r["module"], r["case"]): r for r in json.loads((RESULTS / "m32_timing.json").read_text())["rows"]}
    table, fails = [], []
    for h in hi:
        low = s[h["module"]][h["case"]]
        rows = compare(h["module"], low, h["outputs"])
        t = tim.get((h["module"], h["case"]), {})
        table.append({"module": h["module"], "case": h["case"], "res_low": RES[h["module"]][0], "res_high": h["res"],
                      "seconds_high": h["seconds"], "predicted_seconds": t.get("predicted_seconds"), "rows": rows})
        for r in rows:
            print(f"{h['module']:4s} {h['case'][:44]:44s} {r['quantity']:22s} {r['low']:.5g} -> {r['high']:.5g} "
                  f"(diff {r['diff']:+.3g}, {r['tolerance']}) {r['verdict']}", flush=True)
            if r["verdict"] == "FAIL":
                fails.append({"module": h["module"], "case": h["case"], **r})
    _save("m32_meep_convergence", {"tolerances": TOLERANCES, "resolutions": RES, "selftest": selftest(),
                                   "table": table, "failures": fails,
                                   "notes": "low = stored production result; high = doubled resolution. FAIL rows "
                                            "are findings; no stored result or report number was changed."})
    print("failures:", len(fails), flush=True)
    return fails


if __name__ == "__main__":
    for part in sys.argv[1:] or ["selftest"]:
        globals()[part]()
