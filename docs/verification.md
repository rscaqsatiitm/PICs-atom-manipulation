# Verification record

This note records the verification practices adopted on 2026-11-20 (rules in `CLAUDE.md`) and the results of
their first application. No stored result and no number in the report or the other docs was changed by this
work. Failing quantities are listed as findings in Sec. 5; how they are resolved is a separate decision.

## 1. Number gate

`report/check_numbers.py` previously listed report numbers not found in `docs/*.md`, `README.md` or
`results/*.json`, for review by hand, and always exited 0. It is now a gate (commit `d73408f`):
- exit 0 if every number in `report/sections/*.tex` is found or allowlisted;
- exit 1 if any number is unmatched;
- exit 2 if an entry of `report/number_allowlist.json` lacks a file, a token or a non-empty reason.
Stale allowlist entries print a warning. The allowlist is empty: 0 numbers were unmatched when the gate was
introduced. `make check-numbers` runs it, `make check` runs it with the tests, `make report` depends on it, and
the GitHub Actions workflow `.github/workflows/ci.yml` runs it with `make test` on every push and pull request.

Deliberate break (`tests/test_check_numbers.py`, run in CI): a copy of `03_platform.tex` with a planted
value in mW (the token is in the test, and must occur nowhere else in the repository) exits 1 and names the
number; the same copy with an allowlist entry for it exits 0; an entry whose
reason is blank exits 2. The unmodified report exits 0 with "unmatched: 0".

A first full test run after this note was written failed `test_planted_number_fails` (exit 0 instead of 1),
because an earlier draft of this note quoted the planted token, so the gate found it in `docs/`. The token was
removed from the note, and a precondition test, `test_planted_token_occurs_nowhere`, now checks that the token
occurs nowhere in the searched corpus. Deliberate break: appending the token to this note fails both tests;
removing it passes both.

## 2. MEEP results: doubled-resolution acceptance

**Tolerances.** Fixed and committed in `9b9b12d` before any run (`scripts/m32_meep_convergence.py`,
`TOLERANCES`); none was changed afterwards.

| module | quantity | tolerance |
|---|---|---|
| m3 grating | directionality_c, up_c, down_c, trans_c | abs 0.01 |
| m3 grating | theta_air_deg | abs 0.1 deg |
| m3 grating | L_1e_power_um | rel 10% |
| m14 facet | R, T, T_within_s0.9, fraction_within_s0.9, NA_1e2 | abs 0.01 |
| m10 star coupler | arm-field fidelity between resolutions | >= 0.99 |
| m10 star coupler | capture | abs 0.01 |

**Resolutions.** m3 and m14: stored production results at 40/um against new runs at 80/um. m10: the production
pair 20 against 40/um was not run (Sec. 3); a cheaper doubled pair was defined and committed before its run
(`ba4b65d`): a new run at 15/um against the stored 30/um result, coarser than production.

**Comparator self-test** (planted truth). For each module, identical inputs give all PASS and a planted shift
larger than the tolerance gives FAIL (`results/m32_meep_convergence.json`, `selftest`; `tests/test_meep_convergence.py`).
All three self-tests detected the planted shift.

The three `results/m32_*.json` files were first written by the script's own writer, without the git revision and
provenance-tagged inputs that `CLAUDE.md` requires. They were rewritten into the `save_result` layout without
rerunning MEEP: the outputs are identical to the first files, `git` records the revision at run time, and a rerun
of `check` on the stored files reproduces the comparison table exactly.

**Results** (`results/m32_meep_convergence.json`; new outputs in `results/m32_meep_highres.json`). Low is the
stored result, high the doubled resolution; for m10, low is the new 15/um run and high the stored 30/um run.

| case | quantity | low | high | difference | verdict |
|---|---|---|---|---|---|
| m3 etch 0.1, BOX 3.0, Si | directionality_c | 0.71535 | 0.70194 | -0.0134 | FAIL |
| | up_c | 0.37823 | 0.35761 | -0.0206 | FAIL |
| | down_c | 0.15051 | 0.15185 | +0.00135 | PASS |
| | trans_c | 0.46421 | 0.48309 | +0.0189 | FAIL |
| | theta_air_deg | 9.651 | 9.6543 | +0.0032 | PASS |
| | L_1e_power_um | 38.986 | 41.122 | +2.14 (5.5%) | PASS |
| m3 etch 0.2, BOX 2.0, Si | directionality_c | 0.75905 | 0.7572 | -0.00185 | PASS |
| | up_c | 0.67147 | 0.65824 | -0.0132 | FAIL |
| | down_c | 0.21314 | 0.21107 | -0.00208 | PASS |
| | trans_c | 0.031896 | 0.035945 | +0.00405 | PASS |
| | theta_air_deg | 8.9426 | 8.9715 | +0.0288 | PASS |
| | L_1e_power_um | 8.7606 | 9.0754 | +0.315 (3.6%) | PASS |
| m3 etch 0.05, BOX 2.65, Si | directionality_c | 0.70237 | 0.68609 | -0.0163 | FAIL |
| | up_c | 0.048295 | 0.072425 | +0.0241 | FAIL |
| | down_c | 0.020466 | 0.033138 | +0.0127 | FAIL |
| | trans_c | 0.93095 | 0.89388 | -0.0371 | FAIL |
| | theta_air_deg | 10.263 | 10.032 | -0.23 | FAIL |
| | L_1e_power_um | 422.7 | 269.59 | -153 (-36%) | FAIL |
| m3 etch 0.05, BOX 1.9, Bragg 2 pairs | directionality_c | 0.98496 | 0.98623 | +0.00127 | PASS |
| | up_c | 0.21235 | 0.21186 | -0.00048 | PASS |
| | down_c | 0.0032429 | 0.0029592 | -0.00028 | PASS |
| | trans_c | 0.78311 | 0.78365 | +0.00054 | PASS |
| | theta_air_deg | 10.111 | 9.9565 | -0.155 | FAIL |
| | L_1e_power_um | 123.71 | 124.06 | +0.349 (0.3%) | PASS |
| m3 etch 0.05, BOX 1.9, Bragg 3 pairs | directionality_c | 0.99742 | 0.99771 | +0.00029 | PASS |
| | up_c | 0.21026 | 0.2172 | +0.00694 | PASS |
| | down_c | 0.00054386 | 0.00049932 | -0.000045 | PASS |
| | trans_c | 0.78771 | 0.78066 | -0.00705 | PASS |
| | theta_air_deg | 10.263 | 9.8809 | -0.382 | FAIL |
| | L_1e_power_um | 126.74 | 122.14 | -4.6 (-3.6%) | PASS |
| m3 etch 0.05, BOX 1.9, metal | directionality_c | 1 | 1 | -1.1e-8 | PASS |
| | up_c | 0.21996 | 0.21425 | -0.00572 | PASS |
| | down_c | 1.0316e-7 | 1.028e-7 | -3.7e-10 | PASS |
| | trans_c | 0.7788 | 0.78414 | +0.00535 | PASS |
| | theta_air_deg | 9.9598 | 9.8809 | -0.0789 | PASS |
| | L_1e_power_um | 120.97 | 124.37 | +3.4 (2.8%) | PASS |
| m14 air | R | 0.1283 | 0.13909 | +0.0108 | FAIL |
| | T | 0.79069 | 0.78883 | -0.00187 | PASS |
| | T_within_s0.9 | 0.77172 | 0.77036 | -0.00136 | PASS |
| | fraction_within_s0.9 | 0.976 | 0.97659 | +0.00059 | PASS |
| | NA_1e2 | 0.96667 | 0.96608 | -0.00059 | PASS |
| m14 glass | R | 0.023256 | 0.025426 | +0.00217 | PASS |
| | T | 0.95907 | 0.95816 | -0.00091 | PASS |
| | T_within_s0.9 | 0.85184 | 0.85181 | -0.00002 | PASS |
| | fraction_within_s0.9 | 0.88819 | 0.88901 | +0.00082 | PASS |
| | NA_1e2 | 1.1581 | 1.1548 | -0.00331 | PASS |
| m10 port 0, 15 vs 30 (coarser than production 20) | arm-field fidelity | 1 | 0.89948 | | FAIL |
| | capture | 0.34674 | 0.36895 | +0.0222 | FAIL |

15 of 48 quantities fail. For context, not as a pre-registered check: between the stored m10 port-0 results the
arm-field fidelity is 0.965 for 15 against 20/um and 0.982 for 20 against 30/um (the latter already stated in
`docs/fpr_design.md`), so none of the available m10 pairs reaches 0.99.

## 3. Timing

Coarse runs at half the production resolution (serial, `results/m32_timing.json`), extrapolated as
resolution^3 for 2D FDTD (commit `21328e5`). The production runs used 4 processes in parallel.

| case | coarse | predicted (serial) | actual (4 in parallel) | actual / predicted |
|---|---|---|---|---|
| m3, six cases | 27.3 to 30.3 s at 20/um | 1747 to 1937 s at 80/um | 2303 to 2815 s | 1.32 to 1.54 |
| m14 air | 8.24 s at 20/um | 527 s at 80/um | 994 s | 1.89 |
| m14 glass | 8.37 s at 20/um | 536 s at 80/um | 1008 s | 1.88 |
| m10 port 0, 15/um | 833 s at 10/um | 2811 s (833 x 1.5^3) | 1672 s | 0.59 |
| m10 port 0, 40/um (production pair) | 833 s at 10/um | 53313 s (about 15 h) | not run | |

The m3 and m14 overruns are consistent with four processes sharing memory bandwidth; the m10 run was faster
than predicted, which suggests that the 10/um coarse run carried a fixed overhead not scaling with
resolution^3. The production m10 pair (20 against 40/um) was not run: at about 15 h on one core it exceeds the
4 h limit, and the user chose the cheaper pair.

## 4. Shared layers of the independent routes

A comparison between two routes cannot certify a layer they share. The layers are named next to each comparison.

| comparison | where | shared layers (not certified by the comparison) |
|---|---|---|
| BPM against 2D FDTD (m10 star coupler) | `docs/fpr_design.md`, "Shared layers of the BPM-FDTD comparison" | geometry builder `FPRLayout`; calibrated lateral indices N_CORE, N_CLAD (from the femwell mode solver and `picsim/materials.py`); launch and projection fields (`narrow_mode`, `tilted`); the 2D effective-index reduction; the wavelength |
| Kirchhoff against 2D FDTD (m14 facet) | `docs/collimator.md` Sec. 3, "Shared layers of the Kirchhoff-FDTD comparison" | the fitted slab core index 1.965; the slab mode (`slab_mode`); the cladding index `N_OX`; the wavelength |

## 5. Findings

These are reported, not resolved. The numbers quoted in the report and docs are unchanged.

**F1. Grating directionality and emission fractions at 40/um (m3).** For the 0.1 um etch on 3 um BOX, the
directionality moves from 0.715 to 0.702 and the up and transmitted fractions by 0.02. For the weakest emitter,
the 0.05 um etch on 2.65 um BOX over Si, every quantity fails: directionality 0.702 to 0.686, 1/e power length
423 to 270 um. Affected text: `report/sections/10_combine.tex` (i), "at best 0.70";
`report/sections/03_platform.tex`, "0.72" for the 100 nm etch; `docs/planar_tile.md` and
`docs/findings.md`, the same values. The Bragg-mirror and metal-mirror directionalities (0.986, 0.998, 1.000)
pass. The rest of the BOX sweep of `results/m30_directionality.json` (17 points over Si, 15 with mirrors) was
not rerun, so the quoted range "0.28 to 0.70" and the emission-rate variation "by up to a factor of eight" are
not checked; in that sweep the 1/e length over Si changes from 163.8 to 422.7 to 237.7 um between BOX 2.60, 2.65
and 2.70 um, so the rate is steep in BOX exactly where the 2.65 um case fails.

**F2. Grating emission angle for the 50 nm etch (m3).** The angle moves by -0.155, -0.23 and -0.382 deg in three
of the four 0.05 um etch cases (tolerance 0.1 deg). In sine this is up to 6.5e-3, two orders of magnitude above
the emission-angle tolerance of 5.2e-5 quoted for the planar tile; that tolerance concerns uniformity along a
grating, not the absolute angle, but the absolute period for a target angle is not fixed by the 40/um runs.
Affected text: `report/sections/03_platform.tex`, "8.9 to 10.0 deg".

**F3. Facet reflection into air (m14).** R moves from 12.8% to 13.9% (difference 0.0108 against 0.01).
Transmission, aperture fraction and NA pass. Affected text: `report/sections/05_rows.tex`, "Facet loss",
"reflects 12.8%"; `docs/collimator.md`, the same value.

**F4. Star-coupler FDTD (m10).** The cheaper pair, 15 against 30/um, gives fidelity 0.899 and a capture
difference of 0.022; the stored 20 against 30/um pair gives 0.982. The FDTD side of the BPM-FDTD comparison is
not converged to the 0.99 level at the resolutions run, and the production pair 20 against 40/um was not run.
The report (`04_engine.tex`) and `docs/fpr_design.md` already attribute the residual BPM-FDTD difference to FDTD
discretisation; this finding quantifies it but does not certify the 20/um FDTD values (0.997, 0.985, 0.93).

**Interpretation (inference, not tested).** The 0.05 um etch spans two cells at 40/um and four at 80/um, so the
grating index contrast is resolved coarsely; the angle shift corresponds to a change of about 6.5e-3 in the
average slab index. The 2.65 um BOX case sits near the minimum of the substrate-interference modulation of the
emission rate, where small phase changes move the rate strongly, and only about 7% of the power radiates over
the 30 um simulated grating, so its 1/e length is poorly conditioned.

**Options for resolution** (none applied): (i) adopt the 80/um values for the failing quantities and check them
at 160/um (2D, about 8 times the 80/um time per case, so about 4 to 6 h per m3 case serial); (ii) keep the 40/um
values and quote the 40 to 80 difference as their discretisation uncertainty; (iii) for the weak-grating cases,
lengthen the simulated grating so that most of the power radiates before the fit. For m10, the production pair
20 against 40/um needs about 15 h on one core.
