# PICs-atom-manipulation

Simulations of the silicon-nitride photonic-integrated-circuit (PIC) elements of the ANRF ARG proposal
ANRF/ARG/2026/009937 (Phase A, Objective 3): the SiN waveguide platform at 894.6 and 1064.4 nm, edge couplers,
MMI splitter trees, MZI switches and the 16-channel 894.6 nm router, thermo-optic phase shifters, and the
optical-phased-array (OPA) steering engine with its relay to the meta-objective and the atoms.

See `docs/findings.md` for the results, `docs/twenty_tweezers.md` for the 20-tweezer design,
`docs/independent_steering.md` for independently steerable tweezers, `docs/fpr_design.md` for the
beam-propagation design and optimisation of that engine, and
`docs/parameters.md` for the provenance of every input.

## Status (first pass complete)

| Module | Script | Main outputs |
|---|---|---|
| 1 Platform | `scripts/m1_platform.py` | design widths, n_eff, n_g, dn/dT, bend loss, substrate leakage, coupling length |
| 2 Couplers, MMI | `scripts/m2_couplers_mmi.py` | taper overlap and tolerance; MMI loss, bandwidth, split-ratio statistics |
| 3 OPA | `scripts/m3_opa.py`, `scripts/m3_grating_fdtd.py` | range vs channels, pitch, gratings (FDTD), multi-tweezer synthesis, noise, update rate, loss |
| 4 Switches, router | `scripts/m4_switches.py` | extinction Monte Carlo, polarisation, actuator dynamics, 16-channel router |
| 5 Heaters | `scripts/m5_heaters.py` | P_pi, rise time, thermal crosstalk; actuator comparison |
| 6 Power | `scripts/m6_power.py` | self-heating, facet intensity, Kerr phase |
| 7 Budget | `scripts/m7_budget.py` | error-budget links, `docs/parameters.md` |
| 8 Twenty tweezers | `scripts/m8_twenty.py` | star-coupler OPA with 20 switched ports vs phase synthesis; `docs/twenty_tweezers.md` |
| 9 Independent steering | `scripts/m9_independent.py` | sliding-window focal-plane OPA vs AOD-emulation routes; `docs/independent_steering.md` |
| 10 Star-coupler BPM | `scripts/m10_fpr.py`, `scripts/m10_fpr_meep.py` | calibrated effective index, wide-angle BPM of the dense-port coupler, FDTD check, least-squares synthesis; `docs/fpr_design.md` |
| 11 Geometry optimisation | `scripts/m11_optimise.py` | oversampling, emitter facet, slab widths, window width |
| 12 Scaling | `scripts/m12_scaling.py` | tweezers per row (packing), aberration vs size, scaling table; `docs/scaling.md` |
| 13 Stacked relay | `scripts/m13_stack.py` | stacked-chip anamorphic relay: facet far field, y fill, 2D row crosstalk (coherent / incoherent), collimator tolerances, layout; `docs/stacked_relay.md` |

## Layout

| path | content |
|---|---|
| `picsim/params.py` | every physical input, with provenance (USER, SOURCE, DERIVED, SIM, CHOICE, RECALLED) |
| `picsim/materials.py` | Sellmeier models: LPCVD Si3N4 (Luke 2015), SiO2 (Malitson 1965), AlN |
| `picsim/modes.py` | vector FEM mode solver wrapper (femwell): straight, bent, coupled and leaky guides |
| `picsim/slab.py` | analytic and finite-difference slab modes |
| `picsim/couplers.py` | inverse-taper edge couplers (mode overlap with Gaussian fibre modes) |
| `picsim/mmi.py` | MMI couplers by effective index + eigenmode expansion |
| `picsim/circuits.py` | MZI transfer matrices, Monte Carlo extinction, actuator dynamics, routing tree |
| `picsim/heaters.py` | 2D thermal FEM of heaters: P_pi, rise time, crosstalk |
| `picsim/opa/` | OPA near field, Fourier-plane slit, pupil relay, multi-tweezer synthesis, phase noise, star coupler |
| `picsim/transport.py` | motional excitation of a moved atom (harmonic trap, staircase trajectories) |
| `picsim/budget.py` | error-budget links from optical specs to atom metrics |
| `picsim/bpm.py` | 2D wide-angle FD-BPM and the calibrated lateral effective-index model |
| `scripts/m*_*.py` | one script per module; each writes `results/<name>.json` (+ `.png`) |
| `tests/` | unit tests against analytic limits |

## Running

```bash
pip install -e .[test]
```

```bash
make test
```

```bash
make all
```

The 2D FDTD of grating antennas uses MEEP (conda-forge only). Create `environment-meep.yml` with micromamba and
pass its interpreter: `make m3g MEEPPY=/path/to/env/bin/python`. Without it, the OPA model still runs; the
grating results are read from `results/m3_grating_fdtd_raw.json` if present.

`gmsh` needs `libGLU` on Linux (`apt-get install libglu1-mesa`).

## Conventions

Lengths in micrometres unless stated. Every `results/*.json` holds the inputs with provenance and the outputs
(type SIM). Simulated values are valid only inside the simulated regime and must be reviewed before they are
quoted in the proposal.
