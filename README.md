# PICs-atom-manipulation

Simulations of the silicon-nitride photonic-integrated-circuit (PIC) elements of the ANRF ARG proposal
ANRF/ARG/2026/009937 (Phase A, Objective 3): the SiN waveguide platform at 894.6 and 1064.4 nm, edge couplers,
MMI splitter trees, MZI switches and the 16-channel 894.6 nm router, thermo-optic phase shifters, and the
optical-phased-array (OPA) steering engine with its relay to the meta-objective and the atoms.

See `docs/findings.md` for the results and `docs/parameters.md` for the provenance of every input.

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
| `picsim/opa/` | OPA near field, Fourier-plane slit, pupil relay, multi-tweezer synthesis, phase noise |
| `picsim/transport.py` | motional excitation of a moved atom (harmonic trap, staircase trajectories) |
| `picsim/budget.py` | error-budget links from optical specs to atom metrics |
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
