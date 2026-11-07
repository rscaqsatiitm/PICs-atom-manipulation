# 3D transport with handoff between the static array and the moving tweezers

> Follow-up (`docs/rearrangement.md`): slower handoff steps (200 us ramps, 400 us lane steps) cut the heating per
> handoff to 0.12 uK at 0.1 / 0.3 mK (1D quantum, handoff and lane axis); with the x move it is about 0.3 uK per
> move (quantum, 0.12 + 0.17 uK), and the 3D classical run gives 0.16 uK.
> The shallow 0.3 mK point is then also the low-heating point. The static array's depth uniformity is not a
> constraint.

Follow-up to `docs/scaling_beyond.md`, whose 1D model left out handoff, the static array and the third
dimension. Code: `picsim/transport3d.py`, `scripts/m18_transport3d.py`; results `results/m18_*.json`,
`results/m18_trap.npz`; tests `tests/test_transport3d.py`. Values are SIM (classical Monte Carlo) and valid in the
regime stated. They should be reviewed before they are quoted.

## 1. Result

A moving tweezer of the adopted design picks an atom out of a static 1075 nm trap, carries it 12 or 40 um and
drops it into an empty static trap. With the right depths and path, it delivers 3,000 of 3,000 simulated atoms
(loss < 1e-3 at 95% confidence) and heats them by 1 to 2 uK per move. Three conditions are needed:

(i) **Move along lanes between the rows, not along the rows.** A tweezer that crosses empty static sites on its
row loses 45 to 95% of atoms at most settings. Only a moving trap five times deeper than the static traps, moving
slowly (400 us), recovers 99.1 to 99.7%. Along a lane 2 um off the row, reached with the row-select y tilt, the
static traps are 1.6e-5 of their depth away and transport is clean. All tweezers of a stack enter and leave the
lanes together.

(ii) **Moving depth at least three times the static hold depth.** At a ratio of 1.5, 5% of atoms stay behind; at
2, 3 to 5% are lost; at 3, 99.5 to 99.9% are delivered with about 10 uK of heating per move; at 5, 99.6 to 100%
with 1 to 2.6 uK.

(iii) **Hold depth of about 20 kT.** Atoms waiting next to a lane feel the passing tweezer's out-of-focus light
(6 to 15% of its peak at 1.5 to 2.5 um, a few um from focus). At a hold depth of 10 kT, 1 to 2% of them are lost
from the hot tail; at 20 kT, none of 3,000.

Operating points (12 um move in 400 us, 100 us ramps, 200 kHz updates, lane 2 um):

| atoms | hold | moving | delivered | heating per move | waiting atoms lost | power per tweezer | 720 tweezers |
|---|---|---|---|---|---|---|---|
| 5 uK | 0.1 mK | 0.3 mK | 3000/3000 | 9.9 uK | 0/3000 | 14 mW | 10 W |
| 5 uK | 0.1 mK | 0.5 mK | 3000/3000 | 1.8 uK | <= 2/3000 | 24 mW | 17 W |
| 10 uK | 0.2 mK | 0.6 mK | 3000/3000 | 10.0 uK | 0/3000 | 29 mW | 21 W |
| 10 uK | 0.2 mK | 1.0 mK | 3000/3000 | 0.9 uK | <= 3/3000 | 48 mW | 34 W |

(Power scales from 57 mW at 1.2 mK, `docs/collimator.md`.) Hence the required depth is set by the atom
temperature, not by the move: about 60 kT for the moving trap, or 100 kT for low heating. Cooling the atoms to
5 uK halves the laser power. The 1D estimate of 0.2 mK in `docs/scaling_beyond.md` is too optimistic once handoff
and waiting atoms are included.

## 2. Model

**Moving tweezer.** The intensity of one tweezer of design D1 (BPM response, least-squares synthesis), through the
Weierstrass + aplanat collimator and the circular NA 0.5 pupil, computed by scalar focusing with defocus on a grid
of 7 x 6 x 20 um (`m18_field`). FWHM 1.73 um (x, along the row), 1.37 um (y), 9.8 um (z); trap frequencies at 1 mK
54, 67 and 22 kHz. Its shape is held fixed while it moves; positions and depth are sampled and held at 200 kHz.

**Static array.** Gaussian beams at 1075 nm, waist 0.85 um (CHOICE), Rayleigh range 2.1 um, 4 um pitch, three
rows of 15 sites, all at one global hold depth (a fixed meta-optic array cannot switch single sites). The two
wavelengths add incoherently.

**Atoms.** Cs, classical, velocity Verlet with 0.25 us steps, gravity along y, 1,500 to 3,000 atoms per case. Each
starts in a thermal distribution (harmonic sampling), truncated to atoms bound in the real Gaussian site.

**Sequence.** Ramp the moving trap up on the source site (sin^2, 100 us); move 2 um in y to the lane (100 us);
move in x (minimum jerk); move back to the row (100 us); ramp down on the target (100 us); hold 100 us. An atom
counts as delivered if it is bound in the target static trap after the moving trap is off.

## 3. Tolerances

Hold 0.1 mK, moving 0.3 / 0.5 mK, 10 uK, lane, 12 um in 400 us (fraction delivered):

| parameter | value | moving 0.3 mK | moving 0.5 mK |
|---|---|---|---|
| none | | 0.9993 | 0.9993 |
| axial offset of the 1064 and 1075 nm foci | 0.5 / 1 / 2 / 4 um | 0.998 / 0.995 / 0.939 / 0.455 | 0.995 / 0.997 / 0.957 / 0.902 |
| lateral offset, y | 0.1 / 0.2 / 0.4 um | 1.000 / 0.998 / 0.999 | 1.000 / 1.000 / 0.998 |
| lateral offset, x | 0.1 / 0.2 / 0.4 um | 0.998 / 1.000 / 0.998 | 1.000 / 1.000 / 0.999 |
| ramp time | 25 / 200 us | 0.994 / 0.998 | 0.993 / 1.000 |
| move | 40 um in 400 / 800 us | 0.998 / 0.998 | 0.998 / 1.000 |
| updates | smooth instead of 200 kHz | 0.999 | 0.999 |
| gravity | off | 0.999 | 0.999 |

Three points follow. (i) **Axial registration** of the two foci is the one tight alignment, to about +/-1 um. A
diffractive objective shifts focus by about 170 um between 1064.4 and 1075 nm (16 um/nm, DERIVED), so the two
paths must be focused independently before they are combined, and held to 1 um. (ii) **Lateral registration** is
loose: 0.4 um costs nothing. The 16 nm row-placement budget of the stacked relay is therefore not a constraint
for transport. (iii) **Speed** is not limiting in this regime: 40 um in 400 us is as safe as 12 um.

## 4. On the row versus along a lane

Along the row, the moving trap passes over empty static sites. The total potential then develops a dip that
moves with it and kicks the atom. Losses depend strongly on the move time, with resonances where the transit time
over a site matches the trap period:

| hold / moving | 200 us | 400 us |
|---|---|---|
| 0.1 / 0.3 mK | 5 to 7% delivered | 54 to 57% |
| 0.1 / 0.5 mK | 42 to 48% | 99.1 to 99.7% (heating 14 uK) |
| 0.2 / 0.5 mK | 11% | 99.3 to 99.6% (heating 24 uK) |

The lane removes this entirely. It costs one shared y tilt per pass (100 us each way), which the row-select
mirror already provides.

## 5. Limits and next steps

The model is classical. Quantum ground-state effects matter when the energy gain is near one trap quantum (about
1 to 3 uK here), so the 1 to 2 uK operating points should be checked quantum-mechanically. Light-assisted processes,
photon scattering (a few per second at 0.5 mK, estimated from the far-detuned rate, not simulated) and the vector light shift are not included. The moving trap's
shape is assumed invariant along the row; least-squares synthesis keeps it within 0.1% (module 10). Its
neighbours in the same row (independent tweezers 4 um away) are not included. The static traps are ideal
Gaussian beams.

Next steps: (i) a rearrangement scheduler using lanes and passes, to count moves per atom (heating accumulates at
1 to 10 uK per move); (ii) a quantum check of the lowest-heating points; (iii) the tolerance of the static array's
own depth uniformity.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| moving-trap field, widths, frequencies | Sec. 2 | SIM | results/m18_field.json, m18_trap.npz |
| scan (hold, moving, ramp, move, path) | Sec. 1, 4 | SIM | results/m18_scan.json (1,500 atoms; review) |
| tolerances | Sec. 3 | SIM | results/m18_tol.json (review) |
| temperature rule | Sec. 1 | SIM | results/m18_temp.json |
| waiting atoms, hold margin, operating points | Sec. 1 | SIM | results/m18_bystander.json, m18_margin.json (3,000 atoms) |
| power per tweezer | Sec. 1 | DERIVED | 57 mW at 1.2 mK (m14_plan) scaled with depth |
| static waist 0.85 um, ramps, lane 2 um, 5 and 10 uK | Sec. 2 | CHOICE | stated assumptions |
| 170 um chromatic focal shift | Sec. 3 | DERIVED | 16 um/nm x 10.6 nm, diffractive objective |
