# Scaling beyond one row of about 40 tweezers

> Follow-up (`docs/megascale.md`): scaling to 1e4, 1e5 and 1e6 sites. A 32 x 4 engine block assembles 1e6 sites
> in 1.5 s and holds them at the loss-limited floor rho ~ 1.5-1.8 (eps_img + T_img/tau); static trapping power,
> the objective's field and the 3 um lane distance are the walls.

> Follow-up (`docs/transport3d.md`): the 3D simulation with handoff refines lever 1. The moving trap must be at
> least three times the static hold depth, the hold about 20 kT, and moves must run along lanes between the rows.
> The minimum moving depth is therefore 0.3 mK for 5 uK atoms (0.6 mK for 10 uK), not 0.2 mK. That gives 10 to 21 W
> for 720 tweezers (17 to 34 W for 1 to 2 uK heating per move).

Follow-up to `docs/scaling.md`, `docs/stacked_relay.md`, `docs/collimator.md` and `docs/row_coherence.md`. Code:
`scripts/m17_scale.py`; results `results/m17_moves.json` (+ .png), `results/m17_table.json`. Transport results
are SIM from a 1D classical model; the scaling table is DERIVED. They should be reviewed before they are quoted.

## 1. Result

The 41 tweezers per row of design D1 are at 3 um spacing. On the 4 um site lattice of the proposal, only 30 are
useful over the same +/-60 um. In general the useful number per row is the row length divided by the site
pitch, so it cannot be raised by packing tweezers more densely. Hence scaling means longer rows and more rows.

The per-tweezer costs do not fall with scale: 8.7 fast actuators and about 57 mW at 1.2 mK per simultaneous
tweezer. Two levers nevertheless remove most of the penalty.

(i) **Shallower transport.** Moves of three sites (12 um) in 200 us survive at 0.2 mK with the 200 kHz update
rate (loss < 1e-3, energy gain < one trap quantum). Moves of ten sites in 200 us survive at 0.5 mK. This
cuts the laser power 2.4 to 6 times: 720 simultaneous tweezers need 6.9 to 17 W instead of 41 W.

(ii) **Time-sharing of rows (row select).** A mirror at a plane conjugate to the pupil in the y path tilts all
rows of a stack rigidly in y. A 16-chip stack (720 tweezers) then serves all 50 rows of a +/-100 um field, 2,250
sites, in four passes. The same tilt moves every row of the stack together in y, so the stack can carry atoms
in 2D: independent along x within each row, common along y.

| configuration | simultaneous tweezers | sites served | passes | chips | fast actuators | power at 1.2 mK / 0.2 mK |
|---|---|---|---|---|---|---|
| one D1 chip (128 ports) | 30 | 30 | 1 | 1 | 255 | 1.7 / 0.3 W |
| one 192-port chip (+/-90 um row) | 45 | 45 | 1 | 1 | 383 | 2.6 / 0.4 W |
| stack of 16 | 720 | 720 | 1 | 16 | 6,128 | 41 / 6.9 W |
| **stack of 16 + row select** | **720** | **2,250** | **4** | **16** | **6,128** | **41 / 6.9 W** |
| two stacks (orthogonal polarisation) + row select | 1,440 | 2,250 | 2 | 32 | 12,256 | 83 / 14 W |

## 2. What limits the number of tweezers

**Per row.** One input port per resolution cell (1.06 um at NA 0.5), so a row of L um needs L/1.06 ports and
holds L/4 + 1 tweezers on 4 um sites. That is 3.8 ports and 7.5 to 8.7 fast actuators per tweezer, independent
of L. The objective's corrected field sets L: the proposal's +/-100 um gives 192 ports and 45 tweezers per row.
Beyond that, the slab grows (5 mm at 512 ports) and the aberration grows (0.13 rad), but neither is a wall below
about 512 ports (`docs/scaling.md`). The wall is the objective.

**Rows.** One chip per row, each with a collimator aligned to 0.1 um (`docs/collimator.md`). Sixteen chips make
an 8 mm stack and 64 um of rows, against a 200 um field.

**Power.** Linear in the number of simultaneous tweezers and in the depth.

## 3. Lever 1: transport depth

The depth of 1.2 mK in the proposal is what fast long moves need. A 40 um move in 100 us survives only at
1.2 mK even with a smooth trajectory. With a 200 kHz staircase it does not survive at any simulated depth,
because each update step (up to 3.75 um) exceeds the 1.43 um waist. Shorter or slower moves need much less depth.
Minimum depth for loss <= 1e-3 and energy gain <= one trap quantum per move, 200 kHz staircase, atoms at 10 uK:

| move | 100 us | 200 us | 400 us |
|---|---|---|---|
| 12 um (3 sites) | 0.8 mK | **0.2 mK** | 0.1 mK |
| 40 um (10 sites) | none <= 1.2 mK | 0.5 mK | 0.2 mK |

(1D classical Monte Carlo, 4,000 atoms, Gaussian potential along the move with w = 1.43 um, minimum-jerk
trajectory.) The excitation is not monotonic in depth: at some depths the trap frequency aliases with the
update steps, as in module 3. A schedule would therefore use the shallowest depth that clears all the
resonances, around 0.2 to 0.3 mK, for moves of a few sites lasting 200 to 400 us. Rearrangement mostly needs
short moves, so a cycle of four passes at about 0.6 ms each (pick-up and release ramps assumed 100 us each,
200 us move, 200 us for the y tilt to settle) takes 2 to 3 ms. Depth can be raised for the occasional long move.

## 4. Lever 2: row select

The pupil of the y relay is conjugate to every chip's facet, and the atom plane to the collimators' focal
plane S (`docs/stacked_relay.md`). A tilt at any plane conjugate to the pupil therefore translates all rows
rigidly, without aberration in the paraxial limit. Shifting the stack's 16 rows by 64 um is a 3.7 mrad tilt at
the objective pupil. At the y telescope's internal Fourier plane, where the beam is 21.6 times smaller, it is
80 mrad optical, i.e. a 40 mrad mechanical rotation. Placing rows to 10 nm needs about 6.5 urad of mechanical
repeatability. That is galvanometer-scanner class (RECALLED; verify against a specific scanner). The mirror is
narrow in y (about 1 mm) and as wide as the x beam (about 20 mm). Translating the whole stack by 8 mm with 1 um
accuracy is the slower alternative (1 um becomes 8 nm at the atoms).

Atoms wait in the static 1075 nm array between passes. The cost is time (four passes for 50 rows) and one
precision scanner. Hardware per site drops by a factor 3.1.

## 5. Other options considered

**A second stack, polarisation-combined.** A polarising beam splitter before the objective combines a second
16-chip stack, orthogonally polarised and offset by 64 um, so the stacks never overlap. This doubles the
simultaneous capacity (32 rows, two passes). It costs a second set of hardware, so the penalty per tweezer is
unchanged. It needs a polarisation-insensitive objective (as in `docs/row_coherence.md`, Sec. 4).

**A larger corrected field.** This is the only way to longer rows. The +/-100 um field is a meta-objective target
in the proposal. Refractive plan objectives of NA 0.5 have flat fields of several hundred micrometres
(RECALLED; verify). With a +/-250 um field, 480 ports would give 110 tweezers per row, at the slab and aberration
costs of `docs/scaling.md` (R ~ 5 mm, about 0.12 rad).

**Planar tiling of surface-emitting chips (speculative, not simulated).** Chips side by side on a common carrier,
each with long grating emitters, imaged in y onto the atoms with demagnification D = tile pitch / 4 um. For a
6 mm tile pitch, D = 1500. The object-side NA in y is then 3.3e-4, which a 2 mm weak grating radiates directly, so no
collimator is needed. Chip placement errors shrink by D (10 um becomes 7 nm). The 0.1 um decentre tolerance is
replaced by an emission-angle tolerance of about 30 urad. That corresponds to delta n_eff ~ 3e-5, trimmable with
a heater (about 1.5 K) or by tilting the chip. The open questions are grating efficiency (up and down emission)
and fitting the engine (slab of about 2 mm) and a 2 mm grating within the tile. This is worth a separate study,
because it removes the stack and the collimators altogether.

**What does not work.** Packing tweezers below the site pitch adds nothing on a 4 um lattice. Several emitter
rows per chip cannot be separated at the atoms: rows a few um apart on the facet land within nanometres of each
other, and a per-row collimator cannot be fitted. Spectral multiplexing of rows through the dispersive star
coupler is excluded by the objective's chromatic focal shift (16 um/nm).

## 6. Recommendation and next steps

The low-penalty path is: 192-port chips (45 per row) filling the +/-100 um field; transport at 0.2 to 0.5 mK with
moves of 200 to 400 us; a 16-chip stack with row select serving 2,250 sites. This gives 720 simultaneous
independent tweezers on 6.9 to 17 W of laser power and 6,128 fast actuators. A second, polarisation-combined
stack doubles the parallelism if the rearrangement time matters more than hardware.

Next steps: (i) a 3D transport simulation (radial and axial motion, handoff to and from the static array) at
0.2 to 0.5 mK; (ii) a rearrangement scheduler that uses independent x moves per row and the common y tilt, to
count passes and time for realistic loading; (iii) an optical design of the row-select scanner plane; (iv) a first
study of the planar surface-emitting tiling.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| minimum depths, loss, excitation | Sec. 1, 3 | SIM | results/m17_moves.json (1D classical; review) |
| configuration table | Sec. 1 | DERIVED | results/m17_table.json |
| 57 mW per tweezer at 1.2 mK | Sec. 1, 2 | DERIVED | results/m14_plan.json |
| ports per row, slab, aberration | Sec. 2, 5 | SIM/DERIVED | results/m12_*.json |
| tilt angles and magnification 21.6 | Sec. 4 | DERIVED | conjugations of docs/stacked_relay.md; f_a = 100 mm, f_FAC = 175 um, beta = 0.7 um |
| atom temperature 10 uK, ramps 100 us, criteria | Sec. 3 | CHOICE | stated assumptions |
| galvanometer repeatability, objective fields | Sec. 4, 5 | RECALLED | verify |
| planar tiling numbers | Sec. 5 | DERIVED | order of magnitude, not simulated |
