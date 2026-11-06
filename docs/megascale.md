# Rearranging 10^4, 10^5 and 10^6-atom arrays at 3-4 um pitch

Follow-up to `docs/scaling_beyond.md`, `docs/transport3d.md`, `docs/rearrangement.md` and `docs/repair.md`. Code:
`scripts/m22_megascale.py` (budget, tiling), `picsim/megasort.py` and `scripts/m23_steady.py` (assembly and steady
state), `scripts/m18_transport3d.py pitch3` (3 um pitch); `picsim/abcd.py`; results `results/m22_*.json`,
`results/m23_steady.json`, `results/m18_pitch3.json`; tests `tests/test_megascale.py`. Values are SIM, DERIVED or
CHOICE as marked, and valid in the regime stated. They should be reviewed before they are quoted.

## 1. Result

At 4 um pitch the transport side scales to 10^6 sites. One engine block of 32 rows x 4 engines (two 16-chip
stacks, 128 engines) assembles 10^6 sites from 85% loading in 1.1 s to 1e-2 defects and 1.5 s to 1e-3, with fast
handoff steps. It then holds the array at 2.3e-4 defects per site in steady state. Smaller arrays reach the same
floor faster: 10^4 sites in 0.1 to 0.3 s, 10^5 sites in 0.2 to 0.6 s. The steady-state floor is not set by the
transport. It follows the loss rates and the imaging cadence,

    rho ~ (1.5 to 1.8) x (eps_img + T_img / tau),

where eps_img is the loss per image, T_img the imaging period and tau the vacuum lifetime. With eps_img = 1e-4,
T_img = 0.1 s and tau = 3000 s this gives 2.3e-4, at every array size.

Three walls come before the transport does:

(i) **Static trapping power.** 4.5 mW per site per mK at the atoms (w0 0.85 um, 1075 nm). 10^6 sites plus
reservoir need 0.5 to 0.9 kW at the atoms at a 0.1 mK hold, and 5 to 9 kW at a 1 mK loading or imaging depth
(18 to 30 kW of laser at 30% delivery).

(ii) **The objective's field.** A single objective needs a corrected field of radius 0.2 to 0.3 mm for 10^4 sites,
0.7 to 0.9 mm for 10^5, and 2.1 to 2.8 mm for 10^6 (about 1.4e7 resolution cells). Zones of +/-250 um would need
2, 8 to 13, or 72 to 127 separate objectives.

(iii) **Pitch.** At 3 um pitch the lane between rows is 1.5 um from the waiting atoms. The moving tweezer then
pulls 55 to 99% of them out, so 3 um needs empty corridor rows, and these cost most of the density gain.

| | 10^4 sites | 10^5 sites | 10^6 sites |
|---|---|---|---|
| side at 4 / 3 um | 0.40 / 0.30 mm | 1.26 / 0.95 mm | 4.0 / 3.0 mm |
| single objective: field radius | 0.28 mm | 0.89 mm | 2.8 mm |
| zones of +/-250 um (objectives) | 2 | 13 | 127 |
| static power at the atoms, 0.1 mK hold / 1 mK (85% loading, incl. reservoir) | 5.3 / 53 W | 53 / 530 W | 0.53 / 5.3 kW |
| engine block (32 x 4): field positions | 10 | 28 | 252 |
| assembly to 1e-3 (85% loading, 32 x 4) | 1.1 s (reload-limited) | 0.12 to 0.35 s | 1.5 s (fast steps) |
| steady-state defects, tau 3000 / 300 s | 2.3e-4 / 6.0e-4 | 2.3e-4 / 6.0e-4 | 2.3e-4 (fast) / 6.0e-3 (slow, 300 s) |
| refills needed per second (tau 3000 s, eps 1e-4 at 10 Hz) | 13 | 134 | 1,335 |

## 2. Optics: one large field or many zones

A square array of side L needs a corrected field of radius L/sqrt(2). In units of the resolution cell (1.06 um at
NA 0.5) that is 372, 1,187 and 3,754 cells per side at 4 um pitch. The proposal's +/-100 um field covers none of
these. 10^4 sites need +/-0.28 mm; 10^5 sites need +/-0.9 mm, which is within the range of refractive NA 0.5
objectives (RECALLED, verify); 10^6 sites need a lithography-class field of +/-2.8 mm at NA 0.5.

Zones (independent arrays, each with its own objective, static array and engine block) avoid the large field.
They behave like the smaller arrays: a +/-250 um zone holds 7,700 sites at 4 um and assembles like the 10^4 case.
But each needs its own objective, about 25 to 30 mm across at a 15 mm working distance, so 127 zones for 10^6
are impractical. Zones suit 10^4 (2 zones) to 10^5 (13). For 10^6, a single large-field objective, or a few
large zones, is unavoidable.

## 3. Static trapping power

From `trap_depth_K` with alpha(1075 nm) = 1120 au and a 0.85 um waist: 4.5 mW per site per mK at the atoms. The
reservoir sites must be counted too. The budget uses the minimum reservoir fraction, 1 - f: 50% of all sites at
50% loading and 15% at 85%. The simulated 3:1 layout at 85% uses 25%, which adds 13% to these powers.
At the 0.1 mK hold of `docs/transport3d.md`, 10^6 sites take 0.53 kW at the atoms with 85% loading; at a 1 mK
loading or imaging depth, 5.3 kW. With 30% delivery (CHOICE) the laser must supply 1.8 to 3 kW at the hold and
18 to 30 kW at 1 mK. This is the first wall at 10^6. Three levers bear on it, all outside this study: (i) loading
and imaging at lower depth (power scales with depth); (ii) a smaller static waist; (iii) an optical lattice in an
enhancement cavity for the static potential (as in the continuously operated lattice array of Gyger et al.,
arXiv:2402.04994, search summary only).

## 4. Transport hardware: engines, tiling and blocks

One engine (192 ports, 383 fast actuators) serves a 180 um row segment: 45 sites at 4 um, 60 at 3 um. Engines can
be tiled in two directions.

**Along y, by stacking chips** (`docs/stacked_relay.md`). One chip per row, with a collimator per chip.

**Along x, several engines on one chip edge.** Each gets its own x lenslet; a common Fourier lens overlaps all of
them in the pupil and maps their offsets to adjacent atom segments. The ABCD check (`m22_tiling`) gives
magnification -60.1 for every engine and all chief rays at the pupil centre. The grating lobes are the catch:
with 1.5 um emitter pitch, each engine's first lobe lands 204 um from its spot at the atoms, inside the
neighbouring segment. A shared slit plane must therefore have opaque gaps at least as wide as the slit zones to
absorb the lobes. That sets an engine pitch of at least 4.4 mm on the chip (four engines per 18 mm edge) and lets
one block position address every other segment. Two scanner positions cover a row; baffles between the lenslet
channels are needed against stray lobe light.

**Block.** 32 rows (two polarisation-combined 16-chip stacks) x 4 engines per row = 128 engines, 49,024 fast
actuators. It holds 5,888 simultaneous tweezers, which need 85 W at 0.3 mK if all are lit (at most during
assembly). In steady state only the few refilling tweezers are lit. A 16 x 1 block (one stack) is too small for
10^6 sites: it never finishes assembling (Sec. 6).

## 5. 3 um pitch needs corridors

The 3D transport model at 3 um pitch, with a 1.5 um lane, 200 us ramps and 400 us lane steps
(`results/m18_pitch3.json`, 3,000 atoms):

| hold / moving | transported atom | waiting atom, same row, lane 1.5 um away | waiting atom, next row, beside the drop |
|---|---|---|---|
| 0.1 / 0.3 mK | 1.000 delivered | 0.335 stay, 0.584 lost | 0.149 stay, 0.851 lost |
| 0.1 / 0.5 mK | 1.000 | 0.001 stay | 0.006 stay |
| 0.2 / 0.6 mK | 1.000 | 0.225 stay | 0.117 stay |
| 0.2 / 1.0 mK | 1.000 | 0.014 stay | 0.012 stay |

The moving tweezer (1.37 um FWHM across the row, with strong out-of-focus light) captures atoms 1.5 um away. At 4 um
pitch the lane is 2 um away, and there waiting atoms are safe at a hold of about 20 kT. Transport at 3 um pitch
must therefore run in empty corridor rows. Consider a repeating layout of computation row, corridor, reservoir
row, at 3 um. Every transport path is then 3 um from occupied sites, but the computation-site density, (1/3)/9 um^-2 = 0.037 um^-2,
is only about 20% above that at 4 um with a 1:1 reservoir layout (0.031 um^-2). Moves across columns would need corridor columns too.
Hence 3 um pitch gains little with this spot size. The `m23` runs at 3 um pitch use the 4 um lane model and are
geometry-only for that reason.

## 6. Assembly and steady state

The simulation (`picsim/megasort.py`, `m23_steady`) works as follows:
- **Layout:** computation rows interleaved with reservoir rows, 1:2 at 50% loading and 3:1 at 85%.
- **Imaging:** every 0.1 s (10 ms each); 1e-4 loss per image and 1e-4 per move; the reservoir is refreshed every 1 s.
- **Visits:** the engine block goes round-robin to tiles with known holes. Each visit runs up to four repair
  passes (common row shift, order-preserving x, one move per refill).

At 4 um pitch:

| sites | loading | block | timing | to 1e-2 | to 1e-3 | steady state, tau 3000 / 300 s |
|---|---|---|---|---|---|---|
| 10^4 | 50% | 16 x 1 / 32 x 4 | slow | 0.34 / 0.12 s | same | 3.2e-4 / 5.8e-4 (16 x 1) |
| 10^4 | 85% | 32 x 4 | slow | 0.12 s | 1.11 s | 2.3e-4 / 6.0e-4 |
| 10^5 | 50% | 16 x 1 / 32 x 4 | slow | 2.0-2.7 / 0.45-0.57 s | 2.0-3.0 / 0.45-0.57 s | 2.3e-4 / 6.0e-4 |
| 10^5 | 85% | 32 x 4 | slow / fast | 0.23 / 0.12 s | 0.34 / 0.12 s | 2.3e-4 / 6.0e-4 |
| 10^6 | 50% | 32 x 4 | slow | 4.5 s (tau 3000) | not in 20 s | 1.1e-3 / 5.5e-3 |
| 10^6 | 85% | 32 x 4 | slow | 3.6 s | not in 20 s | 1.2e-3 / 6.0e-3 |
| 10^6 | 85% | 32 x 4 | **fast** | **1.1 s** | **1.5 s** | **2.3e-4** (tau 3000) |
| 10^6 | 50% / 85% | 16 x 1 | slow | never | never | 3-13% |

(Slow: 200 us ramps, 400 us lane steps, 0.2 to 0.3 uK per move. Fast: 50 us ramps, 100 us lane steps, 200 us per
40 um, about 1 to 10 uK per move from `docs/rearrangement.md`.)

Three findings:

(i) **Up to 10^5 sites, the floor is set by losses.** It is the same for every block and timing:
2.3e-4 (tau 3000 s) and 6e-4 (tau 300 s).

(ii) **At 10^6, latency matters.** With slow steps, a sweep of 252 block positions takes longer than the imaging
period. New holes then wait several sweeps, and the floor rises five-fold. Fast steps restore the
loss-limited floor; they cost heating (about 1 to 10 uK per move, but each atom moves once).

(iii) **Each refill is one move**, from a nearby reservoir row into the hole, so heating does not accumulate.

**The floor.** A loss becomes visible only at the next image, and an image itself loses eps_img per atom. Hence
rho ~ k (eps_img + T_img / tau), with k = 1.5 to 1.8 from the simulations. Imaging more often shortens the
latency but costs more imaging loss, so the eps_img term is a floor of its own:
- 1e-3 defects: eps_img + T_img/tau <= 6e-4, met by eps_img = 1e-4 and tau >= 200 s at 10 Hz imaging;
- 1e-4 defects: eps_img <= 3e-5 and tau >= 3,300 s at 10 Hz.

For comparison, Manetsch et al. (arXiv:2403.12021, abstract read 2026-11-06) report a 23 min lifetime in a
room-temperature apparatus and 99.98952% imaging survival (eps_img = 1.05e-4) for a 6,100-atom Cs tweezer array. That
places current technology at the 1e-3 to 3e-4 level.

**Supply.** In steady state at 10^6 sites the reservoir is consumed at about 1,300 atoms/s (tau 3000 s). Chiu et
al. (arXiv:2506.20660, abstract read 2026-11-06) report 30,000 initialized qubits per second from lattice conveyor
belts, maintaining more than 3,000 atoms for over 2 hours. That rate would cover the steady-state demand with a
margin of 20.

## 7. Walls in order, and what removes each

| scale | first wall | second | third | transport (this design) |
|---|---|---|---|---|
| 10^4 | field +/-0.28 mm (or 2 zones) | static power 50 W at 1 mK | losses: 2-6e-4 floor | 1 stack is enough; < 0.4 s assembly |
| 10^5 | field +/-0.9 mm or 13 zones | static power 0.5 kW at 1 mK | losses | 32 x 4 block: 0.1-0.6 s assembly |
| 10^6 | static power 5 kW at 1 mK | field +/-2.8 mm (SBP 1.4e7) | vacuum and imaging loss (rho >= 2e-4 now); supply 1.3k atoms/s | 32 x 4 block with fast steps; 1.5 s assembly |

**Recommendations.**
- **10^4:** one large-field objective or two zones, one 16-chip stack, 4 um pitch.
- **10^5:** a 32 x 4 engine block, 4 um pitch, a field of +/-0.9 mm or about 13 zones, and about 85% loading
  (it cuts the reservoir from two thirds of the sites to a quarter, and the static power by 40%).
- **10^6:** feasible on the transport side with the 32 x 4 block and fast steps. The static potential needs a
  different source (lattice or much lower depths) and a lithography-class objective, and the loss rates must
  improve to reach 1e-4 defects.

## 8. Limits and next steps

The steady-state model assumes perfect detection apart from imaging loss, an external reservoir supply, a single
objective whose field reaches every site, and one move per refill. It does not schedule the reservoir refresh
physically (conveyor timing) or include photon scattering. Heating per move for the fast steps is taken from the
1D quantum and 3D classical results (1 to 10 uK) and is not tracked per atom here. The x-tiling optics is first
order (no baffle or lenslet aberration design). The 3 um corridor layouts are argued from the 4 um transport
results, not simulated.

Next steps: (i) a corridor layout at 3 um pitch in the 3D and steady-state models; (ii) a lattice-based static
potential and its power budget at 10^6; (iii) an optical design of the large-field objective or of the zone
layout; (iv) per-atom heating in the steady-state model, with fast and slow steps mixed (slow for refills, fast
for assembly).

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| geometry, SBP, zones, static power, blocks, loss rates | Sec. 1-4 | DERIVED | results/m22_budget.json |
| x-tiling lobe geometry and ABCD | Sec. 4 | DERIVED | results/m22_tiling.json |
| 3 um pitch transport | Sec. 5 | SIM | results/m18_pitch3.json (3,000 atoms; review) |
| assembly and steady state | Sec. 1, 6 | SIM | results/m23_steady.json (review) |
| rho ~ 1.5-1.8 (eps + T_img/tau) | Sec. 1, 6 | SIM/DERIVED | fit to m23 floors |
| 23 min lifetime, 99.98952% imaging survival, 6,100 atoms | Sec. 6 | SOURCE | Manetsch et al., arXiv:2403.12021, abstract read 2026-11-06 |
| 30,000 qubits/s, > 3,000 atoms for > 2 h | Sec. 6 | SOURCE | Chiu et al., arXiv:2506.20660, abstract read 2026-11-06 |
| 2,088 sites at 6 K, 828-atom rearrangement | context | SOURCE | Pichard et al., arXiv:2405.19503, abstract read 2026-11-06 |
| continuous lattice array > 1,000 atoms | Sec. 3 | SEARCH SUMMARY | Gyger et al., arXiv:2402.04994 (abstract not read) |
| eta 0.3, hold 0.1 mK, load 1 mK, w0 0.85 um, eps_img, p_move, T_img, t_reload | Sec. 3, 6 | CHOICE | stated assumptions |
| refractive NA 0.5 objective fields | Sec. 2 | RECALLED | verify |
