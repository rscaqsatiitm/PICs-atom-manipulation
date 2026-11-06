# Rearranging 10^4, 10^5 and 10^6-atom arrays at 3-4 um pitch

> Follow-up (Sec. 9): a 3 x 4 um lattice (3 um along the rows, rows 4 um apart) transports as cleanly as 4 x 4 um
> and has 33% more sites. This supersedes the corridor conclusion of Sec. 5. A free-space lattice needs more power
> than tweezers; only a cavity lattice (buildup ~100) needs less. At 10^6 sites, fast handoff steps are needed in
> steady state as well, at a cost of at most 10 uK per refilled atom.
>
> Follow-up (Sec. 10): engines tiled along x can sit at 3.0 mm pitch (ghost below 2e-7 of a tweezer peak), so six
> fit on an 18 mm chip edge. The 3 x 4 um floor of 8e-4 in Sec. 9 was an edge artefact; it is 2.4e-4 at image, as on
> 4 x 4 um. The defect fractions quoted in Sec. 1-9 are sampled just after an image and count imaging losses
> twice; the time average is about 25% lower (1.8e-4 against 2.4e-4 at the nominal point). Intermediate steps
> (2 uK per move) reach a time average of 2.2e-4 at 10^6 sites with the 128-engine block if the imaging period is
> raised from 0.1 to 0.3 s.
>
> Follow-up (Sec. 11): stacks combine in space as well as by polarisation. A slotted 45 deg comb mirror in the
> 184 um gaps between chip beams at the collimator plane interleaves two 16-chip stacks built at 1.0 mm pitch, at a
> worst-row cost of 0.45% in power and, in 2D, a full-row crosstalk of 2.3e-4 at 4 um (1.8e-4 without the
> comb). Two such pairs and a PBS form a 64-row block,
> so the 64 x 4 block (256 engines on 12 mm chip edges) is buildable.
>
> Follow-up (Sec. 12): with the real (glued, traced) collimator the comb results hold. But a cylindrical
> collimator at NA 0.9 defocuses the light of tweezers away from the row centre (equivalent facet focus error
> -99 s_x^2 um, s_x = x_a / 288 um). With the best common focus offset a row is usable over |x_a| <= 45 um, half
> of the 180 um segment. Earlier power and crosstalk figures hold near the row centre only; a collimator that is
> insensitive to the x angle (for example all-reflective) is needed for full segments.

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
one block position address every other segment (Sec. 10 supersedes this pitch: the wave-optics check allows
3.0 mm, six engines per edge). Two scanner positions cover a row; baffles between the lenslet
channels are needed against stray lobe light.

**Block.** 32 rows (two polarisation-combined 16-chip stacks, or comb-interleaved, Sec. 11) x 4 engines per row = 128 engines, 49,024 fast
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

## 9. Follow-up: denser lattice, static potential and step timing

**3 um along the rows, 4 um across (`results/m18_aniso.json`, 3,000 atoms per case, 200 us ramps and 400 us lane
steps).** Only the closest approach to waiting atoms matters, and that is set by the lane (2 um from the rows), not
by the pitch along the row. Hence a 3 x 4 um lattice transports as cleanly as 4 x 4 um:

| layout | hold / moving | transported | heating | waiting atoms lost: in-row neighbour at 3 um, passed by the lane, next row beside the drop |
|---|---|---|---|---|
| 3 x 4 um, lane 2 um | 0.1 / 0.3 mK (5 uK) | 3000/3000 | 0.16 uK | 0 / 0 / 0 |
| 3 x 4 um, lane 2 um | 0.1 / 0.5 mK (5 uK) | 3000/3000 | 0.05 uK | 0 / 0 / 3 of 3000 |
| 3 x 4 um, lane 2 um | 0.2 / 0.6 mK (10 uK) | 3000/3000 | 0.59 uK | 0 / 0 / 0 |
| 3 um, corridor rows 6 um apart, lane 3 um | all three | 3000/3000 | 0.3-1.5 uK | 0 / 0 / 0 |

The 3 x 4 um lattice has 33% more sites per area than 4 x 4 um (0.063 against 0.047 computation sites per um^2
with a 3:1 reservoir layout), at no transport cost. The corridor layout is safe but less dense. This replaces
the conclusion of Sec. 5: 3 um pitch works along the rows. Only the row spacing must stay at 4 um.

**Static potential: tweezers or a lattice (`results/m22_lattice.json`, DERIVED).** The figures are laser power at
30% delivery, for N sites plus a 25% reservoir at 4 x 4 um, at depths of 0.1 / 1 mK. The lattice uses two
orthogonal retro-reflected beams (incoherent sum, peak 8 I per beam; coherent four-beam: 16), held within -20% of
the peak over the array.

| static potential | 10^4 sites | 10^5 sites | 10^6 sites |
|---|---|---|---|
| tweezer array (w0 0.85 um) | 20 / 200 W | 0.2 / 2.0 kW | 2.0 / 20 kW |
| free-space lattice, Gaussian beams | 0.5 / 5.0 kW | 5 / 50 kW | 50 / 500 kW |
| free-space lattice, flat-top beams (incoherent / coherent) | 85 / 42 W ; 0.85 / 0.42 kW | 0.85 / 0.42 ; 8.5 / 4.2 kW | 8.5 / 4.2 ; 85 / 42 kW |
| cavity lattice, Gaussian, buildup 10 / 100 | 50 / 5 W ; 0.5 / 0.05 kW | 0.5 / 0.05 ; 5 / 0.5 kW | 5 / 0.5 ; 50 / 5 kW |

(Each pair is 10 / 100 buildup or incoherent / coherent; then ";" separates 0.1 mK and 1 mK.) A tweezer puts its
power into about pi w0^2/2 = 1.1 um^2 per site. A Gaussian lattice beam that keeps the corners within 20% spreads
it over about 15 a^2 / G = 31 um^2. Hence a free-space lattice loses by about 25 with Gaussian beams and by 2 to 4
with flat-top beams. Only a cavity with a buildup near 100 wins, by a factor of about 4. That cavity would need a
waist near 10 mm for 10^6 sites, which is demanding. Two levers do more than the choice of potential: depth (a
factor of 10 between 1 and 0.1 mK), and a higher delivery efficiency for the simpler lattice optics (up to about
2.7 if 0.8 instead of 0.3). For 10^6 sites, the static potential needs about 2 kW at the hold depth with tweezers,
or about 0.5 kW with a buildup-100 cavity lattice. Loading and imaging must happen at a lower depth, or in a
smaller zone.

**Step timing in steady state, with heating per atom (`results/m23_mixed.json`, 32 x 4 block, 85% loading,
tau 3000 s).** Heating per move: 0.3 uK with slow steps and 10 uK with fast steps (Sec. 6).

| sites, pitch | timing | to 1e-3 | steady-state defects | atoms moved | heating: mean / max |
|---|---|---|---|---|---|
| 10^5, 4 x 4 | slow | 0.35 s | 2.5e-4 | 16% | 0.05 / 0.3 uK |
| 10^5, 3 x 4 | slow | 0.23 s | 2.5e-4 | 16% | 0.05 / 0.3 uK |
| 10^6, 4 x 4 | slow | not in 20 s | 1.2e-3 | 17% | 0.05 / 0.3 uK |
| 10^6, 4 x 4 | fast | 1.45 s | 2.5e-4 | 17% | 1.7 / 10 uK |
| 10^6, 4 x 4 | adaptive (fast until 1e-2, then slow) | 2.47 s | 1.2e-3 | 17% | 1.4 / 10 uK |
| 10^6, 3 x 4 | fast | 1.45 s | 8.1e-4 | 17% | 1.7 / 10 uK |

Every atom moves at most once (refills come straight from a reservoir row), so no atom gains more than one move's
heating. Up to 10^5 sites, the slow, cold steps reach the loss-limited floor. At 10^6 sites the refill load is
too high for slow steps: switching to them after assembly brings back the latency-limited floor of 1.2e-3. Fast
steps are therefore needed in steady state too, at a cost of at most 10 uK per refilled atom. The 3 x 4 um lattice
at 10^6 sites carries 33% more refills per block position and rises to 8e-4. A larger block (more engines) or
intermediate steps (about 2 uK per move at 200 us lane steps, from the quantum check) would restore 2.5e-4.
(Correction, Sec. 10: the 8.1e-4 on 3 x 4 um was an edge artefact of the block placement, not latency.)

## 10. Follow-up: tiling wave optics, remedies at 10^6 sites, imaging cadence

**Engines tiled along x: wave optics (`results/m24_tiling_wave.json`).** Sec. 4 required opaque gaps as wide as
the slit zones, from a first-order lobe argument, and hence an engine pitch Pe >= 4.4 mm. We now propagate the
emitter near field of one engine (D1, LinearSynth tweezer at x_a = 0, 60 and 90 um, i.e. up to the segment edge)
over f1 = 3.12 mm to its lenslet (clear aperture = fill x Pe), apply an ideal lens phase, propagate over f1 to the
shared slit plane Sigma and measure the intensity in the neighbouring engine's slit zone. Sigma is conjugate to the
atoms, so this is the ghost a neighbour segment sees. We keep only directions inside the pupil acceptance
A/(2 f1), since steeper light misses the objective; an exact ray trace through the neighbouring lenslets (all
angles that miss the engine's own lenslet) checks the same question geometrically.

| engine pitch Pe | worst ghost peak / tweezer peak (accepted; raw) | ghost power / own slit power | accepted rays through neighbouring lenslets |
|---|---|---|---|
| 8.0 mm | 4.5e-13 ; 1.7e-9 | 2.6e-13 | none |
| 4.43 mm | 3.6e-10 ; 2.5e-7 | 1.4e-10 | none |
| 3.5 mm | 1.7e-8 ; 3.1e-6 | 1.0e-8 | none |
| 3.0 mm | 1.7e-7 ; 7.7e-6 | 6.3e-8 | none |
| 2.6 mm | 1.6e-6 ; 1.3e-5 | 7.5e-7 | none |
| 2.4 mm | 6.7e-5 ; 5.6e-5 | 7.2e-5 | none |

(Worst case over x_a = 0, 60, 90 um at fill 1; it is the segment-edge tweezer, x_a = 90 um, except at 3.5 mm, where
x_a = 60 um is marginally worse. A lenslet fill of 0.8 at 4.43 mm raises the accepted ghost to 2e-9.) The first-order lobes are not the limit: the lenslet aperture and the
pupil acceptance remove them, and no accepted ray reaches another slit zone at any pitch scanned. What sets the
pitch is the lenslet aperture needed to pass the engine's own slit zone, 2 f1 tan(asin s_cut) = 2.37 mm. Below
2.6 mm the lenslet clips the own beam and its diffraction tail enters the neighbour zone. At Pe = 3.0 mm the
accepted ghost is below 2e-7 of a tweezer peak, under the coherent criterion of 1e-6 (`docs/row_coherence.md`).
Hence six engines fit on an 18 mm chip edge instead of four. Two scanner positions per row are still needed:
the segment pitch at the atoms, Pe f_obj / f_c = 0.092 Pe = 277 um at Pe = 3.0 mm (m22 combiner, f_c 188 mm,
f_obj 17.3 mm), exceeds the 180 um segment, so one block position addresses every other segment. The model is
1D (x only), uses an ideal thin lens and no surface scatter; baffles between lenslet channels remain the defence against scattered light, which this
calculation does not include.

**Intermediate steps on 3 x 4 um (`results/m18_midlane.json`, 3,000 atoms per case).** Ramps 100 us, lane steps
200 us, 300 us moves along the lane, transport over 9 um:

| hold / moving depth (atom temperature) | transported | heating per move | bystanders at 3 um, 6 um, next row (lost / displaced) |
|---|---|---|---|
| 0.1 / 0.3 mK (5 uK) | 3000/3000 | 3.0 uK | 0 / 0 / 0 |
| 0.2 / 0.6 mK (10 uK) | 3000/3000 | 0.8 uK | 0 / 0 / 0 |

Hence intermediate steps are clean on the dense lattice, at 1 to 3 uK per move. We use 2 uK per move below
(the quantum check gave a similar figure).

**Three defect fractions.** Until now we quoted the defect fraction sampled just after each image. That value
counts imaging losses twice: the image removes a fraction eps_img of the atoms after detecting them, so those holes
stay unseen for a full period, and at the next image they are still there when the new losses are added. The
simulator now also records the value just before each image and the time average (`Sim.mean_defect`, a trapezoid sum over visit, image and idle intervals; vacuum losses are uniform in time, so
the rule is exact for them in expectation, while each visit's fills are booked at its start, which biases the
average low by at most one visit time, 1 to 2 ms against T_img = 100 ms). The time average is what a
computation started at a random moment sees. With fast steps at 10^6 sites the three are 2.4e-4 (at image),
1.4e-4 (before image) and 1.8e-4 (average): the at-image value is close to 2 eps_img + (T_img + t_img)/tau =
2.4e-4, and the average is close to eps_img + (T_img + t_img)/(2 tau) = 1.2e-4 plus the repair latency.

**Remedies at 10^6 sites (`results/m23_remedy.json`, 85% loading, tau 3000 s, eps_img 1e-4, T_img 100 ms, rows
4 um apart).** Blocks of 128 to 256 engines, three step timings, heating per move 10 / 2 / 0.3 uK (fast /
intermediate / slow). Defect fraction at image / time average over the last 3 s of a 20 s run, and the time from
the first image to 1e-3 (at image):

| block (engines) | pitch along rows | fast | intermediate | slow |
|---|---|---|---|---|
| 32 x 4 (128) | 4 um | 2.5e-4 / 2.0e-4, 1.67 s | 5.3e-4 / 4.8e-4, 3.14 s | 1.2e-3 / 1.1e-3, not in 20 s |
| 32 x 8 (256) | 4 um | 2.4e-4 / 1.8e-4, 1.01 s | 3.9e-4 / 3.4e-4, 1.58 s | 8.5e-4 / 8.0e-4, 4.18 s |
| 48 x 4 (192) | 4 um | 2.3e-4 / 1.8e-4, 1.01 s | 4.1e-4 / 3.6e-4, 2.25 s | 8.7e-4 / 8.2e-4, 3.06 s |
| 64 x 4 (256) | 4 um | 2.4e-4 / 1.8e-4, 0.90 s | 3.4e-4 / 2.9e-4, 1.47 s | 6.6e-4 / 6.1e-4, 2.72 s |
| 32 x 4 (128) | 3 um | 2.4e-4 / 1.9e-4, 1.34 s | 4.3e-4 / 3.8e-4, 2.35 s | 1.1e-3 / 1.0e-3, 6.77 s |
| 32 x 8 (256) | 3 um | 2.4e-4 / 1.8e-4, 0.79 s | 3.0e-4 / 2.5e-4, 1.46 s | 6.4e-4 / 5.9e-4, 4.40 s |
| 48 x 4 (192) | 3 um | 2.4e-4 / 1.8e-4, 0.79 s | 3.7e-4 / 3.2e-4, 1.91 s | 7.5e-4 / 7.0e-4, 4.29 s |
| 64 x 4 (256) | 3 um | 2.4e-4 / 1.8e-4, 0.68 s | 3.2e-4 / 2.7e-4, 1.46 s | 5.8e-4 / 5.3e-4, 2.39 s |

Four results follow. First, the 8.1e-4 reported in Sec. 9 for 3 x 4 um was an artefact. The last row of block
positions ended short of the array edge, and the holes in rows below it had no reservoir row in any block that
reached them (563 of 716 remaining holes sat there). `Block.positions` now places the last block flush with the
edge, and 3 x 4 um holds 2.4e-4 with the 32 x 4 block and fast steps, as 4 x 4 um does. Second, fast steps reach
the loss floor with any block from 128 engines; more engines only shorten assembly (1.67 to 0.90 s at 4 um).
Third, intermediate steps (2 uK per move) come within a factor of 1.4 to 1.9 of that floor with 256 engines (time
average 2.5e-4 to 3.4e-4) but not with 128. Fourth, slow steps (0.3 uK) do not reach the floor with any block
tried; their best is 5.3e-4. The 3 um rows do slightly better than 4 um rows with the slower timings. We
attribute this to fewer block positions per sweep (a 180 um segment holds 61 sites at 3 um against 46 at 4 um, so
the 32 x 4 block has 222 positions against 252); we did not isolate this cause.

The block shape is a hardware choice. A 32-row block is two 16-chip stacks combined by polarisation (Sec. 4);
48 or 64 rows need three or four stacks and a combining method beyond polarisation (designed in Sec. 11: a
slotted comb mirror, combined with the PBS for 64 rows).
32 x 8 needs eight engines per chip edge: at the 3.0 mm pitch allowed above that is a 24 mm edge. At 2.6 mm
(21 mm edge) the accepted ghost rises to 1.6e-6, just above the 1e-6 criterion. Hence our preferred block for
10^6 sites is 32 x 8 on a 24 mm chip edge (98,048 fast actuators), with fast steps if the 10 uK per refilled atom
is acceptable and intermediate steps (time average 2.5e-4 on 3 x 4 um) if it is not.

**Imaging cadence and loss (`results/m25_imaging.json`, `results/m25_imaging_fit.json`).** We scan the imaging
period T_img (0.03, 0.1, 0.3, 1 s), the loss per image eps_img (1e-5, 1e-4, 1e-3) and the vacuum lifetime tau
(300, 3000, 30000 s) at 10^5 sites (4 x 4 um, 85% loading, 32 x 4 block, intermediate steps, t_img 10 ms, 8 s
runs, averages over the second half), add t_img = 3 and 20 ms, and take three points at 10^6 sites per timing.
Selected time-averaged defect fractions at 10^5 sites:

| eps_img, tau | T_img 0.03 s | 0.1 s | 0.3 s | 1 s |
|---|---|---|---|---|
| 1e-5, 3000 s | 2.0e-5 | 3.1e-5 | 6.7e-5 | 1.8e-4 |
| 1e-4, 300 s | 2.9e-4 | 4.1e-4 | 7.8e-4 | 1.9e-3 |
| 1e-4, 3000 s | 1.6e-4 | 1.3e-4 | 1.8e-4 | 3.1e-4 |
| 1e-4, 30000 s | 1.4e-4 | 1.2e-4 | 1.2e-4 | 1.5e-4 |
| 1e-3, 3000 s | 3.1e-3 | 1.4e-3 | 1.2e-3 | 1.3e-3 |

The values at image fit rho_image = 2.42 eps_img + 0.96 (T_img + t_img)/tau (median error 12%, worst 33%), close
to the 2 eps_img + (T_img + t_img)/tau expected from the double counting. This replaces the rule
rho ~ 1.5-1.8 (eps_img + T_img/tau) of Sec. 1 and 6, which was fitted to at-image values over a narrower range.
The time average follows a hole-lifetime model. An imaging loss is unseen for one period T = T_img + t_img and
is then repaired after a mean delay d; a vacuum loss waits T/2 on average for the next image, and then d. Hence

rho_mean = eps_img (T + d)/T + (T/2 + d)/tau.

One delay, d = 25 ms, fits all 40 points at 10^5 sites with a median error of 10% (d = 0 underestimates by a
median 16%). The worst point (61%) is eps_img = 1e-3 at T_img = 30 ms, where the repair saturates (below).
Minimising rho_mean over T gives the optimum imaging period, T* = (2 eps_img d tau)^(1/2) when d << T. For
d = 25 ms and tau = 3000 s this is 40 ms, 0.12 s and 0.39 s for eps_img = 1e-5, 1e-4 and 1e-3, in agreement with
the minima of the scan (0.03, 0.1 and 0.3 s). At the optimum, rho_mean ~ eps_img + (2 eps_img d / tau)^(1/2) + d/tau:
for eps_img = 1e-4 and tau = 3000 s that is 1.5e-4 (simulated 1.3e-4 at T_img = 0.1 s). Hence the loss per image dominates the floor
once tau exceeds a few thousand seconds. Lowering eps_img, not imaging more often, is the lever. The imaging time
matters only through the dead time at short periods: at T_img = 30 ms, t_img = 3 / 10 / 20 ms gives
1.35 / 1.55 / 1.70e-4.

At 10^6 sites the delay is not a constant, because the block needs a sweep over all its positions to reach a hole
(252 positions for a 32 x 4 block on 4 x 4 um; a one-pass visit takes 0.9 to 1.0 ms with fast steps and 1.35 to
1.5 ms with intermediate steps, from `Timing.pass_time`). Time averages for eps_img = 1e-4 and
tau = 3000 s:

| timing (32 x 4 block) | T_img 0.03 s | 0.1 s | 0.3 s |
|---|---|---|---|
| fast | 1.2e-3 | 1.8e-4 | 2.0e-4 |
| intermediate | 2.3e-3 | 4.9e-4 | 2.2e-4 |

At T_img = 30 ms both timings saturate: the loss rate, eps_img / T + 1/tau = 2.8e-3 per site per second, needs
2,800 fills/s. At rho = 2e-4 a block position holds about 0.8 holes (2e-4 x 10^6 / 252), so a sweep of 252 visits
(0.25 to 0.4 s) refills only about 200 atoms, and the imaging dead time takes 25% of the time. The holes
therefore pile up until each visit has enough work. By Little's law, the mean
hole lifetime is rho_mean / (eps_img / T + 1/tau): 0.15 s for fast steps at T_img = 0.1 s, and 0.40 and 0.34 s
for intermediate steps at 0.1 and 0.3 s. Intermediate steps at T_img = 0.3 s reach 2.2e-4, within 20% of fast
steps at 0.1 s, with the 128-engine block and at one fifth of the heating per refilled atom (2 against 10 uK).
Hence, at 10^6 sites, a slower imaging cadence is an alternative to a larger block. The imaging period should
match the sweep time of the block rather than be as short as possible.

## 11. Follow-up: combining stacks into blocks of 32 to 64 rows

**Why stacks can be combined in space.** The plane S, formed by the back focal planes of the per-chip
collimators, is conjugate to the atoms with magnification 1/D (D = 125 at 4 um rows; `docs/stacked_relay.md`,
Sec. 2). At S, each chip's beam is a stripe of +/- f_FAC NA_FAC = +/-158 um on the 0.5 mm chip pitch, so
neighbouring beams leave 184 um gaps. A slotted mirror at 45 deg (a "comb") placed in these gaps interleaves two
stacks built at twice the chip pitch. Stack A passes through the slots; stack B, turned by 90 deg, is reflected
by the teeth. Each beam then sees a window of half-width a <= 250 um. The comb uses no polarisation, and the
downstream relay is the same as for a single stack of all the rows. In addition, it doubles the chip pitch to
1.0 mm, which leaves room for the fibre and electrical access that Sec. 7 of `docs/stacked_relay.md` left open
(383 fast actuators per engine, 1,532 lines per four-engine chip). Our choice of a 45 deg comb is a compromise
between two considerations. First, it keeps the optical path of every row of stack B equal to that of stack A:
with B's facets in the plane y = -z_c, the path to any plane z = Z is Z for all rows. However, the comb then
spans the combined stack height along z, so rows at the far end meet it farther from S.

**The penalty is set by the light the objective accepts (`results/m26_edge.json`).** We propagate the field of
one chip in y from the collimator exit (aperture f_FAC before S, sine-condition lens, beta 0.7 um, NA_FAC 0.9)
to the comb at a distance z from S, apply the window, propagate back to S, and continue to the atoms through
the pupil (`picsim/opa/stack.RelayY`). Without a comb this reproduces the m13 values (transmission 0.973,
crosstalk at 4 um 1.09e-4 against 1.06e-4). Single tweezer, hard collimator edge, slot window of half-width a:

| z from S | a = 250 um: loss / peak / crosstalk at 4 um / at 3 um | a = 225 um: loss / peak / crosstalk at 4 um |
|---|---|---|
| no comb | 0 / 1 / 1.1e-4 / 2.9e-4 | 0 / 1 / 1.1e-4 |
| 10 mm | < 1e-3 / 0.998 / 1.0e-4 / 2.8e-4 | 6e-4 / 0.992 / 1.1e-4 |
| 20 mm | 2.5e-3 / 0.992 / 1.5e-4 / 4.0e-4 | 7.2e-3 / 0.985 / 1.8e-4 |
| 30 mm | 8.5e-3 / 0.987 / 2.7e-4 / 7.1e-4 | 1.6e-2 / 0.981 / 3.1e-4 |
| 45 mm | 1.9e-2 / 0.973 / 6.0e-4 / 1.7e-3 | 3.4e-2 / 0.956 / 7.9e-4 |
| 60 mm | 3.6e-2 / 0.953 / 1.3e-3 / 2.9e-3 | 5.3e-2 / 0.931 / 1.6e-3 |

The hard edge of the collimator is not the problem. Its Fresnel tail spreads into the gap, but it consists of
directions steeper than the pupil accepts, and clipping it costs nothing. A graded collimator edge
(raised cosine 0.8-0.9) therefore does not help: it raises the baseline crosstalk to 1.5e-4 and gains no
distance. (A Fresnel-tail estimate that ignores the pupil predicts 0.7% loss at 10 mm; the simulated loss there
is below 0.1%.) What matters is the light inside the pupil acceptance, which leaves S at half-angle NA/D = 4 mrad.
It reaches the window edge at

z_max = (a - f_FAC NA_FAC) D / NA,

which is 23 mm for a = 250 um and 17 mm for a = 225 um. This agrees with the scan: at z = 20 mm and a = 250 um
the worst row loses 0.25% and the crosstalk at 4 um rises from 1.1e-4 to 1.5e-4, while at 30 mm it reaches
2.7e-4. A single edge (block combining, rows 0-15 from one stack and 16-31 from the other) behaves similarly. It
also shifts the boundary row by 1.7 nm at z = 20 mm (a = 250 um) and by 12 nm at a = 200 um, a static offset to
calibrate. The rows are mutually incoherent (`docs/row_coherence.md`), so this crosstalk adds in intensity.

**2D check (`results/m26_twod.json`).** The 1D model leaves out the circular pupil, which couples x and y, and the
full rows. We therefore pass the comb-windowed y pupil of both chips through the 2D model of `docs/stacked_relay.md`
(D1 x response, circular NA 0.5 pupil). Chip m forms a full row of 41 tweezers at 3 um, and chip m+1 forms
11 tweezers offset by 1.5 um; one tweezer of chip m then moves by 0.25 um. Without a comb this reproduces the m13
values. Rows incoherent, neighbour row at 4 um (at 3 um in brackets):

| comb (both chips windowed) | single tweezer: peak / crosstalk at 4 um | full row: crosstalk at 4 um, max over x | neighbour, static: depth / position | neighbour, moving tweezer: depth / position |
|---|---|---|---|---|
| none | 1 / 1.1e-4 | 1.8e-4 | 1.7e-4 / 0.002 nm (2.0e-4) | 2.4e-5 / 0.04 nm (6.2e-5) |
| z 15 mm, a 250 um | 0.995 / 1.1e-4 | 1.8e-4 | 1.8e-4 / 0.02 nm (2.2e-4) | 2.9e-5 / 0.03 nm (7.4e-5) |
| z 23 mm, a 250 um (worst row, bold layouts) | 0.990 / 1.5e-4 | 2.3e-4 | 2.1e-4 / 0.08 nm (2.6e-4) | 3.9e-5 / 0.03 nm (9.6e-5) |
| z 23 mm, a 225 um | 0.982 / 1.6e-4 | 2.4e-4 | 2.2e-4 / 0.12 nm (2.2e-4) | 3.9e-5 / 0.03 nm (8.0e-5) |
| z 30 mm, a 250 um (beyond z_max) | 0.984 / 2.1e-4 | 2.9e-4 | 2.5e-4 / 0.14 nm (3.1e-4) | 5.3e-5 / 0.03 nm (1.2e-4) |
| single edge, z 23 mm, a 250 um (block boundary) | 0.995 / 1.3e-4 | 2.0e-4 | 1.9e-4 / 0.07 nm (2.2e-4) | 3.4e-5 / 0.03 nm (8.4e-5) |

Hence 2D does not change the conclusions. At the worst row the full-row crosstalk rises by a factor 1.3, and the
effect of a moving tweezer on the next row by a factor 1.55 to 1.65 (to 3.9e-5 in depth at 4 um, 9.6e-5 at 3 um). This is
still more than an order of magnitude below the in-row figures of D1 (0.11%, 3.5 nm). Tweezers at x = 0 and at the
row end (60 um) give the same values to 0.1%, because the comb acts on y alone and the circular pupil couples the
two only weakly. The 1D single-tweezer figure (1.8e-4 at z = 23 mm) slightly overstated the 2D one (1.5e-4). The
block-combining edge moves the boundary row by 7 nm in 2D (peak position), a static offset to calibrate. With
coherent rows the comb would raise the static depth change at 4 um from 2.6% to 2.9%; rows must stay incoherent,
as they must without a comb.

**Comb roll and edge straightness.** A comb edge that is not parallel to x makes the window depend on x. Beyond the
slit plane the x relay is telecentric (chief rays leave the slit plane parallel to the axis), so each tweezer
meets the comb at its slit-plane position X = f1 s. Four engines at 3.0 mm pitch span X = +/-5.6 mm. One tweezer's
footprint at the comb is at most +/-0.78 mm (17 mm beyond the slit plane, half-angle (A/2)/f1 = 0.046). Hence, for a
roll of 2 mrad, the edge moves by only 3 um across one tweezer, against the tens of um over which the result
changes. We therefore treat each tweezer with its local window a(X), without a non-separable propagation. In 2D,
at z = 23 mm, the peak falls from 0.989 to 0.972 as a goes from 250 to 205 um (about 0.04% per um), while the
crosstalk at 4 um stays at 1.4 to 1.6e-4. A roll therefore mainly makes the depth vary along the row: about +/-0.4%
for a +/-10 um change of a, which a 1.8 mrad roll gives over +/-5.6 mm. The engines set each tweezer's power in any
case, so this is a calibration, not a limit. The same holds for edge straightness of +/-10 um over 12 mm.

**Layouts (`results/m26_layout.json`; first comb point z_c = 7 mm, just beyond the x slit plane at 6.24 mm,
CHOICE).** The worst row is the one at the far end of the comb:

| rows | layout | chips per stack, pitch | comb z range | worst row (a = 250 um): loss / peak / crosstalk at 4 um | polarisation-insensitive objective |
|---|---|---|---|---|---|
| 32 | one stack | 32, 0.5 mm | | (no comb) | no |
| 32 | 2 stacks, PBS | 16, 0.5 mm | | (no comb) | yes |
| 32 | **2 stacks, comb** | **16, 1.0 mm** | **7-23 mm** | **0.45% / 0.991 / 1.8e-4** | **no** |
| 48 | 4 stacks, comb + PBS | 12, 1.0 mm | 7-19 mm | 0.21% / 0.992 / 1.4e-4 | yes |
| 48 | 2 stacks, comb | 24, 1.0 mm | 7-31 mm | 0.93% / 0.987 / 2.9e-4 | no |
| 64 | **4 stacks, comb + PBS** | **16, 1.0 mm** | **7-23 mm** | **0.45% / 0.991 / 1.8e-4** | **yes** |
| 64 | 2 stacks, comb | 32, 1.0 mm | 7-39 mm | 1.3% / 0.984 / 3.9e-4 | no |
| 64 | 2 stacks, PBS | 32, 0.5 mm | | (no comb) | yes |

With a = 225 um (about what a 50 um tooth plate, t / (2 sqrt 2) = 18 um per edge, and +/-10 um edge
placement leave) the worst rows of the
two bold layouts lose 0.86% and reach 2.0e-4. A second comb level (stacks at 2.0 mm, combined pairwise and then
again) would put the second comb beyond z_max, so we do not use it. The single taller stack (32 or 64 chips at
0.5 mm) needs no combiner. In first-order optics it only enlarges the y-relay field: +/-8 or +/-16 mm at the
first telescope lens, and 3.7 or 7.4 mrad at the pupil. Its costs are mechanical: 16 to 32 mm of facets flush to
about 100 um, and all the I/O at 0.5 mm pitch.

**Tolerances (DERIVED from the m13 conjugations).** A comb tilt theta about x deflects stack B by 2 theta at S.
This shifts its facet image in the pupil by 2 theta F_y (F_y = 2.17 m). Matching the shift to the m13 decentre
results gives 1% peak loss at theta = 0.16 mrad and 3% at 0.28 mrad. The same tilt moves B's rows by
2 theta z / D: 0.32 nm per urad at z = 20 mm, and 50 nm at 0.16 mrad. This static offset alternates between rows
(A and B rows interleave), so it must be calibrated once against the static array. Keeping the alternating
pitch error below 3 nm needs a tilt stability of about 10 urad (drift, not alignment). Moving the comb by dz
along z changes B's path by dz, which tolerates 100 um (m13 chip piston, peak >= 0.992). It also shifts the
windows by dz in y, which should stay within about +/-10 um (the a scan). Light that misses a tooth leaves at
90 deg to a dump. Light from stack A that hits the back of a tooth is reflected towards stack B and needs an
absorbing back face.

Hence the blocks of Sec. 10 can be built as follows. (i) 32 rows: two 16-chip stacks at 1.0 mm pitch and one
comb. This needs no polarisation optics, and it gives alternate rows from different stacks, so the passive
orthogonal-polarisation decorrelation of `docs/row_coherence.md` comes for free if the objective accepts it.
(ii) 64 rows: two such comb pairs, combined by a PBS. Four 16-chip stacks with four engines per chip at 3.0 mm
pitch (12 mm edge) then form the 64 x 4 block: 256 engines, 98,048 fast actuators. With intermediate steps at
10^6 sites it reaches a time average of 2.9e-4 on 4 x 4 um and 2.7e-4 on 3 x 4 um (Sec. 10), without the 24 mm
chip edge that 32 x 8 needs. The worst rows lose 0.9 to 1.5% of their peak depth to the comb (a = 250 / 225 um),
which is recovered by the per-tweezer power the engines set in any case. The edge scan is 1D in y and the 2D check uses the
same windowed pupil; both assume an ideal thin collimator and the window of each slot at one distance z; the 45 deg tilt places the two edges of a 0.5 mm slot
0.25 mm before and after that distance, which changes z by at most 4% at z_c = 7 mm. Tooth-edge roughness,
scatter from the comb, and a lens design of the larger y-relay field are not included.

## 12. Follow-up: the real collimator

Sec. 11 used an ideal thin collimator in air. The adopted collimator (`docs/collimator.md`) is a fused-silica
Weierstrass block glued to the facet, followed by a 10 um gap and a 250 um aplanat of index 1.8. We repeat the
comb study with this traced lens (`scripts/m26_stack_combine.py lens`, `results/m26_lens.json`). We then ask how
the cylindrical collimator treats light that leaves the chip at an angle along x (`skew`, `skew_offset`;
`results/m26_skew.json`, `results/m26_skew_offset.json`). The first changes little. The second limits the whole
stacked relay.

**The comb with the real collimator.** The exit field is built from the traced lens map (exit height, wavefront
and amplitude per ray) at a plane 60 um beyond the last vertex. S lies at the traced back focal plane,
0.44 mm from the facet. Without a comb this reproduces the m14 wave model: crosstalk 8.0e-5 at 4 um against 8.0e-5,
and 3.8e-4 against 4.0e-4 at 3 um. The real geometry also moves the comb. The traced x lens of m14 puts the x slit
plane 7.37 mm from the facet, the collimator glass adds 0.15 mm in x (t (1 - 1/n) per section), and we keep
0.5 mm of clearance. The comb therefore starts 7.6 mm from S, and the far row of a 32-row comb sits at 23.6 mm.
Slot window a = 250 um (225 um in brackets):

| z from S | real collimator: loss / peak / crosstalk at 4 um | thin collimator (Sec. 11) |
|---|---|---|
| 15 mm | 5e-4 / 0.994 / 7.6e-5 (3.4e-3 / 0.986 / 7.5e-5) | 5e-4 / 0.994 / 1.2e-4 |
| 20 mm | 2.6e-3 / 0.991 / 9.2e-5 (7.6e-3 / 0.985 / 9.7e-5) | 2.5e-3 / 0.992 / 1.5e-4 |
| 23.6 mm (far row) | 5.2e-3 / 0.991 / 1.1e-4 (9.6e-3 / 0.985 / 1.2e-4) | 4.7e-3 / 0.991 / 1.8e-4 |
| 30 mm | 9.0e-3 / 0.988 / 1.9e-4 (1.7e-2 / 0.983 / 2.0e-4) | 8.5e-3 / 0.987 / 2.7e-4 |

The immersed facet puts more light near the edge of the beam at S, but loss and peak follow the thin-lens values
within 12%, and the crosstalk is lower. The rule z_max = (a - f_FAC NA_FAC) D / NA therefore still holds. In 2D at
the far row (z = 23.6 mm, a = 250 um), the full-row crosstalk at 4 um is 1.6e-4 (1.75e-4 without the comb). A
moving tweezer changes the next row by 2.8e-5 in depth at 4 um (1.2e-5 without the comb) and by 7.0e-5 at 3 um
(3.6e-5). Hence the comb conclusions of Sec. 11 hold for the real collimator, at the slightly larger distances of
the real geometry.

**Light that leaves the chip at an angle along x.** A tweezer at x_a along the row is formed by light that leaves
the emitters with x direction sine s_x = x_a / x_scale (x_scale = 288 um for D1). That is up to 0.31 at the end
of a 180 um segment and 0.355 at the slit edge. The collimator is cylindrical. Its surface normals have no x
component, so k_x is conserved, and the y-z problem is exactly the 2D trace with every index n replaced by
sqrt(n^2 - s_x^2). A unit test checks this against 3D vector refraction. All the earlier y models (m13, m14, and
Sec. 11 above) used s_x = 0. For s_x > 0 the effective index ratios of the lens rise: 1.45 -> 1.50 for the block
and 1.80 -> 1.89 for the aplanat at s_x = 0.355. The lens is then no longer aplanatic for the facet. Peak relative
to an ideal collimator at the same s_x (1D, 2D in brackets), equivalent facet focus error, and crosstalk at 4 um:

| x_a (s_x) | peak | equivalent facet defocus | residual after refocusing the facet | crosstalk at 4 um |
|---|---|---|---|---|
| 0 | 1.000 (1.000) | 0 | 0 | 8.4e-5 |
| 15 um (0.052) | 0.999 (0.997) | -0.27 um | 0.001 waves | 2.2e-4 |
| 30 um (0.104) | 0.971 (0.942) | -1.10 um | 0.005 waves | 3.0e-3 |
| 45 um (0.156) | 0.794 (0.725) | -2.46 um | 0.014 waves | 8.2e-3 |
| 60 um (0.208) | 0.481 (0.440) | -4.36 um | 0.029 waves | 6.9e-3 |
| 90 um (0.312) | 0.231 | -9.67 um | 0.096 waves | 9.0e-3 |

The error is almost purely a focus shift of the facet, dz = -99 s_x^2 um. A thin glass lens in air would give
-f (n + 1) s_x^2 / (2 n) = -148 s_x^2 um, the same order. Refocusing the facet for one tweezer restores its peak
fully; the residual after refocusing stays below 0.1 waves up to 90 um. But all tweezers of a row share the facet
and the collimator. The phase sits at S, which is conjugate to the atoms, so refocusing the atom plane does not
help either. In 2D the best atom-plane focus recovers the peak up to x_a = 30 um (1.02, with a 1.5 um y-focus
shift, outside the +/-0.5 um axial registration of `docs/transport3d.md`), and not beyond. Beyond about 60 um the
2D spot breaks up (peak 0.11 at 75 um), because the y pupil is then filled to its rim and the circular pupil cuts
it differently for each x pupil column.

**What a common focus offset recovers.** A single facet offset for the whole row (the length of the glued block,
set at fabrication) trades the centre against the ends. With dz = -1.25 um the peak stays at or above 0.96 out to
|x_a| = 45 um, and the static crosstalk at 4 um stays at or below 3.6e-3. That is within the criteria of
`docs/row_coherence.md` (others' depth change below 3e-3 when a tweezer moves), since for incoherent rows the
moving-tweezer effect was about 0.13 of the static crosstalk in m13. We did not recompute it for these cases.
Beyond 50 um no offset helps: 0.90 at 50 um and 0.66 at 60 um. No offset brings the worst peak over
|x_a| <= 90 um above 0.43.

Hence a refractive cylindrical collimator at NA 0.9 serves about half of the 180 um segment (|x_a| <= 45 um,
about 23 sites at 4 um), not all of it. Every power and crosstalk figure of m13, m14 and Sec. 11 holds near the
row centre and is optimistic for the outer half. For the large arrays this doubles the number of engines per row
at fixed coverage: 90 um segments instead of 180 um. Alternatively, the collimator must become insensitive to s_x.
Serving |x_a| <= 90 um within the +/-0.5 um focus tolerance, with the best common offset, needs
|dz/s_x^2| <= 10 um: ten times less than now. We see four routes, none of them designed here.

(i) An all-reflective collimator in one medium. The law of reflection and a uniform medium are unchanged by the
effective-index map, so W = 0 for every s_x; a flat exit face normal to the collimated beam keeps it so (DERIVED).
It does not obey the sine condition and needs a folded geometry within the 0.5 mm chip pitch.

(ii) A cemented refractive design. Glass-glass interfaces are less sensitive than glass-air ones (relative power
change 0.41 s_x^2 for 1.45 | 2.0 against 0.84 s_x^2 for 1.45 | air, DERIVED), but a factor of two is not enough on
its own.

(iii) A corrector at the x slit plane, where x encodes s_x. This fails at large s_x. At x_a = 90 um the error acts at
S like a y lens of 4 mm focal length (2.9 waves at 158 um), so the beam has passed a focus before the slit plane,
7 mm later.

(iv) A larger vertical facet mode (lower facet NA). This helps by less than a factor of two, because the optimal
pupil fill, and with it f_FAC, grows with the mode.

The downstream cylindrical lenses of the y relay see x angles below about 0.05 and are not modelled here.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| geometry, SBP, zones, static power, blocks, loss rates | Sec. 1-4 | DERIVED | results/m22_budget.json |
| x-tiling lobe geometry and ABCD | Sec. 4 | DERIVED | results/m22_tiling.json |
| 3 um pitch transport | Sec. 5 | SIM | results/m18_pitch3.json (3,000 atoms; review) |
| assembly and steady state | Sec. 1, 6 | SIM | results/m23_steady.json (review) |
| rho ~ 1.5-1.8 (eps + T_img/tau) | Sec. 1, 6 | SIM/DERIVED | fit to m23 floors (at-image values; superseded by Sec. 10) |
| 23 min lifetime, 99.98952% imaging survival, 6,100 atoms | Sec. 6 | SOURCE | Manetsch et al., arXiv:2403.12021, abstract read 2026-11-06 |
| 30,000 qubits/s, > 3,000 atoms for > 2 h | Sec. 6 | SOURCE | Chiu et al., arXiv:2506.20660, abstract read 2026-11-06 |
| 2,088 sites at 6 K, 828-atom rearrangement | context | SOURCE | Pichard et al., arXiv:2405.19503, abstract read 2026-11-06 |
| continuous lattice array > 1,000 atoms | Sec. 3 | SEARCH SUMMARY | Gyger et al., arXiv:2402.04994 (abstract not read) |
| eta 0.3, hold 0.1 mK, load 1 mK, w0 0.85 um, eps_img, p_move, T_img, t_reload | Sec. 3, 6 | CHOICE | stated assumptions |
| refractive NA 0.5 objective fields | Sec. 2 | RECALLED | verify |
| anisotropic and corridor transport | Sec. 9 | SIM | results/m18_aniso.json (3,000 atoms; review) |
| tweezer and lattice power | Sec. 9 | DERIVED | results/m22_lattice.json (interference gain 8 / 16, -20% uniformity, eta 0.3, cavity buildup CHOICE) |
| mixed timing and heating | Sec. 9 | SIM | results/m23_mixed.json (review) |
| ghost in the neighbour slit zone vs engine pitch | Sec. 10 | SIM | results/m24_tiling_wave.json (1D, ideal thin lenslet, no surface scatter; review) |
| segment pitch at the atoms 0.092 Pe | Sec. 10 | DERIVED | f_obj / f_c from results/m22_tiling.json |
| intermediate-step transport on 3 x 4 um | Sec. 10 | SIM | results/m18_midlane.json (3,000 atoms; review) |
| remedies at 10^6 sites, edge-artefact correction | Sec. 10 | SIM | results/m23_remedy.json (one seed per case; review) |
| imaging map, rho_image fit, latency d = 25 ms | Sec. 10 | SIM | results/m25_imaging.json, results/m25_imaging_fit.json (one seed per case; review) |
| T* = (2 eps d tau)^(1/2), Little's-law hole lifetimes, holes per position | Sec. 10 | DERIVED | from the hole-lifetime model and the SIM values |
| coherent ghost criterion 1e-6 | Sec. 10 | DERIVED | docs/row_coherence.md |
| 24 mm chip edge for 32 x 8 | Sec. 10 | open | not designed |
| comb edge penalty vs distance and window | Sec. 11 | SIM | results/m26_edge.json (1D, thin collimator, both slot edges at one z; review) |
| z_max = (a - f_FAC NA_FAC) D / NA | Sec. 11 | DERIVED | in-band half-angle NA/D at S |
| layouts, worst-row penalty, y-relay field | Sec. 11 | SIM/DERIVED | results/m26_layout.json; scripts/m13_stack.design |
| comb tilt tolerance, row shift per urad | Sec. 11 | DERIVED | m13 decentre tolerance (results/m13_tolerance.json) and conjugations |
| z_c = 7 mm, tooth plate 50 um, +/-10 um edge placement | Sec. 11 | CHOICE | stated assumptions |
| 2D comb check: full row, neighbour static and moving | Sec. 11 | SIM | results/m26_twod.json (m13 Focus2D; review) |
| roll and straightness by a local window a(X) | Sec. 11 | DERIVED | telecentric slit space; footprint (A/2)/f1 |
| comb fabrication, edge quality, scatter; mount stability ~10 urad | Sec. 11 | open | [REF?] |
| comb with the traced collimator, real slit-plane distance | Sec. 12 | SIM | results/m26_lens.json (1D and 2D; review) |
| skew-ray detuning of the collimator vs x_a | Sec. 12 | SIM | results/m26_skew.json (effective-index trace, exact for cylindrical surfaces; downstream relay ideal; review) |
| common facet offset | Sec. 12 | SIM | results/m26_skew_offset.json (first-order focus phase, 1D; review) |
| thin-lens skew estimate, glass-glass sensitivity, corrector focal length | Sec. 12 | DERIVED | formulas in the text |
| moving-tweezer / static crosstalk ratio 0.13 | Sec. 12 | SIM | m13 (results/m13_crosstalk.json); not recomputed for skewed tweezers |
| skew-insensitive collimator (reflective) | Sec. 12 | open | not designed |
