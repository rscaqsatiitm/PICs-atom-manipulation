# Independently steerable tweezers: an integrated equivalent of the AOD

Design study for Objective 3 (proposal ANRF/ARG/2026/009937, draft v97). Code: `scripts/m9_independent.py`,
`picsim/opa/star.py`; results `results/m9_*.json`. Simulated values are SIM outputs, valid in the regime
simulated, and must be reviewed before they enter the proposal.

> Follow-up: the star coupler has since been designed by beam propagation and optimised
> (`docs/fpr_design.md`). The optimised design (192 arms, Rowland-type port arc, 1.2 um emitter facets,
> least-squares port synthesis) needs about 43 mW per 1.2 mK tweezer (0.86 W for 20) and keeps the
> independence figures below. The numbers in this note are for the first, paraxial design.

## 1. Result

A sliding-window focal-plane OPA moves 20 tweezers independently along a row, switches each on and off on
its own, and does so at one optical frequency. A dense row of 128 input ports, one per resolvable spot, feeds
a star coupler, 128 arms and 128 emitters; port k forms a spot at x_k = k lambda/(2 NA) = 1.064 um at the
atoms. A tweezer is a Gaussian-weighted group of about five ports, and it moves by sliding the window. In this
(position) basis different tweezers use different ports, so their fields add linearly and do not interfere
except through their tails. With 20 tweezers at random positions over +/-60 um (window width sigma = 0.9 ports,
calibrated):

| quantity | value |
|---|---|
| moved tweezer, 0.1 to 4 um moves: position error / peak change | <= 0.9 nm / <= 0.10% |
| the other 19 during the move: position / peak change | <= 1.5 nm / <= 0.19% |
| tweezer switched off: residual at its site / change of the others | 8.9e-6 / <= 0.08% |
| static peak spread of the 20 (depth equalised by the tree) | 1.4% |
| two tweezers in phase quadrature at 2.5 um: pull / peak change | 3.5 nm / 0.16% |
| spot FWHM along the row | 1.66 um |
| fibre power per tweezer at 1.2 mK: centre / field mean / field edge | 48 / 64 / 102 mW |
| 20 tweezers spread over the field | 1.3 W |

The cost is in actuators and their accuracy: 127 fast tunable splitters and 128 port phase shifters per chip,
and 10 nm position noise needs per-port amplitude and phase errors of about 1% and 10 mrad rms (Sec. 4).

## 2. What makes AOD tones independent

Two physical ingredients, and the integrated device must supply an equivalent of at least one (standard AOD
physics; [REF?] for a citable source):

(i) frequency tagging: each RF tone shifts its diffracted light by its own frequency, so tweezers from
different tones add incoherently on average, and each tone's amplitude sets one tweezer's depth;
(ii) a time-delay medium: the acoustic transit time across the aperture converts each tone's frequency into
its own angle; the number of resolvable spots is the RF bandwidth times the transit time.

The independence is per axis. A single AOD moves its spots independently along one line; crossed AODs give a
grid of row and column tones, not independent 2D positions [17]. Hence the fair target here is a 1D AOD along
a row. Our earlier shared-arm OPA (`docs/twenty_tweezers.md`) has neither ingredient: every tweezer uses
every arm, so phase-only synthesis couples them all.

## 3. Routes considered

| route | how independence arises | efficiency for 20 | show-stopper or cost |
|---|---|---|---|
| E sliding-window focal-plane OPA | tweezers on separate ports (position basis), linear superposition | 0.40 at centre, 0.19 at field edge (as the star coupler of m8) | 127 fast splitters + 128 phases; ~1% per-port accuracy |
| A' complex synthesis on the arms | linear superposition of 20 phase ramps (amplitude and phase per arm) | 0.75 (peak-to-mean 1.34, phases optimised) | 128 amplitude and phase modulators; moving one tweezer changes the optimum, so power headroom is needed |
| F AOD emulation on the arms (one RF tone per tweezer) | frequency tags + per-tone phase ramp | 0.44 (peak-to-mean 2.28 over arms and time) | 256 RF-driven modulators; modulator nonlinearity makes intermodulation ghosts, as in an AOD |
| D wavelength-addressed dispersive OPA (one laser per tweezer) | optical frequency tags + delay lines | high (passive chip) | a diffractive meta-objective shifts focus by 14 um/nm, so only +/-71 pm (+/-19 GHz) is allowed; even 5 GHz per site then needs 6.7 m of delay lines to cover the field |
| G integrated surface-acoustic-wave deflector in the slab | a literal AOD | unknown (photoelastic response of SiN not established) | 128 spots need 640 MHz of RF bandwidth on a 1 mm aperture, 128 MHz on 5 mm (v_SAW ~5 km/s, recalled) |
| T time-multiplexed OPA (one tweezer at a time) | time-averaged potentials | 1 per slot | 17 MHz pattern rate, 60 ns slots, 20x peak intensity |

Route E is the only one that is single-frequency, uses the components already in the proposal (SiN star
coupler, MZI trees, PZT or AlN shifters), and needs no new physics. A' is a useful fallback with fewer
optical elements but needs amplitude modulation on every arm and power headroom. D is excluded by the
chromatic focal shift unless the meta-objective is made achromatic over several nanometres around 1064 nm.

## 4. Route E: design and simulated performance

**Geometry.** 128 input ports, 1.7 um wide at 2.0 um pitch on the input arc; FPR radius
R = N p_a p_port / lambda_s = 886 um, which makes one port per resolvable direction; 128 arms (1.7 um wide,
2.0 um pitch at the FPR); 128 end-fire emitters at 1.5 um pitch; the slit-filtered relay and NA 0.5
meta-objective of `docs/findings.md`. Adjacent ports at 0.3 um gaps have a coupling length of 625 um, so the
~20 um over which they converge on the FPR transfers below -20 dB (40 um would still be within -20 dB).

**Window width.** The Gaussian window in port space sets the spot apodisation (uncalibrated values):

| sigma (ports) | FWHM | position nonlinearity | peak ripple | neighbour / peak at 4 um |
|---|---|---|---|---|
| 0.8 | 1.51 um | 18 nm | 3.4% | 1.1e-4 |
| 0.9 | 1.66 um | 8.6 nm | 1.5% | 2.0e-5 |
| 1.0 | 1.82 um | 3.7 nm | 0.55% | 2e-6 |
| 1.2 | 2.16 um | 0.5 nm | 0.07% | 7.6e-5 |

The nonlinearity and ripple are deterministic functions of the window position and are removed by a lookup
table, together with the fall-off of delivered power from 0.40 at the centre to 0.19 at +/-64 um. We use
sigma = 0.9 as a compromise between spot size (power per depth) and the residual that calibration must remove.

**Independence.** Moving one of 20 tweezers by up to 4 um changes the other 19 by at most 1.5 nm and 0.19%;
removing its window leaves 8.9e-6 at its site. Two tweezers approach without interference if their relative
phase is set to pi/2, since the cross term of two real fields then vanishes: at 2.5 um the pull is 3.5 nm and
the peak change 0.16%, against 146 nm and 9.5% in phase. This is a control an AOD lacks, where approaching
tones beat at their frequency difference.

**Actuator accuracy.** Each tweezer is set by about five ports, so per-port errors average little:

| per-port amplitude / phase error (rms) | position noise | peak noise |
|---|---|---|
| 0.3% / 3 mrad | 3.1 nm | 0.5% |
| 1% / 10 mrad | 9.4 nm | 1.9% |
| 3% / 30 mrad | 28 nm | 5.6% |

This is about 15 times tighter than the arm-side OPA (0.14 rad for 10 nm at N = 64, `m3_noise`). PZT
hysteresis must therefore be corrected, for example by closed-loop calibration against on-chip taps.

**Speed.** A moving tweezer is updated through the tree. The staircase result of `m3_update_rate` applies:
>= 200 kHz updates at 1.2 mK. PZT shifters (0.4 MHz [MONT]) meet this; thermo-optic shifters (34 us rise)
do not. The node settings follow in closed form from the desired port amplitudes (each node splits in the
ratio of the powers of its two subtrees), so an FPGA can compute them in real time.

**Power.** Facet -0.37 dB, 7-stage tree -0.39 dB, propagation -0.3 dB (assumed), meta-objective 0.7
(assumed), and the simulated FPR, arms, emitters, slit and pupil. A 1.2 mK tweezer needs 48 mW at the fibre at
the field centre and 102 mW at the edge; equalising depth over the field costs a mean factor of 1.34, so 20
tweezers spread over +/-60 um need about 1.3 W. Restricting operation to +/-40 um lowers the mean factor.

**On/off.** A tweezer is ramped on or off by scaling its window through the tree (adiabatic, >= 100 us), with
the residual set by the tree extinction; the simulated residual with the window removed is 8.9e-6.

## 5. Where this leaves the 894.6 nm addressing beams

Raman gates act on atoms held in static traps, and the addressed Raman leg must keep a fixed frequency
relative to the global leg, so frequency tagging is not available at 894.6 nm. The switched-port star
coupler of `docs/twenty_tweezers.md` (20 ports, two cascaded AlN MZIs each, neighbour intensity 3.5e-6)
remains the recommended addressing engine. Route E at 894.6 nm would give steerable addressing, but its
off-state would be set by the extinction of the tree nodes rather than by a cascaded switch per beam.

## 6. Limits and next steps

The FPR model is 2D and paraxial with isolated port and arm modes; dense ports at 0.3 um gaps couple, so the
port array and FPR need EME or 3D beam-propagation design, as for AWG star couplers. The relay is scalar and
aberration-free. Steering is along one row; the second axis needs a second engine or an OPA-AOD hybrid.
Next steps: (i) EME/BPM design of the dense-port FPR; (ii) a closed-loop calibration scheme for 255
actuators; (iii) a rearrangement-sequence simulation (atom survival per move, including the actuator noise
spectrum and the >= 200 kHz update).

## 7. Suggested wording for the proposal

Hypothesis H3b: "A focal-plane OPA, in which each tweezer is a sliding window over a dense array of input
ports, moves 20 tweezers independently along a row and switches each on and off, as the tones of an AOD do,
but at one optical frequency, so that tweezers can approach within 2.5 um without beating or intermodulation."
Table 1, OPA row: "independently steerable, individually switchable tweezers >= 20 on one row; update rate
>= 200 kHz; position noise <= 10 nm rms".

## Verification notes

No new literature values were used except the recalled surface-acoustic-wave velocity (route G, order of
magnitude only). [17] is the proposal's reference for crossed AODs; [MONT] is arXiv:2601.15695 (read in
session). The AOD physics in Sec. 2 is stated without a citation: [REF?].

| value | location | type | pointer |
|---|---|---|---|
| independence, switch-off, approach, spread | Sec. 1, 4 | SIM | results/m9_E_independence.json (review before submission) |
| window table, field fall-off 0.40 to 0.19 | Sec. 4 | SIM | results/m9_E_window.json (review) |
| noise table | Sec. 4 | SIM | results/m9_E_independence.json, amplitude_phase_noise (review) |
| 48 / 64 / 102 mW, 1.3 W | Sec. 1, 4 | SIM/DERIVED | results/m9_E_power.json (review) |
| R = 886 um, 1.064 um per port | Sec. 4 | DERIVED | m9_E_window.json; tests/test_star.py |
| coupling length 625 um (1.7 um ports, 0.3 um gap) | Sec. 4 | SIM | femwell, run in session (not stored) |
| A' 0.75, F 0.44 | Sec. 3 | SIM | results/m9_papr.json (review) |
| D: 14 um/nm, +/-71 pm, 6.7 m; G bandwidths; T rates | Sec. 3 | DERIVED | results/m9_routes.json |
| v_SAW ~5 km/s | Sec. 3 | RECALLED | verify for an AlN/SiN/SiO2 stack |
| port pitch 2.0 um, port width 1.7 um, sigma 0.9, propagation 0.3 dB/cm, meta-objective 0.7 | Sec. 4 | CHOICE | stated assumptions |

Open placeholders: [REF?] AOD multi-tone independence and intermodulation; [VALUE?] PZT hysteresis and
long-term drift after closed-loop calibration; [VALUE?] dense-port FPR excess loss from EME/BPM.
