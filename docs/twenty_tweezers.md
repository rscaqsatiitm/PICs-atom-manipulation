# Twenty steerable, individually switchable tweezers at 1064.4 and 894.6 nm

Design study for Objective 3 (proposal ANRF/ARG/2026/009937, draft v97). Code: `picsim/opa/star.py`,
`scripts/m8_twenty.py`; results `results/m8_*.json`. All simulated values are SIM outputs, valid in the
regime simulated, and must be reviewed before they enter the proposal.

## 1. Result

Twenty tweezers that are steerable and individually switchable can be delivered on the 400 nm SiN platform
by a multi-port star-coupler OPA with 128 emitters, one chip per wavelength. Each tweezer has its own input
port and its own MZI switch; the ports illuminate a common 128-arm phased array through a free-propagation
region (FPR); a phase ramp on the arms moves all lit tweezers together. On simulated designs, with the 20
ports placed on a row of 4 um sites:

| | 1064.4 nm (transport tweezers) | 894.6 nm (addressing beams) |
|---|---|---|
| spot FWHM along the row | 1.52 um | 1.27 um |
| intensity at the 4 um neighbour / peak | 1.0e-6 | 3.5e-6 (neighbour error 8.5e-6) |
| peak spread, all 20 lit, random port phases | <= 0.73% | <= 0.89% |
| stray at unlit sites (10 of 20 lit) | <= 1.2e-5 | <= 2.1e-5 |
| static position error of a spot | <= 6 nm | <= 5 nm |
| switch per port | one trimmed thermo-optic MZI | two cascaded AlN MZIs |
| extinction needed | >= 31 dB (1 uK residual at 1.2 mK) | >= 50 dB (idle error 2.5e-5) |
| fibre power | 33-43 mW per tweezer at 1.2 mK; 0.74 W for 20 | <= 0.14 mW per site; 2.4 mW for 20 |

The price is that steering is collective: all lit tweezers move by the same displacement. Which sites carry
tweezers is set by which ports are lit. This is the primitive an AOD provides for parallel rearrangement
(select a subset of a row, move it), but at one optical frequency, so two tweezers brought together do not
beat. Fully independent steering of 20 tweezers through one aperture requires phase-only synthesis on a
shared OPA, and that architecture is a poor on/off engine (Sec. 4).

## 2. Architecture

One chip per wavelength, both on the shared 400 nm LPCVD SiN process:

(i) edge-coupled input (two fibres at 1064.4 nm, to keep each facet below 20 MW/cm^2);
(ii) a 1x32 tunable-splitter tree (five stages of thermo-optic MZIs with 3 um 2x2 MMIs), which sends power
only to lit ports, so the power per tweezer does not depend on how many are lit;
(iii) 20 port switches, thermo-optic at 1064.4 nm (tweezer on/off ramps of >= 100 us), two cascaded AlN MZIs
at 894.6 nm (gate pulses, ZVD-shaped drive, <= 100 ns);
(iv) 20 input tapers, 1.4 um (1064.4 nm) or 1.0 um (894.6 nm) wide, spaced 1.89 um and 1.64 um at the FPR;
(v) an FPR of radius 223 um (1064.4 nm) or 197 um (894.6 nm);
(vi) 128 arms, 1.7 um wide at the FPR edge (0.3 um gaps), each with a phase shifter, routed to
(vii) 128 end-fire emitters at 1.5 um pitch;
(viii) the relay of Sec. 6 of `docs/findings.md`: cylindrical collimation, Fourier-plane slit, anamorphic
telescope, meta-objective (NA 0.5).

The two chips are combined before the meta-objective by a dichroic. Their port geometries map to the same
4 um sites; the residual registration is set by the common phase ramp of each chip.

Each port illuminates all 128 arms with the same envelope, the far field of its taper mode, which acts as a
Gaussian-like apodisation of the array. Port j imposes a linear phase k_s d_j p_a / R per arm (k_s, the slab
wavenumber; d_j, port offset; p_a, arm pitch), which steers its beam to sin(theta_j) = n_s d_j p_a / (R p_out).
The port pitch that places beams 4 um apart at the atoms follows from the relay scale
x = sin(theta) N p_out / (2 NA): a 4 um step is sin(theta) = 0.0208 for N = 128 and p_out = 1.5 um.

## 3. Design choices and trade-offs

**Apodisation.** The arm-power envelope sets the spot tails, and the tails set two things: the light that a
lit tweezer leaves on its neighbours, and the interference between lit ports, whose relative phases drift
with the router temperature. With the envelope truncated at 1/e^2 at the outer arms (as for a standard star
coupler), the neighbour intensity is 1.2e-3 and, with all 20 lit and random port phases, peak intensities
vary by up to 19%. Narrowing the envelope removes both at the cost of a wider spot:

| envelope at outer arm (nominal) | FWHM at 1064.4 nm | neighbour / peak | peak spread, 20 lit |
|---|---|---|---|
| e^-2 | 1.10 um | 1.2e-3 | 19% |
| e^-4 | 1.24 um | 4.2e-4 | 11% |
| e^-6 | 1.38 um | 9.6e-5 | 5.8% |
| e^-8 | 1.52 um | 1.3e-5 | 2.0% |

(Input taper 0.6 um.) A narrower envelope needs a smaller FPR radius, which brings the ports closer together
(1.0 um at e^-8). Widening the input taper restores the radius: with 1.4 um tapers the port pitch is
1.89 um, and the taper mode, being closer to Gaussian than the 0.6 um guide mode, also lowers the tails
further. Two 1.4 um tapers 0.49 um apart have a coupling length of 1.33 mm at 1064.4 nm, so the ~20 um over
which they converge on the FPR transfers well below -20 dB. The final designs were chosen from a scan of
taper width and envelope (`results/m8_taper_scan.json`): 1.4 um tapers at 1064.4 nm, 1.0 um tapers at
894.6 nm, both at a nominal e^-6 envelope. Because the wide-taper modes are not Gaussian, the realised
envelope at the outer arm is 2.7e-5 (1064.4 nm) and 7.2e-6 (894.6 nm), and a nominally narrower envelope
(e^-8) puts the outer arms on the side lobes of the port far field (realised 1e-3 to 2.5e-3, worse tails).
The designs are therefore specified by FPR radius and taper width, not by the nominal edge level.

We keep the cleaner design at 1064.4 nm although a wider envelope (e^-5) gives 9% more depth per watt,
because its 3% peak spread and 4e-5 neighbour light would load the static atoms and make tweezer depths
depend on drifting port phases.

**Emitter pitch.** 1.5 um, the smallest pitch at which a 50 um straight end-fire section keeps
emitter-to-emitter transfer below -20 dB (`m3_pitch_efficiency`). It raises the slit efficiency to 0.65
from 0.52 at 2 um.

**Channel count.** With N = 96 the steering half-range is 51 um (1064.4 nm) or 43 um (894.6 nm), so a row of
20 sites (+/-38 um) leaves only 5 to 13 um for collective moves; N = 128 gives +/-68 um and +/-57 um.

## 4. Why not phase synthesis on one OPA

On a shared 128-channel OPA, 20 tweezers can be computed with independent positions (`m3_multibeam`), but a
tweezer is switched off only by re-synthesis. Removing one of 20 (warm-started refinement with that site
among the dark points) leaves 1.9 to 3.3% of the tweezer intensity at its site (-15 to -17 dB), changes the
19 kept tweezers by up to 12%, and along a linear phase path between the two patterns the kept tweezers
deviate by up to 19%. Stray light at untargeted sites stays at about 2%. Hence it cannot gate 894.6 nm light
at 50 dB, and it perturbs the other atoms whenever one tweezer is switched. It remains the route to fully
independent steering at 1064.4 nm if that is ever required.

## 5. Steering

All lit tweezers move by the displacement set by a linear phase ramp across the 128 arms. With all 20 lit
(peak change quoted at 1064.4 nm; at 894.6 nm it is within 0.6% of these values):

| collective shift | position error (static, max) | peak change (min / max) |
|---|---|---|
| 4 um | 25 nm (1064.4), 5 nm (894.6) | -2.5% / +2.3% |
| 8 um | 25 nm, 5 nm | -5.1% / +4.3% |
| 16 um | 25 nm, 6 nm | -10.7% / +7.8% |
| 24 um | 26 nm (1064.4 only) | -16.7% / +10.3% |

The position error is a static offset of the outer spots, caused by the element-factor slope across them; it
changes by less than 2 nm over a 24 um sweep and is calibrated out. The peak change follows the emitter
element factor and can be pre-compensated through the port powers. At 894.6 nm the outer spots leave the
+/-57 um field beyond a 19 um shift; subsets of the row can move farther.

Because the collective move is one linear ramp, its fast part needs a single waveform: the 128 arm shifters
take fixed gains n = 0 ... 127 from one analog ramp, with slow per-arm offsets for calibration. A continuous
analog ramp has no staircase, which removes the update-rate heating of `m3_update_rate` (40 quanta per
100 us move at 100 kHz and 1.2 mK). The shifters must be fast and low-power: PZT on SiN (0.4 MHz,
V_pi about 5.6 V at 1064.4 nm by wavelength scaling of [MONT]) rather than thermo-optic (34 us rise).

## 6. Switching

At 1064.4 nm a tweezer is ramped on or off adiabatically over >= 100 us, so a thermo-optic MZI (34 us rise)
suffices; >= 31 dB keeps the residual below 1 uK at 1.2 mK, which a trimmed single MZI exceeds at phase
errors up to 10 mrad (5th percentile 38 dB, `m4_extinction`; computed at 894.6 nm, to be repeated at
1064.4 nm). At 894.6 nm the switch gates Raman pulses. The proposal's 50 dB gives an idle-qubit error of
2.5e-5; the switch leak equals the 3.5e-6 optical tail of a lit neighbour at 54.6 dB. Two cascaded AlN MZIs
give a 5th percentile of 66 to 88 dB, so the switch is not the limit; ZVD-shaped drive reaches -50 dB within
49 ns for a 2% error in the actuator resonance (`m4_dynamics`, with an assumed resonance).

## 7. Power

At 1064.4 nm, per tweezer, from the fibre: facet -0.37 dB, router -0.28 dB, switch -0.06 dB, propagation
(1 cm, 0.3 dB/cm assumed) -0.30 dB, FPR, arms, emitters, slit and pupil 0.40 to 0.52 of the port power
depending on the port, meta-objective 0.7 (assumed): -6.3 to -7.4 dB to the atom plane. For a 1.2 mK depth
(1160 a.u., Gaussian 1.04 um along the cross axis) this needs 33 to 43 mW at the fibre per tweezer, 0.74 W
for 20, against 1.4 W in the proposal. At 894.6 nm, 20 sites at the intensity of 20 uW in a 1.04 um Gaussian
need 2.4 mW in total.

## 8. Implementation load

| item | 1064.4 nm chip | 894.6 nm chip |
|---|---|---|
| router MZIs (thermo-optic) | 31 | 31 |
| port switches | 20 thermo-optic | 40 AlN (two per port, high voltage) |
| arm phase shifters | 128 PZT (one fast ramp + slow offsets) | 128 (slow; registration only) |
| FPR radius | 223 um | 197 um |
| electrical channels | about 180 | about 200 |

About 380 channels across the two chips exceed the 32-64-channel OPA and 16-channel router of the proposal,
and the 32-channel high-voltage driver of Sec. 8 covers only the 894.6 nm switches if one AlN MZI per port is
replaced by a trimmed thermo-optic stage (a trimmed single MZI reaches 50 dB in 89% of devices at 3 mrad).
Thermo-optic power, at half of P_pi on average, is about 42 mW per heater without trenches (`m5_heaters`),
i.e. about 2.1 W for 51 heaters on the 1064.4 nm chip; trenches with undercut bring this to about 0.1 W.

## 9. What this does not provide, and limits of the model

It does not move different tweezers by different amounts in one step; a sorting sequence is built from
parallel moves, as with AODs. It steers along one row; the second axis needs a second chip, a row-select
element, or an OPA-AOD hybrid (`docs/findings.md`, Sec. 6). The FPR model is 2D and paraxial (flat arcs, no
Rowland aberration, isolated arm modes), and the arm phases from fabrication are assumed calibrated. The
relay is scalar and aberration-free; a meta-objective with a Strehl of 0.8 adds a halo that is not included
in the crosstalk figures.

## 10. Suggested wording for the proposal

Hypothesis H3b: "An OPA fed through a multi-port star coupler forms up to 20 tweezers, each switched on and
off by its own integrated MZI, and moves all lit tweezers together by setting a single phase ramp across its
128 emitters. All tweezers share one optical frequency, hence ... ." Table 1, OPA row: "simultaneous,
individually switchable tweezers >= 20 on one row; collective move >= +/-16 um with all 20 lit; neighbour
intensity <= 1e-4 at 894.6 nm".

## Verification notes

Sources read in session: as listed in `docs/findings.md` ([MONT] arXiv:2601.15695 for the PZT figures).
No new literature values were introduced in this study; all new numbers are SIM or DERIVED.

| value | location | type | pointer |
|---|---|---|---|
| FWHM, neighbour, spread, stray, position error, delivery | Sec. 1, 3, 5 | SIM | results/m8_final.json (review before submission) |
| apodisation table | Sec. 3 | SIM | results/m8_apod.json (review) |
| taper-width scan, realised envelope 2.7e-5 / 7.2e-6 | Sec. 3 | SIM | results/m8_taper_scan.json, m8_final.json (review) |
| coupling length 1.33 mm (0.49 um gap, 1.4 um tapers) | Sec. 3 | SIM | femwell, run in session; noted in m8_final.json |
| re-synthesis residual 1.9-3.3%, 12%, 19% | Sec. 4 | SIM | results/m8_resynth.json (review) |
| collective-shift table | Sec. 5 | SIM | results/m8_final.json, steer (review) |
| 33-43 mW, 0.74 W, 2.4 mW | Sec. 1, 7 | SIM/DERIVED | results/m8_final.json, power (review) |
| 31 dB, 50 dB, 54.6 dB | Sec. 6 | DERIVED | results/m8_final.json, switching |
| emitter pitch 1.5 um; arm pitch 2.0 um; arm width 1.7 um; propagation 0.3 dB/cm; meta-objective 0.7 | Sec. 2, 7 | CHOICE | stated assumptions |
| 1.2 mK, 1160 a.u., 1.04 um, 4 um sites, NA 0.5, 20 uW per site | throughout | USER | proposal draft v97 |

Open placeholders: [VALUE?] PZT shifter footprint and hysteresis for 128 arms; [VALUE?] AlN resonance and Q
at 894.6 nm; [VALUE?] measured arm phase errors after fabrication (calibration range).
