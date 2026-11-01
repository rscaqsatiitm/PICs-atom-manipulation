# Phase-A simulations of the SiN PIC elements: first results

Proposal ANRF/ARG/2026/009937 (draft v97), Objective 3 and the Phase-A error budget.
Code: `picsim/`, scripts `scripts/m1` to `m7`, results in `results/*.json` (inputs with provenance, outputs).
Every simulated value below is a SIM output, valid only within the regime simulated, and must be reviewed
before it enters the proposal.

## 1. Summary

The passive SiN elements meet the proposal's targets with margin on a 400 nm LPCVD film. Five results bear
directly on the proposal text.

(i) **OPA update rate.** A tweezer moved by a staircase of OPA updates is excited when the update rate
aliases the move spectrum onto the trap frequency. At the 1.2 mK transport depth (radial trap frequency
84 kHz for w0 = 1.04 um), a 100 kHz update adds about 40 motional quanta to a 100 us, 4 um move and about
2,900 quanta to a 50 us move, against 0.005 and 0.2 for a smooth trajectory. At 200 kHz and above the
staircase adds nothing measurable. Hence the >= 100 kHz target of Table 1 holds for 0.3 mK tweezers
(42 kHz) but should read >= 200 kHz for moves at 1.2 mK.

(ii) **OPA field and tweezer count.** The number of resolvable positions along the steering axis equals the
channel number N, so the range at the atoms is N lambda/(4 NA F), with F the pupil fill. At NA 0.5 and full
fill, 64 channels steer over +/-33 um (17 sites at 4 um pitch), not +/-100 um, and cannot place more than
about 15 tweezers on one row. More than 20 simultaneous tweezers on one row need 128 channels; 128 channels
at 0.75 fill cover +/-89 um at 1.33 times the power per tweezer.

(iii) **Switch architecture.** A single MZI with the simulated 2x2 MMI spread reaches 50 dB in only 49% of
devices at 3 mrad phase accuracy. A cascaded pair, as in [ZHAO], reaches >= 50 dB in every trial
(5th percentile 66 to 88 dB). The binding constraints are then polarisation and free-space crosstalk, not
the interferometer.

(iv) **Addressing crosstalk is optical, not on-chip.** With a uniformly filled NA 0.5 pupil, the 894.6 nm
spot puts 1.0e-3 of its peak intensity on the 4 um neighbour, which gives a neighbour error of 2.5e-3,
above the 1e-3 budget. A Gaussian pupil fill truncated at 0.7 of the pupil radius brings this to 1.2e-4
(error 2.9e-4) at a spot FWHM of 1.13 um instead of 0.92 um.

(v) **Position noise must be specified as a spectral density.** The error-budget statement "position noise
<= 10 nm rms keeps transport survival >= 99%" does not fix the heating rate. If 10 nm rms is white over
1 MHz, a 1.2 mK tweezer heats at about 7,700 quanta/s (about 0.8 quanta per 100 us move, against a depth of
about 300 quanta). Heating below 10 quanta/s needs S_x(84 kHz) <= 1.3e-25 m^2/Hz.

The end-to-end loss of the end-fire OPA from fibre to atoms, for 21 tweezers at N = 128, is 5.9 dB (7.1 dB
for the best grating emitters), against the ~10 dB assumed in the proposal; the per-tweezer on-chip power
falls from 66 mW to about 26 mW. This estimate excludes sidewall scattering and assumes a meta-objective
transmission of 0.7.

## 2. Scope, methods and limits

We model the elements named in Sec. 4 (Phase C) of the proposal. The waveguide platform is stoichiometric
LPCVD Si3N4, 400 nm thick (the CoE-CPPICS limit), with 200 and 300 nm as foundry alternatives, in SiO2.
Material dispersion is from Luke et al. [LUKE] and Malitson [MAL], and the thermo-optic coefficients
(Si3N4 2.45e-5 /K, SiO2 1e-5 /K) from [TRC]. The working wavelengths are 894.6 nm (addressing) and
1064.4 nm (transport).

The methods are: vector FEM mode solving (femwell 0.1.12, second-order Nedelec elements) for straight, bent,
coupled and leaky guides; effective-index reduction with guided-mode eigenmode expansion for MMIs; scalar
overlap with Gaussian fibre modes for edge couplers; 2D FDTD (MEEP 1.35, 40 px/um) for grating antennas,
with the Si substrate included; 2x2 transfer matrices with Monte Carlo for switches; 2D thermal FEM for
heaters; 1D scalar Fourier optics for the OPA (emitter near field, Fourier-plane slit, aplanatic
objective); and an exact harmonic-oscillator result for transport heating. The mode, slab, MZI, OPA and
transport models are checked against analytic limits or published examples (`tests/`, 16 tests); the FDTD
closes its energy balance to better than 0.5%, and the heater result changes by 1.4% when the domain is
doubled.

The principal limits are as follows. The MMI and grating models are 2D. Edge-coupler numbers are mode
overlaps without Fresnel or taper conversion loss. The OPA model is scalar and 1D along the steering axis,
with an aberration-free objective. Sidewall scattering loss and fabrication roughness are not modelled. The
AlN actuator is a single mechanical mode with assumed f0 and Q.

## 3. Platform (Module 1)

The design width is the widest width (0.1 um grid) for which TE1 is cut off or within 0.015 of the cladding
index and the TE0/TM0 splitting is at least 0.015. The second condition matters: on 400 nm film at 894.6 nm,
a 0.4 um square guide is polarisation-degenerate (TE0 = TM0 = 1.693), which would couple polarisations in
bends. Hence we use 0.5 um there.

| Film | lambda (nm) | width (um) | n_eff | n_g | dn_eff/dT (1/K) | R for < 0.01 dB/90 deg, TE / TM (um) |
|---|---|---|---|---|---|---|
| 400 nm | 1064.4 | 0.6 | 1.708 | 2.14 | 2.36e-5 | 8 / 10 |
| 400 nm | 894.6 | 0.5 | 1.745 | 2.17 | 2.41e-5 | 6 / 6 |
| 200 nm | 1064.4 | 0.9 | 1.613 | 1.98 | 2.02e-5 | 20 / 60 |
| 200 nm | 894.6 | 0.7 | 1.651 | 2.05 | 2.16e-5 | 20 / 30 |

Substrate leakage on 400 nm film is below 1e-5 dB/cm for a buried oxide of 2 um or more at both
wavelengths; on 200 nm film at 1064.4 nm the oxide must be at least 3 um (TE 7e-8 dB/cm; at 2 um, TM leaks
0.33 dB/cm). The coupling length of two 1064.4 nm guides on 400 nm film is 1.8 mm at 1.6 um pitch and
29 mm at 2.1 um pitch. At 1.4 W the Kerr phase is 0.05 rad/cm and the peak in-guide intensity
420 MW/cm^2.

On 400 nm film the TM bend loss is only about ten times the TE loss at the same radius, hence bends cannot
strip TM light. On 200 nm film at 1064.4 nm a 20 um bend radiates 1.5 dB per 90 deg of TM against
5e-4 dB of TE, which makes a compact TM filter possible; at 894.6 nm the same film gives only 0.03 dB per
90 deg of TM.

## 4. Edge couplers and splitters (Module 2)

Inverse tapers on 400 nm film reach a mode overlap of -0.37 dB with a 2.5 um lensed fibre at 1064.4 nm
(150 nm tip, +/-0.6 um lateral tolerance at -1 dB) and -0.48 dB with HI1060 (100 nm tip, +/-1.55 um). At
894.6 nm the best values are -0.46 dB (lensed) and -1.1 dB (PM780). Several optima lie at the narrowest tip
swept (100 nm), so narrower tips may do better if lithography allows. These overlaps exclude Fresnel and
taper conversion loss and are therefore lower bounds on the 1 to 2 dB per facet of the proposal. At 1.4 W
through one 400 nm facet the intensity is 45 MW/cm^2; three input fibres keep it below 20 MW/cm^2. We
did not find a CW facet-damage threshold for SiN at 1 um in this session.

A 4 um wide 1x2 MMI has 0.044 dB excess loss at both wavelengths (lengths 14.8 um at 1064.4 nm and 17.7 um
at 894.6 nm), so the 6-stage, 64-way tree costs 0.27 dB, against the ~1 dB assumed. A symmetric 1x2 MMI is
balanced by symmetry; imbalance comes only from asymmetric errors, 0.13 dB for a 20 nm lateral offset of the
input guide. With 1-sigma errors of 10 nm (width), 5 nm (thickness) and 5 nm (offset), the 2x2 split ratio
varies by 0.0024 (894.6 nm) and 0.0034 (1064.4 nm). Our first-pass 894.6 nm 2x2 design sits at 0.490
rather than 0.500; this offset is removable by redesign and is carried separately below.

## 5. Switches and router (Modules 4 and 5)

| Architecture (894.6 nm, split spread 0.0024) | phase error 3 mrad: median / 5th pct (dB) | P(>= 50 dB) | phase error 10 mrad: median / 5th pct (dB) | P(>= 50 dB) |
|---|---|---|---|---|
| single MZI | 49.8 / 41.4 | 0.49 | 45.3 / 38.8 | 0.21 |
| single MZI, trimmed input splitter | 55.1 / 48.8 | 0.89 | 44.6 / 38.3 | 0.18 |
| two cascaded MZIs | 100 / 88 | 1.00 | 92 / 81 | 1.00 |

Arm amplitude imbalance (0.5% rms) is included. With the 0.010 design offset added, a single MZI falls to
P = 0.15 while the cascade keeps P = 1.00 (5th percentile 66 dB). An incoherent stray-light floor of -70 dB
caps the cascade at 70 dB. A single balanced MZI with a phase error delta has ER = -20 log10(delta/2), so
50 dB needs delta <= 6.3 mrad; for the AlN half-wave voltage at 894.6 nm, which we bracket between 84 V
(lambda-scaled from 74.7 V at 795 nm [ZHAO]) and 132 V (interpolated to 200 V at 1013 nm [ZHAO]), the
driver noise must stay below 0.17 to 0.27 V rms.

Polarisation sets a separate floor. TM light is not nulled when the TE path is. For an input polarisation
extinction of 20 dB, 50 dB at the output needs 24 to 30 dB of TM suppression after the chip; for 25 dB
input, 19 to 25 dB. On 400 nm film bends cannot provide this (Sec. 3), hence a polariser in the free-space
delivery path or a dedicated on-chip TE-pass element is required.

The settling of the switch is set by mechanical ringing, not by the 26 ns rise [ZHAO]. For an assumed single
mode at f0 = 20 MHz and Q = 10, an unshaped step needs 0.98 us to reach -50 dB after switch-off. A
two-impulse shaper (ZV) gives 24 ns but 0.43 us with a 2% error in f0; a three-impulse shaper (ZVD) gives
49 ns with a 2% error and 0.19 us with 5%, at a 28 ns optical rise. Hence <= 100 ns needs ZVD shaping with
f0 known to about 2%. These f0 and Q values are assumptions to be replaced by measured ones.

The 16-output router (four stages of thermo-optically tuned MZI splitters) leaks at most 1.4e-3 into the
worst unaddressed output (95th percentile) for a 30 mrad trim error, so with the cascaded output switch the
idle-site extinction exceeds 90 dB. Its insertion loss is 5.0 dB, of which 3 dB is the assumed free-space
delivery, so 16 sites at 20 uW need about 1 mW at the input.

Thermo-optic shifters with a 2 um wide heater 2 um above the core need 84 mW for pi at 1064.4 nm (68 mW at
894.6 nm), with a 34 us rise and 16% thermal crosstalk to a guide 10 um away; isolation trenches reduce
these to 54 mW and 2%, and an additional undercut to 3 mW, with a 1.5 ms rise. They suit the slow router
trims but not OPA updates at >= 100 kHz. For the OPA, the PZT shifters of [MONT] (V_pi 2.8 V at 532 nm,
about 5.6 V at 1064.4 nm by wavelength scaling, 0.4 MHz bandwidth) settle in about 2 us, within a 5 us
update period.

## 6. Optical phased array (Module 3)

**Range and resolution.** The atom-plane field is a scaled copy of the slit-filtered OPA far field.

| N | fill | usable half-range (um) | sites on one row (4 um) | spot FWHM along x (um) | power for equal depth |
|---|---|---|---|---|---|
| 32 | 1.0 | 16 | 9 | 0.95 | 1.00 |
| 64 | 1.0 | 33 | 17 | 0.95 | 1.00 |
| 64 | 0.5 | 68 | 35 | 1.89 | 1.99 |
| 128 | 1.0 | 67 | 33 | 0.94 | 1.00 |
| 128 | 0.75 | 89 | 45 | 1.26 | 1.33 |
| 256 | 1.0 | 134 | 67 | 0.94 | 1.00 |

"Usable" means the steered peak keeps at least half its on-axis value. Uniform channel amplitudes give 4.7%
first sidelobes along x; a Gaussian channel apodisation (1/e^2 at the aperture edge) lowers them to 0.9% at
a FWHM of 1.09 um.

**Emitter pitch.** A Fourier-plane slit passes only the central diffraction zone, so the efficiency into it
falls with pitch: 0.65 at 1.5 um, 0.52 at 2 um, 0.34 at 3 um (for the 0.6 um emitter guide; a wider emitter
mode at the facet would raise these). The minimum pitch is set by coupling between neighbouring emitters
over their parallel length, for < -20 dB power transfer:

| Emitter | parallel length (um) | minimum pitch (um) | upward fraction | emission angle |
|---|---|---|---|---|
| end-fire facet | 50 (assumed) | 1.46 | 1 (in plane) | 0 deg, wavelength independent |
| grating, 200 nm etch, 2 um oxide | 20 | 1.30 | 0.76 | 8.9 deg |
| grating, 100 nm etch, 3 um oxide | 90 | 1.56 | 0.72 | 9.7 deg |
| grating, 50 nm etch, 3 um oxide | 384 | 1.82 | 0.68 | 10.0 deg |

The grating directionality oscillates with oxide thickness through the substrate reflection (0.42 to 0.76
over 2 to 4 um), and the emission angle moves by -0.10 deg/nm. The 200 nm etch reflects 6% back into the
guide. Hence the end-fire array is the better first choice: denser pitch for short parallel sections, no
directionality loss, and no wavelength dependence of the pointing. The grating array becomes attractive if
a bottom reflector or two-layer grating is available.

**Simultaneous tweezers.** With weighted Gerchberg-Saxton alone, light at untargeted lattice sites reaches
11 to 21% of the tweezer intensity. A refinement that also penalises stray light and off-centre spots
holds it at 3 to 7%:

| N | K | max stray at untargeted sites | max (Imax - Imin)/(Imax + Imin) | power in tweezers / emitted |
|---|---|---|---|---|
| 64 | 8 | 0.056 | 0.062 | 0.44 |
| 64 | 12 | 0.053 | 0.061 | 0.46 |
| 128 | 16 | 0.054 | 0.106 | 0.45 |
| 128 | 21 | 0.057 | 0.106 | 0.45 |

(Four random target sets per case, one row along the steering axis.) A 5.7% stray at a static site
corresponds to about 70 uK under a 1.2 mK transport tweezer.

**Phase accuracy.** For independent channel phase errors sigma, the position jitter at the atoms is
sigma sqrt(12/N) lambda/(4 pi NA) (Monte Carlo agrees to 2%). Hence 10 nm needs sigma <= 0.14 rad at
N = 64 (0.10 rad at 32, 0.19 rad at 128); 8-bit phase resolution (7 mrad) is sufficient, and a PZT
hysteresis error must stay below about 4% of the applied phase at N = 64.

**Loss budget (end-fire, N = 128, K = 21, 2 um pitch).** Edge coupler -0.37 dB, MMI tree -0.27 dB,
propagation 1 cm at 0.3 dB/cm (assumed) -0.30 dB, emitters with slit, pupil and synthesis -3.43 dB,
meta-objective transmission 0.7 (assumed) -1.55 dB: total -5.9 dB. The grating case adds -1.2 dB of
directionality. At 1.5 um pitch the slit term improves by about 1 dB.

## 7. Error-budget links and power handling (Modules 6 and 7)

The proposal's footnote 4 is reproduced: an idle-qubit error of (pi^2/4) 10^(-ER/10), 2.5e-5 at 50 dB when
one Raman beam is addressed. The trap depth of 6.6 mW in a Gaussian of w0 = 1.04 um with 1160 a.u. is
1.01 mK, not 1.2 mK; 1.2 mK needs 7.8 mW, or w0 = 0.96 um. The OPA spot along the steering axis is 22%
brighter per watt than that Gaussian, which brings 6.6 mW to about 1.24 mK.

Intensity noise must satisfy S_I(2 nu) <= 1.4e-11 /Hz at 1.2 mK for a parametric heating time of 1 s. At
1.4 W into the input guide, with 0.3 dB/cm fully absorbed, the guide warms by 3.8 K and accrues
5.3 rad/cm of common-mode phase before the split; after the 64-way split each channel accrues at most
0.08 rad/cm, which the per-channel calibration absorbs.

## 8. Suggested changes to the proposal

(i) Table 1, OPA update rate: ">= 200 kHz" for transport at 1.2 mK. Alternatively, keep >= 100 kHz and
lower the transport depth to <= 0.4 mK during moves (trap frequency 48 kHz) with move times of >= 100 us,
for which a 4 um move adds 0.025 quanta (0.008 for a smooth trajectory).
(ii) Table 1, OPA tweezers: state ">20 tweezers with 128 channels" or ">=15 with 64 channels on one row",
and state the steering field (+/-33 um at N = 64, F = 1) separately from the +/-100 um meta-objective
field.
(iii) Sec. 4, "~10 dB of chip and OPA loss": the model gives about 6 dB (end-fire) to 7 dB (grating); keep
10 dB as the conservative value but note the modelled one.
(iv) Sec. 4 error budget: replace "position noise <= 10 nm rms" by a spectral density at the trap frequency,
S_x(nu_trap) <= 1.3e-25 m^2/Hz for < 10 quanta/s at 1.2 mK, and add the 894.6 nm pupil apodisation needed
for neighbour crosstalk <= 1e-4.
(v) Sec. 4(ii): 6.6 mW gives 1.0 mK for w0 = 1.04 um; quote 7.8 mW or the OPA spot.
(vi) Phase C: specify cascaded MZIs (or trimmed single MZIs with <= 3 mrad phase accuracy), an output
polariser of >= 30 dB on the 894.6 nm path, and shaped (ZVD) drive for <= 100 ns settling.

## 9. Next steps

(i) 3D FDTD or EME of the inverse taper, MMI and end-fire facet to replace the 2D and overlap estimates;
(ii) end-fire emitter with an expanded facet mode to raise slit efficiency at 1.5 um pitch;
(iii) mechanical FEM of the AlN actuator to replace the assumed f0 and Q;
(iv) a 2D model of the OPA spot through the meta-objective design of Phase B, including its aberrations,
since a Strehl of 0.8 adds a halo not included in the crosstalk estimate;
(v) GDS layout of the 64-channel OPA and 16-channel router on one 10 x 10 mm die.

---

## Verification notes

**References read in this session.**
[LUKE] K. Luke et al., Opt. Lett. 40, 4823 (2015), doi:10.1364/OL.40.004823 (Sellmeier via refractiveindex.info).
[MAL] I. H. Malitson, JOSA 55, 1205 (1965), doi:10.1364/JOSA.55.001205 (via refractiveindex.info).
[TRC] arXiv:2606.05673, introduction: dn/dT of Si3N4 and SiO2 at 193 THz.
[ZHAO] M. Zhao et al., arXiv:2508.09920: 71.4 dB, 26 ns, V_pi 74.7 V (795 nm) and 200 V (1013 nm),
cascaded MZIs, -68 dB on-chip and -50.8 dB free-space crosstalk.
[MONT] N. Montifiore et al., arXiv:2601.15695: PZT coil MZM, V_pi 2.8 V at 532 nm, 0.4 MHz, 21.5 dB.
[OHARA] K. M. O'Hara et al., arXiv:physics/0003049: Gamma = pi^2 nu^2 S_I(2 nu), Q_dot = 4 pi^4 M nu^4 S_x(nu).
AlN index: Pastrnak and Roskovcova, Phys. Stat. Sol. 14, K5 (1966), via refractiveindex.info.

**Recalled, not verified (check before use).** Weighted Gerchberg-Saxton: R. Di Leonardo et al., Opt.
Express 15, 1913 (2007) (publisher page blocked); fibre MFDs (HI1060 6.2 um, PM780 5.3 um); Si index 3.55;
Kerr index 2.5e-19 m^2/W; thermal conductivities and heat capacities; Cs mass.

**Provenance ledger (quoted values).**

| value | location | type | pointer |
|---|---|---|---|
| 40 and 2,900 quanta at 100 kHz; 200 kHz threshold | Sec. 1(i), 8(i) | SIM | results/m3_update_rate.json (review before submission) |
| 0.025 and 0.008 quanta at 0.4 mK | Sec. 8(i) | SIM | picsim.transport.n_bar_staircase, run in session (review) |
| 84 kHz, 42 kHz trap frequencies | Sec. 1(i) | DERIVED | picsim.transport.radial_freq_Hz, w0 = 1.04 um (USER) |
| +/-33, 67, 89 um ranges; 17, 33, 45 sites | Sec. 1(ii), 6 | SIM | results/m3_field_tradeoff.json (review) |
| 49%, 89%, P = 1.00, 66-88 dB | Sec. 1(iii), 5 | SIM | results/m4_extinction.json (review) |
| 1.0e-3, 1.2e-4 neighbour intensity; 2.5e-3, 2.9e-4 error | Sec. 1(iv) | SIM | results/m7_budget.json (review) |
| 7,700 quanta/s; 1.3e-25 m^2/Hz | Sec. 1(v), 8(iv) | DERIVED | [OHARA] formula, results/m7_budget.json |
| -5.9 dB, -7.1 dB, 26 mW | Sec. 1, 6 | SIM | results/m3_loss_budget.json (review) |
| design widths, n_eff, n_g, dn/dT, bend radii | Sec. 3 | SIM | results/m1_*.json (review) |
| -0.37 to -1.1 dB overlap, 45 MW/cm^2 | Sec. 4 | SIM | results/m2_couplers.json (review) |
| 0.044 dB per MMI, 0.27 dB tree, 0.0024/0.0034 spread | Sec. 4 | SIM | results/m2_mmi.json (review) |
| 84-132 V, 0.17-0.27 V rms | Sec. 5 | DERIVED | from [ZHAO] values, results/m4_dynamics.json |
| 0.98 us, 24 ns, 49 ns, 0.19 us | Sec. 5 | SIM | results/m4_dynamics.json (f0, Q assumed; review) |
| 1.4e-3 leak, > 90 dB, 5.0 dB, 1 mW | Sec. 5 | SIM | results/m4_router.json (review) |
| 84, 54, 3 mW; 34 us, 1.5 ms; 16%, 2% | Sec. 5 | SIM | results/m5_heaters.json (review) |
| min pitch 1.30-1.82 um; slit efficiency 0.65/0.52/0.34 | Sec. 6 | SIM | results/m3_pitch_efficiency.json (review) |
| directionality 0.42-0.76; -0.10 deg/nm; 6% reflection | Sec. 6 | SIM | results/m3_grating.json, m3_grating_fdtd_raw.json (review) |
| stray 3-7%; uniformity to 0.106 | Sec. 6 | SIM | results/m3_multibeam.json (review) |
| sigma 0.10/0.14/0.19 rad; 4% hysteresis | Sec. 6 | SIM/DERIVED | results/m3_noise.json (review) |
| 1.01 mK, 7.8 mW, 0.96 um, 1.22 | Sec. 7 | DERIVED/SIM | results/m7_budget.json |
| 3.8 K, 5.3 rad/cm, 0.08 rad/cm | Sec. 7 | SIM | results/m6_power.json (review) |
| 6.6 mW, 1.2 mK, 1.04 um, 1160 a.u., 10 dB, 1.4 W, >= 50 dB, <= 100 ns, 10 nm, NA 0.5, 4 um | throughout | USER | proposal draft v97 |
| meta-objective transmission 0.7; 0.3 dB/cm; 50 um end-fire section; f0 = 20 MHz, Q = 10; fab sigmas | Sec. 5, 6 | CHOICE | stated assumptions, not claims |

**Open placeholders.** [VALUE?] CW facet-damage threshold of SiN/SiO2 at 1 um; [VALUE?] PZT hysteresis
magnitude; [VALUE?] AlN actuator f0 and Q at 894.6 nm; [VALUE?] measured V_pi at 894.6 nm; [REF?] verified
citation for weighted Gerchberg-Saxton.
