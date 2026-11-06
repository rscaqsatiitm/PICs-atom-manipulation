# Phase-to-intensity conversion of the decorrelating phase noise

Follow-up to `docs/row_coherence.md`, which replaces per-chip frequency shifters by common Gaussian phase noise
(3 rad rms, correlation rate B = 2 GHz) on the seed. Code: `scripts/m16_pm_am.py`; results
`results/m16_pm_am.json` (+ .png). Values are SIM, valid for the element models and parameter ranges stated
(CHOICE ranges, not measured components), and should be reviewed before they are quoted.

## 1. Result

Phase noise alone leaves the intensity constant: the simulated control gives -376 dB/Hz, the numerical floor.
Intensity noise appears only where an element filters the spectrum. Two elements matter: narrow band-pass
filters and multipath (any delayed copy, from two reflections). Gain tilt, fibre dispersion and polarisation
crosstalk in PM fibre are harmless. The criterion is a parametric heating rate Gamma = pi^2 nu^2 S(2 nu) of
at most 0.1 s^-1 at nu = 84 kHz, i.e. a relative intensity noise (RIN) density at 168 kHz of
S <= 1.4e-12/Hz (-118 dB/Hz). A representative chain gives -134.5 dB/Hz (Gamma = 2.5e-3 s^-1), 16 dB inside the
criterion. The design rules that follow are:
(i) no filter narrower than about 100 pm (20 times the 1.35 GHz rms modulation width) in the light path;
(ii) the product of any two reflections in the chain, R1 R2, below 1.7e-3 (amplitude r <= 0.042), i.e. no pair of
uncoated glass surfaces in series.

## 2. Elements

RIN density in the 100-300 kHz band, record 131 us at 16 GS/s, worse of two relative phases for copies:

| element | case | RIN rms (all frequencies) | RIN at 168 kHz | heating rate |
|---|---|---|---|---|
| amplifier gain tilt | 0.05 / 0.5 / 5 dB/nm | 6e-5 / 6e-4 / 6e-3 | -237 / -216 / -182 dB/Hz | < 1e-6 s^-1 |
| band-pass, FWHM 1 nm, offset 250 pm | ASE filter | 7e-3 | -173 dB/Hz | < 1e-6 s^-1 |
| band-pass, FWHM 100 pm, centred / offset 25 pm | | 1.0e-2 / 7.0e-2 | -130 / -133 dB/Hz | 7e-3 / 3e-3 s^-1 |
| band-pass, FWHM 30 pm, centred / offset 7.5 pm | | 9.6e-2 / 0.21 | **-110 / -113 dB/Hz** | **0.70 / 0.37 s^-1** |
| multipath, R1 = R2 = 1e-2, delay >= 0.3 ns | two -20 dB reflections | 1.4e-2 | -130 dB/Hz | 7e-3 s^-1 |
| multipath, R1 = R2 = 1e-2, delay 0.1 / 0.01 ns | on chip, mm scale | 7e-3 / 1e-4 | -133 / -171 dB/Hz | 3e-3 / 1e-6 s^-1 |
| multipath, R1 = R2 = 1e-3 or 1e-4, any delay | | <= 1.4e-3 | <= -149 dB/Hz | <= 8e-5 s^-1 |
| PM fibre, PER -20 dB, 5 / 20 m (DGD 1.2 ps/m) | then the chip TE guide as polariser | 3e-4 / 5e-3 | -161 / -137 dB/Hz | 6e-6 / 1.4e-3 s^-1 |
| dispersion, 20 ps^2/km, 5 m to 1 km | | <= 8e-4 | < -280 dB/Hz | 0 |
| chain (Sec. 3) | | 6.9e-3 | -134.5 dB/Hz | 2.5e-3 s^-1 |

The structure follows from the physics. Gain tilt and dispersion respond to the instantaneous frequency, the
time derivative of the phase, whose spectrum vanishes towards low frequency. So their intensity noise lies at
GHz, not at 168 kHz. A delayed copy of amplitude r produces an intensity term 2 r cos(phi(t) - phi(t - tau)). Once
tau exceeds the correlation time (about 0.3 ns), this term is a random telegraph-like signal of variance 2 r^2
spread over about B, so S ~ 4 r^2 x 0.40 / B (the 0.40/B is the low-frequency density found in module 15). For
tau far below the correlation time the term is quadratic in the phase increment and much smaller. A band-pass
narrower than the spectrum is the multipath effect again: it is a delay distributed over its impulse response,
about 1/(pi x FWHM) = 40 ps for 30 pm.

**Dependence on B** (r = 1e-2): at tau = 1 ns, S is -125, -130 and -138 dB/Hz for B = 0.5, 2 and 8 GHz, because the
noise spreads as 1/B. At tau = 0.1 ns, B = 0.5 GHz is the quietest (-170 dB/Hz): there the copy is still inside the
correlation time. A higher B therefore helps against long-delay reflections and hurts against short ones. B = 2
GHz is a reasonable middle.

## 3. Representative chain

Amplifier tilt 0.5 dB/nm; ASE filter of 1 nm FWHM, offset 250 pm; PM fibre 5 m with -20 dB polarisation
crosstalk; two -30 dB connector reflections 5 m apart (tau = 25 ns); an on-chip etalon between a -30 dB input
facet and the 2.3% glued emitter facet (module 14), 190 ps round trip; 5 m of dispersion. Result: RIN rms 6.9e-3
(at GHz), -134.5 dB/Hz at 168 kHz, Gamma = 2.5e-3 s^-1 (energy e-folding about 7 min). The on-chip etalon and the
connector pair dominate.

## 4. Comparison with no modulation

Without phase modulation, the same reflections do not produce noise at 168 kHz. Instead they give a slow
intensity offset 2 r cos(Phi) that follows the drift of the path phase Phi: 2% for r = 1e-2, at Hz to kHz rates.
The modulation converts this slow, uncontrolled offset into broadband noise of the same rms. Spread over GHz, it
deposits only 1e-13 to 1e-12 /Hz near the trap frequencies. Phase noise therefore helps against parasitic
etalons as well as against row crosstalk, provided rule (ii) holds.

## 5. Limits and what to measure

The elements are linear and treated separately. Gain dynamics of the amplifier (which can suppress or shape
slow intensity noise) and nonlinear effects (self-phase modulation, stimulated Brillouin scattering) are not
modelled. The criterion uses the radial frequency; the axial frequency is lower and its parametric resonance
lies where the RIN density is the same, so the conclusion holds. Measurement: delayed self-heterodyne or direct
RIN of the modulated, amplified light after the full chain, with the modulation on and off, at 10 kHz to
1 MHz. The prediction to test is a floor below -130 dB/Hz with a 1 nm ASE filter and AR-coated components.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| all RIN densities, heating rates | Sec. 1-4 | SIM | results/m16_pm_am.json (review) |
| element parameter ranges | Sec. 2, 3 | CHOICE | stated ranges, not component data |
| PM-fibre DGD 1.2 ps/m, beta2 20 ps^2/km | Sec. 2 | RECALLED | verify for the chosen fibre |
| Gamma = pi^2 nu^2 S(2 nu) | Sec. 1 | SOURCE | O'Hara et al., arXiv:physics/0003049 (module 3) |
| 0.1 s^-1 criterion | Sec. 1 | CHOICE | energy doubling time ~7 s |
