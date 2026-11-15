# The planar surface-emitting tile: an alternative to the chip stack and its collimator

Module 30. Code: `scripts/m30_planar_tile.py`, `scripts/m3_grating_fdtd.py` (extended with a reflector stack under
the BOX and a metal substrate); results `results/m30_coupling.json`, `results/m30_directionality.json`,
`results/m30_imaging.json`, `results/m30_power.json`. Values are SIM or DERIVED, valid in the regime stated, and
should be reviewed before they are quoted.

## 1. Result

In the tile proposed in `docs/scaling_beyond.md` (Sec. 5), each engine lies flat. Its emitters continue as weak
gratings that radiate upwards, and the grating aperture is imaged onto the atoms' y axis with demagnification
D = tile pitch / 4 um = 1500 (6 mm tiles). No fast-axis collimator is needed, and tile placement errors shrink by D.
We simulated the four quantities the proposal left open:
- (i) the grating directionality;
- (ii) the y imaging and its angle tolerance;
- (iii) the coupling between neighbouring long emitters;
- (iv) the power per tweezer.

(i) **Directionality is solved by a mirror under the BOX** (2D FDTD, 30 um weak grating, about 10 deg emission).
- **Si substrate alone:** at best 0.70 at a 2.65 um BOX, and between 0.28 and 0.70 as the BOX varies, so the
  BOX thickness must be controlled to about +/-25 nm.
- **Bragg mirror, two Si/SiO2 pairs:** 0.965 to 0.985 for a 1.8 to 2.2 um BOX, except 0.60 at 2.1 um.
- **Bragg mirror, three pairs:** 0.993 to 0.997, except 0.945 at 2.1 um.
- **Perfect metal mirror:** 1.000, but its absorption is not modelled.
- **Deeper etch.** A 0.1 um etch gives the same directionality (0.72 / 0.986 / 0.998 / 1.000).

(ii) **The grating must be about 4 mm long, not 2 mm.** The atoms see the grating's near field demagnified 1500
times, through a pupil of object-side NA 0.5 / 1500 = 3.3e-4. A near field that fills this pupil like the edge
emitter's facet image has a 1/e^2 radius of about D x 0.97 um = 1.45 mm.
- **2 mm grating:** truncation leaves crosstalk at 4 um of 3e-3 to 7e-3 for every profile, 30 to 90 times the edge
  emitter's 8e-5.
- **4 mm grating, Gaussian apodization w = 1.7 mm:** 1.2e-4 crosstalk at a y FWHM of 1.30 um, the shortest design
  that meets it.
- **Angle tolerance of that design:** 5.2e-5 in sine for peak >= 0.97, 3.0e-5 for 0.99, and 3.6e-5 for crosstalk
  <= 2e-4. This confirms the 30 urad estimate. Object defocus is irrelevant at this NA: 1.6 mm changes the peak by less than 1e-7
  (2e-4 um at the atoms).

(iii) **Neighbouring emitters couple, and this sets the y angle per tweezer.**
- **Mechanism.** On a 4 mm grating the emitters form a coupled array. A tweezer at x_a is close to one Bloch mode,
  whose effective index, and so whose y emission sine, depends on x_a. The swing across the field is
  4 kappa / k0.
- **Uniform array.**
  - At the D1 pitch of 1.5 um the swing is 4.7e-3, 90 times the tolerance.
  - Pitch 2.0 um: 1.6e-4. Pitch 2.3 um: 2.2e-5, which meets the tolerance; 2.2 um is about the limit.
  - A wider pitch narrows the grating-lobe-free field to lambda / (2p): 0.65 of D1 at 2.3 um, so each engine
    serves 35% fewer sites, or needs 1.5 times more ports.
- **Alternating widths.** Widths of 1.0 / 1.4, 0.9 / 1.5 or 0.8 / 1.6 um at 1.5 um pitch cut the swing to
  2.1e-4, 1.4e-4 and 1.0e-4: 22 to 45 times, but still 2 to 4 times over the tolerance. The two emitter types also
  radiate at different angles unless their grating periods differ.
- **Method check.** These come from 8-guide array supermodes. The two-guide extraction fails for mismatched guides
  (Sec. 3); for uniform arrays it agrees with the 8-guide result (4.65e-3 against 4.66e-3).

(iv) **Power per 1.2 mK tweezer:** 46.5 mW with a three-pair Bragg mirror (47.1 mW with two pairs, 66 mW with the
substrate alone), against 57.7 mW for the edge emitter with the reflective collimator. This counts the y axis only;
the x axis is that of the edge emitter, at the pitch D1 uses.

**At the wider pitch the coupling needs, the gain is lost (provisional).** The 35% loss of x field at 2.2 um pitch
may be recoverable by changing the x-relay magnification, since the number of resolvable spots equals the number of
emitters at any pitch. The price is the fill factor (1.2 um emitters: 0.80 at 1.5 um, 0.55 at 2.2 um). The fraction
of one emitter's far-field power inside the grating-lobe-free range |s| <= lambda / (2p) falls from 0.71 at 1.5 um to
0.54 at 2.2 um (0.51 at 2.3 um), a relative efficiency of 0.76 (DERIVED from the femwell emitter mode; the x relay is
not simulated). The tile at 2.2 um would then need about 61 mW per tweezer, no better than the edge emitter.

Hence the tile removes the stack, the collimator and its +/-0.1 um glue tolerance. It is about 20% more efficient
only at the D1 pitch, where coupling rules it out; at the pitch coupling allows it is about equal to the edge
emitter. It brings its own costs:
- a 4 mm apodized grating per emitter, with the emission rate controlled along it;
- a Bragg mirror under the BOX;
- an emitter pitch of about 2.2 um, or a phase-mismatched emitter design that has not been found;
- a 96 mm y-relay object field for 16 tiles.
The coupling problem is the open item. It does not arise for the edge emitter, whose emitters end at the facet.

## 2. Directionality (MEEP)

2D FDTD of the partially etched SiN grating:
- 400 nm film, 50 nm etch, duty 0.5, 30 um long;
- top oxide 2 um;
- 40 px/um;
- normalised to a straight-guide reference.
The original cases are reproduced exactly (0.715 at 0.1 um etch, 3 um BOX).

| reflector | best BOX | directionality (0.05 / 0.1 um etch) | range over the BOX values simulated |
|---|---|---|---|
| Si substrate | 2.65 um | 0.702 / 0.722 | 0.28 to 0.70 (2.0 to 2.8 um, 50 nm steps) |
| Bragg, 2 pairs (Si 75 nm / SiO2 184 nm) | 1.9 um | 0.985 / 0.986 | 0.60 to 0.985 (1.8 to 2.2 um) |
| Bragg, 3 pairs | 1.9 um | 0.997 / 0.998 | 0.945 to 0.997 |
| perfect metal | 1.8 to 2.2 um | 1.000 / 1.000 | 1.000 |

The BOX also sets the emission rate. The 1/e power length at 50 nm etch runs from 115 um to 960 um with the metal
mirror, and from 124 um to 2.7 mm with the Bragg mirror. An apodized 4 mm grating therefore needs both its etch
(about 15 nm for a decay length of order 1.5 mm, scaling as etch^-2 from these data) and the BOX controlled along
its length. Back-reflection into the guide is below 1e-3.

## 3. Coupling of long emitters

Supermodes of 1.2 x 0.4 um SiN emitters at 1064.4 nm (femwell). For a uniform array the two-guide splitting gives
kappa, and the Bloch band spans 4 kappa / k0 = 2 x (splitting):

| pitch (gap) | splitting | coupling length | band swing in sine | field lambda / 2p (fraction of D1) |
|---|---|---|---|---|
| 1.5 um (0.3) | 2.3e-3 | 0.23 mm | 4.7e-3 | 0.355 (1.00) |
| 1.7 um (0.5) | 6.0e-4 | 0.89 mm | 1.2e-3 | 0.313 (0.88) |
| 2.0 um (0.8) | 8.0e-5 | 6.7 mm | 1.6e-4 | 0.266 (0.75) |
| 2.3 um (1.1) | 1.1e-5 | 48 mm | 2.2e-5 | 0.231 (0.65) |
| 2.6 um (1.4) | 1.5e-6 | 0.35 m | 3.1e-6 | 0.205 (0.58) |

For alternating widths the two-guide extraction fails. The splitting at 1.5 um pitch falls below the index
difference of the isolated guides (by 5e-5 to 8e-4), because each guide's index is pulled by its neighbour's
dielectric, and unequally. We therefore read the band width from the spread of the four upper-band supermodes of an
8-guide array. For a chain of m guides, n_j = n0 + c cos(j pi / (m + 1)), and the infinite band spans 2c. The same
method on uniform 8-guide arrays reproduces the two-guide values (4.65e-3 at 1.5 um, 1.60e-4 at 2.0 um).

## 4. Imaging and power

1D scalar Fourier optics:
- near field A(u) over the grating, demagnified by D = 1500;
- object-side pupil |s| <= 3.3e-4;
- normalised so that the line integrates to the fraction of the emitted power that passes the pupil, as in the
  wave model of module 13.

| profile | crosstalk at 4 um | peak per emitted power | y FWHM |
|---|---|---|---|
| 2 mm, uniform weak grating (Le 1 mm) | 3.3e-3 | 0.73 | 1.03 um |
| 2 mm, Gaussian w 1.2 mm | 3.6e-3 | 0.80 | 1.02 um |
| 4 mm, Gaussian w 1.45 mm | 3.9e-4 | 0.78 | 1.22 um |
| 4 mm, Gaussian w 1.7 mm (adopted) | 1.2e-4 | 0.73 | 1.30 um |
| 5 mm, Gaussian w 1.7 mm | 4e-5 | 0.69 | 1.38 um |

Tolerances of the adopted profile:
- **Emission-angle error:** 20 / 40 / 60 / 100 urad give peaks of 0.996 / 0.983 / 0.962 / 0.896 and crosstalk of
  1.4e-4 / 2.2e-4 / 3.9e-4 / 1.1e-3.
- **Decentre:** the tile's y placement moves the row by 1/1500 of the error.

**Power.** The edge chain delivers a line peak of 0.59 per unit guided power, from the facet transmission into
glass, the pupil fill and a silvered M1 with a TIR M2 (m14, m27). The tile delivers 0.733 x directionality. The
power per tweezer scales inversely, which gives the figures of Sec. 1(iv).

**Layout.** 16 tiles at 6 mm make a 96 mm y-relay object field at NA 3.3e-4. Its Lagrange invariant matches the
atoms (48 mm x 3.3e-4 = 16 rows x 4 um / 2 x 0.5). A 4 mm grating takes two thirds of the 6 mm tile pitch along
the grating, so the engine (FPR slab of about 2 mm) has to sit beside the emitter array in x rather than in line.

## 5. Limits

- **Dimensionality.** The FDTD is 2D (infinitely wide emitters) with a uniform 30 um grating. Directionality is
  taken as independent of length and apodization. Metal absorption is not modelled.
- **Coupling.** The coupling is computed for straight guides, not for the etched grating region, which changes the
  modes slightly. The tolerance applies to the swing across the whole field; partial compensation by trimming each
  emitter's grating period was not studied.
- **Imaging.** The imaging model is scalar and 1D in y; the x axis and the circular pupil are not combined here.
- **Engine.** The engine itself (FPR, actuators) is unchanged and not re-simulated at a wider emitter pitch.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| directionality, emission rate, back-reflection | Sec. 1, 2 | SIM | results/m30_directionality.json (2D FDTD; review) |
| coupling, band swing, field fraction | Sec. 1, 3 | SIM | results/m30_coupling.json (review) |
| imaging profiles, angle and focus tolerance | Sec. 1, 4 | SIM | results/m30_imaging.json |
| power per tweezer | Sec. 1, 4 | DERIVED | results/m30_power.json |
| capture fraction vs pitch, relay-magnification alternative | Sec. 1 | DERIVED | results/m30_power.json (x_relay_alternative; provisional) |
| tile pitch 6 mm, D = 1500, grating length choices | Sec. 1 | CHOICE | docs/scaling_beyond.md Sec. 5 |
| n(Si) 3.55, n(SiO2) 1.4496 at 1064.4 nm | Sec. 2 | SOURCE / CHOICE | picsim/params.py, picsim/materials.py |
| edge chain 57.7 mW | Sec. 1 | SIM | results/m27_wave.json |

Open placeholders: [REF?] fabrication of SiN gratings with a Si/SiO2 Bragg mirror under the BOX; [VALUE?] absorption
of a practical metal mirror at 1064 nm under oxide.
