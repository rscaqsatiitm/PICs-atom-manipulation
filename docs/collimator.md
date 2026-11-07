# Fast-axis collimator, facet loss, relay lens and the power plan

> Follow-up (`docs/megascale.md`, Sec. 12): this module and `docs/stacked_relay.md` treat the y axis for light with
> no x component (s_x = 0). A tweezer at x_a is made of light with s_x = x_a / 288 um, and the cylindrical collimator
> then acts with effective indices sqrt(n^2 - s_x^2). The adopted collimator defocuses the facet by -99 s_x^2 um. The
> y figures below therefore hold near the row centre: the peak falls to 0.97 at x_a = 30 um and 0.48 at 60 um.
> With the best common focus offset a row is usable over |x_a| <= 45 um.

Follow-up to `docs/stacked_relay.md`, which assumed an ideal collimator (focal length 175 um, NA 0.9, Abbe sine
condition). Code: `picsim/raytrace.py`, `picsim/opa/collimator.py`, `scripts/m14_collimator.py`,
`scripts/m14_facet_meep.py` (2D FDTD); results `results/m14_*.json`. Values are SIM or DERIVED, valid in the regime
stated. They should be reviewed before they are quoted.

## 1. Result

A collimator that meets the requirement exists in a classical form. It is an aplanatic (Weierstrass) cylinder
on a fused-silica block glued to the facets, followed 10 um later by a two-asphere aplanat of index 1.8. Ray
traced, it collimates to NA 0.9 with the sine condition exact and no residual wavefront error. Rays meet every
surface at 25 deg or less, and the whole lens is 0.40 mm long and 0.32 mm tall. In the wave model it performs
like the ideal lens in glass. Two findings change the system numbers:

(i) **The emitter facet loses about a quarter of the light into air.** A 2D FDTD of the facet cross-section
gives 79% transmission into air (12.8% reflected back into the chip) and 96% into glass (2.3% reflected). Within
the collimator aperture (|s| < 0.9) the figures are 77% (air) and 85% (glass). This loss was not in any earlier
power figure. Gluing the collimator to the facet recovers about 10% and cuts the back-reflection sixfold.

(ii) **Power per 1.2 mK tweezer rises to 57 mW** with the adopted collimator, from 47 mW (`docs/stacked_relay.md`,
facet loss omitted). With an ideal lens in air it would be 61 mW.

| quantity | adopted collimator (Weierstrass block + aplanat) |
|---|---|
| focal length / NA / sine-condition error over the aperture | 175 um / 0.9 / < 4e-5 |
| length (facet to last vertex) / clear height | 395 um / +/-158 um |
| largest incidence angle on any surface | 25 deg |
| peak per guided power vs ideal lens in air | 1.10 (2D check: 1.107) |
| row crosstalk at 4 um (1D line; 2D single tweezer 1.7e-4) | 8.4e-5 |
| decentre / focus for <= 3% loss (1D line; in 2D: 0.988 / 0.953 at 0.5 / 1 um focus, 0.968 at 0.1 um decentre) | +/-0.1 um / +/-1 um |
| surface form error for <= 2% loss | <= 20 nm rms |
| laser power per 1.2 mK tweezer | 57 mW |

The decentre tolerance (+/-0.1 um) does not depend on the lens design. It is set by the small vertical mode of
the facet. The first x lens of the relay, by contrast, needs no special accuracy, because calibration absorbs
its static aberration (Sec. 6). One HI1060-matched input fibre per chip carries the 2.6 to 3.4 W a row needs, at
8 to 11 MW/cm^2 (Sec. 8).

## 2. Why NA 0.9 at 175 um is hard

A collimator must bend rays from sin(theta) = 0.9 (64 deg in air) to the axis within a lens no taller than the
0.5 mm chip pitch. Single refracting surfaces fail:

(i) a plano-hyperbolic lens with its curved side to the facet collimates exactly, but only up to
sin(theta) = sqrt(1 - 1/n^2): 0.72 for fused silica and 0.83 for n = 1.8;
(ii) a plano-asphere with its flat side to the facet reaches 0.9. In fused silica its exit surface slopes at 81 deg
and rays meet it at 43 deg, at the edge of total internal reflection. With n = 1.8 the slope falls to 58 deg, but the
flat surface then sees 64 deg incidence (51% bare transmission), which an AR coating cannot easily recover;
(iii) one elliptical surface on a block glued to the facet (source inside glass, eccentricity 1/n) collimates
exactly at any NA up to sqrt(n^2 - 1) = 1.05. However, it violates the sine condition: the ray at s = 0.9 exits at
68% of the height an aplanatic lens would give.

The Weierstrass construction resolves this. A circular surface (radius R) on a glass block, with the facet R/n
from its centre, is aplanatic: every ray appears to come without aberration from a virtual point, and the NA
drops by n^2 (0.9 -> 0.43 for fused silica). A second lens of focal length f n^2 = 368 um then collimates NA 0.43.
We design it as a two-asphere aplanat in 2D (Wassermann-Wolf construction: exact collimation and exact sine
condition, `collimator.aplanat`). A scan of R, the air gap, the aplanat thickness and its index
(96 combinations, 29 reaching NA 0.9; `m14_scan`) selected R = 80 um, a 10 um gap, a 250 um thick aplanat and n = 1.8. That
combination gives the lowest largest incidence angle (25.4, 18.2 and 19.2 deg on the three surfaces) and slopes
of 64, 7 and 36 deg.

| design | sine condition y_S/(f s) at s = 0.3 / 0.6 / 0.9 | largest incidence | surfaces | Fresnel transmission at the edge, bare |
|---|---|---|---|---|
| Weierstrass + aplanat (glued) | 1.000 / 1.000 / 1.000 | 25 deg | 3 (one spherical) | 0.74 |
| ellipse (glued or printed) | 0.954 / 0.834 / 0.675 | 43 deg | 1 | 0.43 |
| plano-asphere, n = 1.8 (air) | 0.988 / 0.960 / 0.971 | 64 deg | 2 | 0.51 |

## 3. Facet transmission

The x model of the chip treated the emitter facet as lossless, and the cross-axis model of module 13 normalised
to the power radiated into air. A Kirchhoff estimate puts the loss into air at 27%. It comes from 8% Fresnel
reflection plus the 20% of the angular spectrum of the 0.39 um vertical mode that lies beyond the light cone.
We checked this by 2D FDTD (MEEP, resolution 40/um) of the vertical cross-section. The SiN core index was set
to 1.965 so that the slab mode index equals that of the 3D mode (1.7995). The guide ends in air, or in fused
silica for a glued lens.

| into | reflected | transmitted | within the collimator aperture, FDTD | same, Kirchhoff |
|---|---|---|---|---|
| air | 12.8% | 79.1% | 77.2% | 72.3% |
| fused silica | 2.3% | 95.9% | 85.2% | 82.9% |

The two methods agree to 3 to 7%. We use the FDTD values to scale the Kirchhoff spectra in the power estimate.
The remainder (8% in air, 2% in glass) leaves through the cladding near the facet. It is a 2D result with the
lateral confinement folded into the core index; a 3D FDTD of the facet would refine it.

## 4. The designs in the wave model

The wave model takes the facet spectrum in air or glass, each lens's traced map (exit height and wavefront vs s,
Fresnel transmission), the pupil fill beta and the NA 0.5 pupil. Peak values are per unit guided power,
relative to the ideal lens in air at beta = 0.7 um, with AR-coated lens surfaces:

| design | beta | peak | y FWHM | crosstalk at 3 / 4 um |
|---|---|---|---|---|
| ideal, air | 0.7 um | 1.00 | 1.31 um | 2.8e-4 / 1.1e-4 |
| **Weierstrass + aplanat** | 0.7 um | **1.10** | 1.35 um | 4.0e-4 / 8.4e-5 |
| ellipse, printed | 0.9 um | 1.03 | 1.35 um | 6.6e-4 / 1.8e-4 |
| ellipse, printed | 0.7 um | 1.22 | 1.13 um | 2.1e-3 / 1.0e-3 |
| plano-asphere n = 1.8 | 0.7 um | 1.03 | 1.27 um | 4.0e-4 / 1.7e-4 |

The aplanatic design reproduces the ideal lens in glass exactly. The ellipse narrows the y spot because its
mapping compresses the outer far field, but the same distortion raises the crosstalk. Opening the pupil fill
to beta = 0.9 um restores 1.8e-4 at 4 um at a cost of 15% in peak. Without AR coatings, all designs lose 10 to
30% (n = 1.8 surfaces reflect 8% each at normal incidence), so coatings are needed in every case.

## 5. Tolerances with the real lenses

| error | ideal (air) | Weierstrass + aplanat | ellipse | plano-asphere n 1.8 |
|---|---|---|---|---|
| decentre 0.1 / 0.2 um | 0.973 / 0.897 | 0.973 / 0.896 | 0.962 / 0.855 | 0.971 / 0.887 |
| focus 0.5 / 1 um | 0.989 / 0.937 | 0.995 / 0.977 | 0.974 / 0.905 | 0.978 / 0.906 |
| form error 20 / 50 nm rms (corr. 20 um), mean | | 0.984 / 0.916 | 0.998 / 0.988 | 0.989 / 0.916 |
| form error 50 nm, 10th percentile | | 0.880 | 0.981 | 0.865 |

(Peak relative to no error, 1D y line, beta = 0.7 um for all designs.)

The tolerance model was checked against exact traces of a displaced source. For decentres of 0.2 and 0.5 um
and a 0.5 um focus error, the exactly traced path change differs from the first-order model by at most 0.003
waves rms (0.03 waves peak-to-valley for the ellipse at 0.5 um, about 2e-4 waves for the aplanatic design).

Three points follow. (i) The decentre tolerance is a property of the facet, not of the lens: +/-0.1 um for every
design. (ii) Focus is looser for the glued designs because the facet radiates into glass (the defocus phase
scales as n_m cos theta): about +/-1 um for 3%. (iii) The ellipse has one surface with a small index step
(n - 1 = 0.45), so it tolerates 50 nm rms form error; the n = 1.8 surfaces need about 20 nm.

## 6. The first x lens of the relay

The x relay's first lens (f = 3.12 mm, NA 0.355 through the slit) sees the whole emitter row, +/-144 um, so it
works off axis. We traced two fused-silica designs from 7 emitter heights. Each off-axis wavefront was applied
to each emitter's angular spectrum in the D1 model (BPM matrix), and tweezers at 0, +/-30 and +/-60 um were
evaluated before and after least-squares recalibration. (Corrected: the first aplanat, with r0 = 1.5 mm, covered
only |s| <= 0.33, short of the slit edge 0.355, and the wavefronts were quoted at the 125 um field point; the
aplanat now uses r0 = 2.0 mm, which covers the slit at every field point, and both lenses are evaluated at
144 um. Its back focal distance is 6.76 mm.)

| lens | wavefront at 144 um field (tilt removed) | edge tweezer, uncalibrated: peak / shift | after recalibration: peak per power / shift | moving tweezer, others' change |
|---|---|---|---|---|
| aplanat (sine condition), r0 2.0 mm | 0.12 waves rms | 0.983 / 21 nm | 0.998 / <= 1.3 nm | 4e-4 / < 1 nm |
| plano-asphere (no sine condition) | 0.76 waves rms | 0.672 / 590 nm | 0.964 / <= 0.5 nm | 4e-4 / 0.7 nm |

The static aberration of the x relay is a fixed linear transformation of the chip's response. The
least-squares synthesis, which already works on the measured response matrix, removes it. Hence the x relay
needs ordinary lens quality: an aplanat costs nothing, and even a lens without the sine condition costs only
3.6%. The same argument does not apply to the collimator, because the y axis has no actuators.

## 7. Assembly and drift

(DERIVED, with RECALLED coefficients marked.)

**Thermal.** For a collimator glued directly onto its own facet, vertical registration is fixed at the facet.
Differential expansion then acts only over the lens's own height (+/-80 um at the glued block): with
Delta alpha ~ 2e-6/K (Si 2.6e-6/K vs fused silica 0.55e-6/K, RECALLED), that is about 0.2 nm/K, negligible. A
monolithic collimator array for the whole stack, referenced at its base, would see the full 8 mm stack height:
16 nm/K at the top chip, so 0.1 um after 6 K. Hence one collimator per chip, or a stack held to +/-1 K.

**Adhesive.** The glue line of the immersed design lies between the facet and the block, along z. Shrinkage
during cure therefore changes the focus, which tolerates about +/-1 um, and only weakly the decentre. The vertical
position must be held during cure to +/-0.1 um. Typical shrinkage and creep of UV adhesives need checking
[REF?].

**Alignment signal.** Decentre moves the facet image across the pupil by delta x 12,400 (`docs/stacked_relay.md`,
Sec. 2): 0.1 um becomes 1.2 mm on the 17 mm pupil. Imaging the pupil, or its conjugate in the y telescope, while
each chip is lit gives a direct and very sensitive alignment signal. A camera with 10 um pixels resolves 1 nm of
decentre in principle. Active alignment against this signal, before and during cure, is the recommended
process. A collimator printed by two-photon lithography directly on the facet would be registered
lithographically; its placement accuracy needs a source [REF?].

## 8. Power and frequency plan for 16 chips

Laser power per 1.2 mK tweezer, field mean, per unit guided power (FDTD-scaled facet, AR-coated lenses):

| collimator | power per tweezer |
|---|---|
| ideal lens in air | 61 mW |
| Weierstrass + aplanat, glued (adopted) | 57 mW |
| ellipse, printed, beta 0.9 | 61 mW |
| plano-asphere n = 1.8 in air | 59 mW |

| tweezers per row | depth | power per chip | 16 chips | input fibres per chip: lensed tip / HI1060-matched tip |
|---|---|---|---|---|
| 45 (4 um) | 1.2 mK | 2.6 W | 41 W | 5 / 1 (8 MW/cm^2) |
| 60 (3 um) | 1.2 mK | 3.4 W | 55 W | 6 / 1 (11 MW/cm^2) |
| 45 | 0.5 mK | 1.1 W | 17 W | 2 / 1 |
| 60 | 0.5 mK | 1.4 W | 23 W | 3 / 1 |

**One input per chip.** The lensed-fibre inverse taper of module 2 has a 3.1 um^2 mode at the chip facet, so a
20 MW/cm^2 cap (CHOICE; no damage threshold was found) allows 0.62 W per fibre and 5 to 6 fibres per chip. Several
inputs would split the 128 ports into sub-trees whose relative optical phases drift with their fibres, and a
tweezer window straddling a sub-tree boundary would then fluctuate. The HI1060-matched taper of module 2 (tip
0.1 um, -0.48 dB, 0.11 dB worse than the lensed tip) has a 30.6 um^2 mode and carries 6.1 W at the same cap.
Hence one fibre per chip is possible and preferred. Its 0.1 um tip width is at the assumed minimum feature
size.

**Frequency offsets.** Superseded by `docs/row_coherence.md`: per-chip shifters are not needed. One phase
modulator on the seed, driven with Gaussian noise (3 rad rms, 1 to 2 GHz), plus fibre delays in 0.3 m steps gives
the same row independence. (The original plan was one acousto-optic shifter per chip feed at 80 + 1.13 m MHz,
m = 0 to 15, pairwise differences 1.13 to 17 MHz, far above 2 nu_trap = 168 kHz, each handling 2.6 to 3.4 W
[VALUE?].)

**Laser.** 41 to 55 W of single-frequency 1064.4 nm light at 1.2 mK, or 17 to 23 W at 0.5 mK. Lowering the transport
depth remains the cheapest lever: power scales linearly with depth, and the lower trap frequency also relaxes
the required update rate (`m3_update_rate`).

## 9. Limits and next steps

The ray trace is 2D (cylindrical lenses), with ideal AR coatings in the wave model. Polarisation effects at the
steep Weierstrass surface (64 deg slope, 25 deg incidence) are limited to s-polarised Fresnel transmission. The
facet FDTD is 2D with an effective core index. The n = 1.8 glass is a class, not a selected catalogue glass.
Tolerances are evaluated one at a time.

Next steps: (i) a tolerance Monte Carlo of the adopted collimator with all errors combined (decentre, focus,
form, thickness, index); (ii) a 3D FDTD of the emitter facet, glued and in air; (iii) a fabrication route for
the Weierstrass block and aplanat at 0.5 mm pitch (moulded glass array, wafer-level optics or two-photon
printing); (iv) alignment trials using the pupil-walk signal.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| designs, limits, geometry | Sec. 1, 2 | SIM | results/m14_designs.json (review) |
| facet R, T, capture | Sec. 1, 3 | SIM | results/m14_facet_meep_{air,glass}.json (review) |
| wave-model table | Sec. 1, 4 | SIM | results/m14_wave.json (review) |
| tolerances, exact checks | Sec. 1, 5 | SIM | results/m14_tolerance.json (review) |
| form Monte Carlo | Sec. 5 | SIM | results/m14_form.json (24 samples per point; review) |
| x lens | Sec. 6 | SIM | results/m14_xlens.json (review) |
| power plan, fibres, offsets | Sec. 1, 8 | DERIVED | results/m14_plan.json |
| parameter scan of the Weierstrass design | Sec. 2 | SIM | results/m14_scan.json |
| 2D check of the adopted design (1.107) | Sec. 1 | SIM | run in session, not stored |
| n(SiO2) 1.4496 | throughout | SOURCE | Malitson (picsim/materials.py) |
| n = 1.8, R 80 um, gap 10 um, T 250 um, AR coatings | Sec. 2, 4 | CHOICE | design choices |
| CTE Si 2.6e-6/K, fused silica 0.55e-6/K | Sec. 7 | RECALLED | verify |
| 20 MW/cm^2 facet cap | Sec. 8 | CHOICE | m6 (no damage threshold found) |

Open placeholders: [REF?] UV-adhesive shrinkage and creep at the 0.1 um level; [REF?] placement accuracy of
two-photon-printed micro-optics on chip facets; [VALUE?] CW power handling of acousto-optic shifters at 3 W,
1064 nm; [VALUE?] available single-frequency 1064 nm power at 40 to 55 W.
