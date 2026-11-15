# Stacked-chip relay: a 2D array of independently steered rows

> Follow-up (`docs/megascale.md`, Sec. 12): this module and `docs/collimator.md` treat the y axis for light with
> no x component (s_x = 0). A tweezer at x_a is made of light with s_x = x_a / 288 um, and the cylindrical collimator
> then acts with effective indices sqrt(n^2 - s_x^2). The adopted collimator defocuses the facet by -99 s_x^2 um. The
> y figures below therefore hold near the row centre: the peak falls to 0.97 at x_a = 30 um and 0.48 at 60 um.
> With the best common focus offset a row is usable over |x_a| <= 45 um. With the all-reflective collimator of
> `docs/collimator.md`, Sec. 10, the s_x = 0 figures hold across the whole row (1D peak >= 0.98 to the slit edge).

Follow-up to `docs/scaling.md`, Sec. 4(v). Code: `picsim/opa/stack.py`, `scripts/m13_stack.py`; results
`results/m13_*.json` (+ .png); tests `tests/test_stack.py`. Each chip is design D1 (`docs/fpr_design.md`): 128
input ports, 192 arms and emitters, one row of independently steered tweezers along x. The model is scalar wave
optics. The x axis uses the beam-propagation response of D1. The y axis uses the femwell facet mode of the D1
emitters. Both axes pass through the circular NA 0.5 pupil. Values are SIM or DERIVED and valid in the regime
stated. They should be reviewed before they are quoted.

## 1. Result

Sixteen D1 chips stacked at 0.5 mm, each with its own fast-axis collimator, can be imaged onto the atoms as 16
rows 4 um apart. The rows then behave as independent as the tweezers within a row, on one condition: every
chip runs at its own laser frequency. Design point, 1064.4 nm, NA 0.5:

| quantity | value |
|---|---|
| demagnification of the stack in y, D = chip pitch / row pitch | 125 |
| collimator (FAC) per chip: focal length / numerical aperture | 175 um / 0.9 |
| y spot FWHM at the atoms (x: 1.71 um) | 1.31 um |
| cross-axis transmission (FAC capture x pupil) / loss to the circular pupil | 0.973 / 3% |
| laser power per 1.2 mK tweezer, field mean (was 43 mW with an assumed Gaussian cross axis) | 47 mW (57 mW with the facet loss, `docs/collimator.md`) |
| crosstalk of a full row (41 tweezers) at the neighbouring row, 4 um / 3 um | 1.8e-4 / 4.3e-4 |
| neighbouring row, chips at different frequencies: depth / position change when a tweezer moves | 2.4e-5 / 0.04 nm |
| the same with all chips at one frequency (coherent) | 2.8e-3 / 4.5 nm (static: 2.6%, 11 nm) |
| FAC alignment for <= 3% loss: decentre / focus / roll | +/-0.1 um / +/-0.5 um / 1 mrad |
| FAC wavefront error for <= 6% loss | <= 0.05 waves rms |
| row position at the atoms per um of chip height error | 8 nm |
| facet-to-pupil length of the first-order layout | 0.38 m |

Hence the stacked relay is optically sound. Its costs are (i) one frequency shifter per chip, (ii)
collimator alignment to about 0.1 um in height, which is set by the small vertical facet mode, and (iii) a
diffraction-limited collimator at NA 0.9. With 192-port chips, 16 rows give 720 (4 um spacing) to 960 (3 um)
independent tweezers. Including the emitter facet loss found later (`docs/collimator.md`), they need 41 to 55 W
at 1.2 mK, or 17 to 23 W at 0.5 mK, spread over 16 chips.

## 2. Optical layout

In x, the relay of `picsim/opa/array.py` images the emitter plane onto the objective pupil, so the atoms see
the far field of the chip. In y, each chip's facet sits in the front focal plane of its own cylindrical
collimator. The collimators' back focal planes form a plane S, which a common y lens group of effective focal
length F_y Fourier transforms onto the pupil. Hence S is imaged onto the atoms with magnification
f_obj/F_y = 1/D, and chip m lands at y_m = m d_c/D. Here, d_c is the chip pitch and f_obj = 17.3 mm the objective
focal length (pupil radius WD tan(asin NA) = f_obj NA at WD = 15 mm). At the same time, the pupil is an image of
the facet in both axes: in x with magnification 60, and in y with 12,400. Because each chip has its own
collimator centred on its facet, all 16 facet images overlap in the pupil, and their row offsets come only from
the tilt y_m/F_y. A single common collimator would not work: it would place the 16 facet images side by side.

The free parameter of the y design is beta, the pupil radius referred back to the facet, beta = NA f_FAC/D. It
plays the role of the pupil fill (Sec. 4). The row pitch then fixes D, and beta fixes f_FAC.

First-order layout (thin cylindrical lenses; ABCD matrices of both axes checked in `m13_design`):

| row pitch | D | F_y | f_FAC | collimated beam half-height | length facet to pupil |
|---|---|---|---|---|---|
| 3 um | 167 | 2.89 m | 233 um | 210 um | 439 mm |
| 4 um | 125 | 2.17 m | 175 um | 158 um | 382 mm |
| 5 um | 100 | 1.73 m | 140 um | 126 um | 347 mm |
| 6 um | 83 | 1.44 m | 117 um | 105 um | 324 mm |

At 4 um, the long effective focal length is realised as a 25x demagnifying y telescope (f = 100 mm and 4 mm)
followed by an 86.6 mm Fourier lens. The x relay is a 4f pair (3.12 mm and 187.7 mm), with the Fourier slit
6.2 mm from the facets. Both pupils coincide 382 mm from the facets. The Lagrange invariant is the same on
both sides (stack half-height 4 mm x collimated NA 0.0022 = atom half-field 30 um x 0.5 x fill); this holds by
construction of the layout, so it is a consistency check, not an independent test. The collimated
beams (158 um half-height) fit within the 250 um half-pitch, also at 3 um row pitch (210 um).

## 3. The facet sets the y optics

The D1 emitter facet (1.2 x 0.4 um SiN in oxide) has a vertical mode of 1/e^2 equivalent radius 0.39 um. It
radiates into air with a far-field 1/e^2 half-angle of sin(theta) = 0.86 (`m13_facet`). A collimator of NA 0.8
captures 96% of the light, and one of NA 0.9 captures 98.8%. Thinning the facet section does not help much:

| SiN at the facet | far-field NA (1/e^2) | equivalent radius | capture at NA 0.8 | coupling length of adjacent emitters (0.3 um gap) |
|---|---|---|---|---|
| 400 nm | 0.86 | 0.39 um | 0.961 | 228 um |
| 200 nm | 0.79 | 0.43 um | 0.971 | 106 um |
| 150 nm | 0.69 | 0.49 um | 0.980 | 72 um |
| 100 nm | 0.51 | 0.67 um | 0.991 | 47 um |

A thinner facet expands the vertical mode, but it also expands the lateral one. At 1.5 um emitter pitch the
coupling length drops to 47 um, which leaves no room for an adiabatic vertical taper. Hence the dense emitter
array and a large vertical facet mode are in conflict. The collimator alignment tolerances below scale as the
vertical mode radius (decentre) and its square (focus). We therefore keep the 400 nm facet and accept the
tolerances.

## 4. Pupil fill and the crosstalk floor

The y profile at the atoms is the facet's far field, clipped by the collimator aperture and scaled by
beta/NA, then low-pass filtered by the pupil. Crosstalk between rows has two sources:
(i) the pupil clips the sharply peaked image of the facet (small beta);
(ii) the collimator aperture cuts the wide far field with a hard edge, and the pupil makes this edge ring (large beta).

FAC NA 0.9, single tweezer, 1D (y only):

| beta | cross-axis transmission | y FWHM | power vs Gaussian w0 = 1.04 um | crosstalk at 3 / 4 / 6 / 8 um |
|---|---|---|---|---|
| 0.5 um | 0.946 | 1.09 um | 0.92 | 3.5e-3 / 1.8e-3 / 7.2e-4 / 3.9e-4 |
| 0.6 um | 0.968 | 1.18 um | 0.97 | 1.2e-3 / 6.2e-4 / 2.5e-4 / 1.3e-4 |
| **0.7 um** | **0.973** | **1.31 um** | **1.07** | **2.8e-4 / 1.1e-4 / 3.2e-5 / 1.5e-5** |
| 0.8 um | 0.973 | 1.49 um | 1.22 | 1.1e-3 / 2.6e-4 / 4.9e-5 / 1.9e-5 |
| 1.0 um | 0.977 | 2.08 um | 1.65 | 2.1e-3 / 6.4e-4 / 1.9e-4 / 9.3e-5 |

The two sources balance near beta = 0.7 um, which we adopt. Raising the collimator NA to 0.95 lowers the
crosstalk at 4 um to 7.4e-5. A graded (raised-cosine) collimator edge removes source (ii), but it costs 7 to 8%
of the light and reaches only 9.4e-5 at 4 um (beta = 0.8 um). None of these options gets below about 5e-5 at
4 um. That is adequate for mutually incoherent rows (Sec. 5) but not for coherent ones. At beta = 0.7 um the y
spot (1.31 um FWHM) is wider than the Gaussian assumed so far (w0 = 1.04 um, 1.22 um FWHM). Together with
the 2.7% cross-axis loss and the 3% loss to the circular pupil, the power per tweezer rises by a factor 1.10,
the same at the field centre and edge.

## 5. Crosstalk between rows

For the 2D test, chip m forms a full row: 41 tweezers at 3 um over +/-60 um, alternating phases 0 and pi/2. The
neighbouring chip forms 11 tweezers offset by 1.5 um in x. We evaluate them at row separations of 3 to 16 um.
"Coherent" means both chips use the same laser frequency, with the worst of 8 relative phases. "Incoherent"
means the chips use different frequencies.

| row separation | full-row crosstalk (max over x) | incoherent: static depth / position | incoherent, row-m tweezer moves 0.25 um | coherent: static depth / position | coherent, moving |
|---|---|---|---|---|---|
| 3 um | 4.3e-4 | 2.0e-4 / 0.13 nm | 6.2e-5 / 0.04 nm | 2.9% / 18 nm | 6.5e-3 / 5.6 nm |
| 4 um | 1.8e-4 | 1.7e-4 / 0.002 nm | 2.4e-5 / 0.04 nm | 2.6% / 11 nm | 2.8e-3 / 4.5 nm |
| 6 um | 6.3e-5 | 4e-6 / 0.02 nm | 3e-6 / 0.01 nm | 0.41% / 8.9 nm | 1.7e-3 / 2.1 nm |
| 8 um | 4.3e-5 | < 1e-6 | < 1e-6 | 0.08% / 1.8 nm | 1.1e-3 / 1.7 nm |
| 12 um | 1.5e-5 | < 1e-6 | < 1e-6 | 0.01% / 5.9 nm | 8e-4 / 0.8 nm |
| 16 um | 1.3e-5 | 2e-6 | < 1e-6 | 0.27% / 1.6 nm | 4e-4 / 0.6 nm |

With incoherent rows, a moving tweezer changes the next row by at most 6e-5 in depth and 0.04 nm in position,
even at 3 um row pitch. This is two orders of magnitude below the in-row figures of D1 (0.11%, 3.5 nm). With
coherent rows, the crosstalk field interferes with the neighbour's field, so the depth changes by about
2 sqrt(epsilon). This interference persists at all distances simulated, because a full row has a tail plateau
of about 1e-5 out to 16 um. Its strength follows the relative optical phase of the two chips, which drifts with
the fibre paths, so the 2.6% and 11 nm at 4 um become slow noise. Hence alternating two frequencies is not
enough. Every chip needs its own frequency, with all pairwise differences of at least 1 MHz, far above twice
the 84 kHz trap frequency. That is 16 offsets (for example one acousto-optic shifter per chip feed) for 16
chips. The beat notes then modulate the depth at MHz frequencies with amplitude about 2 sqrt(epsilon) (3% at
4 um). This lies far from the parametric resonance at 2 nu_trap, and the atoms see the time average. The
tweezers within one row stay coherent, as route E requires.

Follow-up (`docs/row_coherence.md`): the chips need not run at different frequencies. Any randomisation of
their relative phase much faster than the trap period suffices. One phase modulator on the common seed, driven
with Gaussian noise, plus different fibre delays per chip does this nearly as well as 16 frequency shifters (within a
factor of 2, far inside the criteria).
Alternating polarisation is a passive partial alternative.

## 6. Alignment tolerances

Per-chip errors, single tweezer at x = 0 (the field edge, x = 60 um, gives the same values within 0.1%). The
table gives the peak at the nominal atom plane relative to no error, and the worst crosstalk at 4 um. None of
these errors moves the row in y by more than 0.01 nm.

| error | value | peak | crosstalk at 4 um | note |
|---|---|---|---|---|
| none | | 1 | 1.1e-4 | y FWHM 1.33 um (2D) |
| FAC decentre (facet vs collimator axis, y) | 0.05 / 0.1 / 0.2 / 0.3 um | 0.992 / 0.969 / 0.881 / 0.750 | 2e-4 / 4e-4 / 1.2e-3 / 2.8e-3 | moves the facet image across the pupil (delta/beta) |
| FAC focus (along z) | 0.25 / 0.5 / 1 / 2 um | 0.993 / 0.973 / 0.893 / 0.641 | 3e-4 / 9e-4 / 2.8e-3 / 4.3e-3 | y focus moves by about 1 um at 1 um error |
| FAC roll about z | 0.5 / 1 / 2 mrad | 0.998 / 0.992 / 0.968 | 1e-4 / 1e-4 / 2e-4 | decentre varies along the 288 um emitter row |
| smile (bow, peak-to-valley) | 0.1 / 0.2 / 0.4 um | 0.998 / 0.992 / 0.967 | 1e-4 / 2e-4 / 4e-4 | |
| FAC wavefront, spherical-like | 0.05 / 0.1 waves rms | 0.942 / 0.783 | 3e-4 / 1.2e-3 | best focus 0.95 / 0.82 |
| FAC wavefront, coma-like | 0.05 / 0.1 waves rms | 0.937 / 0.767 | 5e-4 / 2.0e-3 | |
| chip piston along z (x axis, D1 weights unchanged) | 10 / 100 / 300 um | >= 0.999 / >= 0.992 / >= 0.908 | | x position <= 1.5 / 0.9 / 6.8 nm, static |

The collimator decentre is the tightest tolerance. A decentre delta shifts the facet image in the pupil by
delta/beta, so the tolerance scales with the vertical facet mode (Sec. 3). It costs power and crosstalk but
does not move the row. The wavefront requirement (about 0.03 waves rms for <= 3%, at NA 0.9) follows from the
same imaging: the collimator is the aperture of the facet-to-pupil image, and light it scatters out of that
image misses the pupil. Defocus-like aberration is harmless after refocusing. A larger beta would relax the
wavefront tolerance, but it raises the crosstalk and the power (Sec. 4).

Row position at the atoms (DERIVED from the conjugations): a height error of chip plus collimator moves its
row by 1/D = 8 nm per um. A collimator-array placement error has the same effect. A chip tilt about x of 1 mrad
moves the row by f_FAC x 1 mrad / D = 1.4 nm. A chip piston of dz defocuses its row in y by dz/D^2 (6 nm per
100 um). A stack built to +/-2 um therefore places rows to +/-16 nm. A static offset of this size must be
matched by the static 1075 nm array, or measured and absorbed into its design. Thermal expansion of a Si stack
(8 mm) moves the outer rows by about 0.2 nm/K at the atoms (DERIVED, CTE 2.6e-6/K, RECALLED).

## 7. Limits and next steps

The model is scalar and uses ideal (sine-condition) thin cylindrical lenses with the error types listed. The
objective's apodisation factor and polarisation are omitted, as in the x model. The facet radiates into air with
Kirchhoff boundary conditions; Fresnel reflection and the angle dependence of transmission are not included. We
have not designed real lens prescriptions: the 0.38 m layout is first order. Gluing the collimator to the facets
with index matching would lower the ray angles inside the lens for the same NA 0.9; it is not simulated. Fibre
and electrical access to chips stacked at 0.5 mm pitch is not addressed. Facets must stay flush to within about
100 um (Sec. 6), so stair-step stacking for wire bonds is excluded, and access must come from the chip sides or
back.

Next steps (i) to (iii) below are addressed in `docs/collimator.md` (ray-traced collimators, assembly and drift,
power and frequency plan) together with the first x lens of the relay. Remaining: a full lens design of the
anamorphic relay over its 8 mm field.

Original list: (i) ray-trace a real collimator (acylinder or two-photon-printed lens) to NA 0.9 and check the
0.03-wave target; (ii) an assembly study of the 0.1 um collimator alignment, including adhesive shrinkage;
(iii) the frequency plan and power distribution for 16 chips (16 shifters, about 2 to 3 W each); (iv) a real lens
design of the anamorphic relay, including its field (16 chips over 8 mm).

## Verification notes

No literature values were retrieved in this module. Recalled values are labelled.

| value | location | type | pointer |
|---|---|---|---|
| facet NA, capture, coupling lengths | Sec. 3 | SIM | results/m13_facet.json (review before quoting) |
| fill table, graded edges | Sec. 4 | SIM | results/m13_fill.json (review) |
| 2D crosstalk, coherent and incoherent | Sec. 1, 5 | SIM | results/m13_crosstalk.json (review) |
| tolerances | Sec. 1, 6 | SIM | results/m13_tolerance.json (review) |
| power factor 1.10, 47 mW, circle loss 3% | Sec. 1, 4 | SIM/DERIVED | results/m13_power.json; 43 mW from docs/fpr_design.md |
| layout table, ABCD checks | Sec. 2 | DERIVED | results/m13_design.json |
| row position per um, tilt, piston | Sec. 6 | DERIVED | conjugations of Sec. 2 |
| chip pitch 0.5 mm, 16 chips, beta 0.7 um, FAC NA 0.9, telescope 25x | Sec. 1, 2 | CHOICE | stated design choices |
| f_obj 17.3 mm | Sec. 2 | DERIVED | WD 15 mm (USER) and NA 0.5 |
| Si CTE 2.6e-6/K | Sec. 6 | RECALLED | verify |
| >= 1 MHz frequency spacing | Sec. 5 | DERIVED | 2 x 84 kHz trap frequency (m3_update_rate) |

Open placeholders: [REF?] achievable wavefront error of NA 0.9 micro-collimators (commercial acylinders,
two-photon-printed lenses); [REF?] alignment accuracy and long-term stability of actively aligned, bonded
collimators; [VALUE?] power handling of fibre-coupled or free-space frequency shifters at 2 to 3 W, 1064 nm.
