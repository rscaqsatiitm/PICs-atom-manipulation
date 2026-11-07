# Independent rows without acousto-optic frequency offsets

Follow-up to `docs/stacked_relay.md` (Sec. 5) and `docs/collimator.md` (Sec. 8), which proposed one acousto-optic
frequency shifter per chip. Code: `scripts/m15_coherence.py` (vector focusing in `picsim/opa/stack.py`,
`Focus2D.field_vector`); results `results/m15_*.json`. Values are SIM or DERIVED, valid in the regime stated, and
should be reviewed before they are quoted.

## 1. Result

Per-chip frequency offsets are not needed. Nearly the same decoupling (within a factor of 2 of fully incoherent
rows, far inside the criteria) is obtained with one electro-optic phase
modulator on the seed laser, driven with Gaussian noise of about 3 rad rms and 1 to 2 GHz bandwidth, ahead of the
1:16 split. Each chip's feed is then given a different fibre length, in steps of at least three correlation
lengths (about 0.3 m at 2 GHz, at most 5 m in total). Within one chip all paths agree to better than a
millimetre, so the light stays coherent and route E is unaffected. Between chips the relative phase becomes a
random process with correlation time below 1 ns, so the atoms see the time average.

A fully passive alternative, orthogonal linear polarisation on alternate chips, meets the dynamic criteria. It
leaves a slow depth wander of 0.4% rms and needs a polarisation-insensitive objective.

Perturbation of a tweezer by all other rows within +/-16 um (both sides; rms over drifting relative phases,
pairs in quadrature; 4 um rows, full rows of 41 tweezers):

| option | added hardware | power | static depth / position (rms) | when a tweezer moves 0.25 um: depth / position (rms) |
|---|---|---|---|---|
| all chips coherent, no measure | none | 1 | 2.3% / 9.2 nm | 2.3e-3 / 5.1 nm |
| per-chip AOM offsets | 16 AOMs at 2.6 to 3.4 W | 1 (+ AOM loss) | 1.8e-4 / 0.09 nm | 2e-5 / 0.05 nm |
| **common phase noise + fibre delays** | **1 EOM on the seed, noise source, <= 5 m fibre** | **1** | **4.0e-4 / 0.17 nm** | **4e-5 / 0.10 nm** |
| alternating x / y polarisation | half-wave-plate strips (passive) | 1 | 4.2e-3 / 2.5 nm | 1.1e-3 / 1.6 nm |
| y apodiser a = 0.65 | graded mirror (passive) | 1.84 | 5.2e-3 / 4.1 nm | 1.0e-3 / 1.8 nm |
| polarisation + apodiser a = 0.8 | passive | 1.55 | 2.8e-3 / 1.8 nm | 5.7e-4 / 0.7 nm |

(Corrected: an earlier version scaled the whole worst-phase change by the residual coherence, which also
removed the incoherent part and put this option below the fully incoherent one. The residual coherence now scales
only the interference term. The residual used is the finite-record upper bound of `m15_prbs`; the ensemble value
is smaller, so this row is an upper bound.)

(Criteria used for the rows within a chip: others' depth change < 3e-3 and position change < 5 nm when a tweezer
moves; proposal budget 10 nm rms position noise.) Leaving all chips coherent breaks the position budget: 9.2 nm
rms of slow wander plus 5.1 nm when neighbours move.

## 2. What the frequency offset actually does

Light from chip m reaches row m+1 at about 1e-4 of the peak intensity (module 13). If the two chips are mutually
coherent, this field interferes with the neighbour's own field, and the depth changes by about 2 sqrt(epsilon)
cos(phi), where phi is the relative optical phase of the two feeds. Separate fibres let phi drift with temperature
and acoustics. A frequency offset delta f makes phi advance at 2 pi delta f. When delta f is far above the trap
frequencies, the atom sees only the time average, and the interference term vanishes. Any process that
randomises phi on time scales much shorter than the trap period does the same. A frequency offset is one such
process; a common random phase seen with different delays is another.

## 3. Common phase modulation with fibre delays

If the seed carries a phase phi(t) and chip m receives it after a delay tau_m, the relative phase of chips m and n
is phi(t - tau_m) - phi(t - tau_n). Its mean contribution to the interference is the field coherence
|<exp(i [phi(t) - phi(t - tau)])>| at tau = tau_m - tau_n.

| modulation | residual coherence, delays 2 to 15 bits | with finite modulator bandwidth (10-90% rise 0.3 bit) | at arbitrary (non-integer) delays |
|---|---|---|---|
| binary 0/pi m-sequence, N = 1023 | 9.8e-4 (= 1/N) | 4.1e-2 (carrier 0.12) | as integer |
| same, depth error +5% / +10% | 5.2e-3 / 2.4e-2 | | |
| Zadoff-Chu polyphase, N = 1021 | < 1e-13 | 4.7e-3 | 6.2e-2 |
| **Gaussian phase noise, 3 rad rms** | 4.2e-3 (record-limited; ensemble exp(-9) = 1.2e-4) | inherent | 7.2e-3 |
| Gaussian phase noise, 4 rad rms | 5.6e-3 (record-limited) | inherent | 8.2e-3 |

Binary modulation is fragile. Every 0-to-pi and pi-to-0 transition passes through +pi/2, so a band-limited
modulator leaves a carrier (0.12 at 0.3-bit rise), and the delays must be exact multiples of the bit length.
Polyphase codes have perfect autocorrelation but also need exact delays. Gaussian phase noise of 3 rad rms needs
neither: the carrier is exp(-sigma^2/2), and any delay longer than about three correlation times decorrelates.
The 4e-3 to 8e-3 above is the statistical floor of a 32,768-bit record. What the atom sees is set by the
low-frequency density of the interference term, 0.40/B. For the coherent term of module 13 (2 sqrt(epsilon) =
2.1e-2), this gives S_I(2 nu_trap) = 8.8e-14/Hz at B = 2 GHz. The parametric heating rate
Gamma = pi^2 nu^2 S_I(2 nu) (O'Hara et al., arXiv:physics/0003049, used in module 3) is then 6e-3 s^-1 at
nu = 84 kHz, an energy e-folding time above 2 minutes. At B = 0.5 GHz it is 2.4e-2 s^-1.

**Spectral width.** With 3 rad of noise at correlation rate B, the optical spectrum is Gaussian with rms width
0.67 B (1.35 GHz, 5 pm, at B = 2 GHz). Its effect on the tweezers comes through the diffractive objective's axial
chromatic shift (16 um/nm, DERIVED for f = 17.3 mm with f proportional to 1/lambda) and the lateral scaling x
proportional to lambda:

| B | optical rms width | axial rms smear | peak | lateral rms smear at 60 um |
|---|---|---|---|---|
| 0.5 GHz | 0.34 GHz | 0.02 um | 0.99997 | 0.07 nm |
| 2 GHz | 1.35 GHz | 0.08 um | 0.9998 | 0.29 nm |
| 10 GHz | 6.7 GHz | 0.41 um | 0.995 | 1.4 nm |

The spectrum is symmetric about the carrier, so the mean light shift at the double-magic wavelength is unchanged
to first order. The second-order term needs the curvature of the differential polarisability [VALUE?].

**Coherence inside a chip.** Path differences within a chip and its relay are at most about 1 mm (about 5 ps). With
a 0.5 ns correlation time (B = 2 GHz) and sigma = 3 rad, the coherence loss is
sigma^2 (tau/tau_c)^2 / 2 ~ 5e-4. The least-squares synthesis of route E therefore sees an essentially coherent
source.

**Hardware.** A fibre-coupled LiNbO3 phase modulator on the low-power seed and a Gaussian noise source
amplified to about 1 V_pi rms (3 rad) over 1 to 2 GHz. Extra fibre is 0, 0.3, 0.6 ... 4.5 m on the 16 chip feeds, before
or after the amplifiers. A seed broadened to about 1 GHz also raises the stimulated-Brillouin threshold of the
high-power amplifier and delivery fibres, which carry 3 W each [REF?, RECALLED: SBS suppression by phase
modulation is standard in high-power single-frequency fibre amplifiers]. If the chosen amplifier already
broadens its seed for this reason, the only addition is the fibre delays. Compared with 16 AOMs, this removes
16 RF drivers and the AOM insertion loss (about 15% each [VALUE?]).

## 4. Passive options

**Orthogonal polarisation on alternate chips.** A half-wave-plate strip in front of every other chip, in the
collimated region after the collimators, rotates its light to y. We computed vector focal fields
(Richards-Wolf, aplanatic NA 0.5 objective) and summed the intensity over Ex, Ey and Ez (scalar light shift):

| separation | polarisations | static: depth / position (worst phase) | moving: depth / position |
|---|---|---|---|
| 4 um | orthogonal | 8.4e-4 / 1.2 nm | 1.6e-4 / 0.3 nm |
| 8 um | parallel | 3.3e-3 / 2.1 nm | 1.1e-3 / 1.4 nm |
| 12 um | orthogonal | 5.0e-5 / 0.3 nm | 4.4e-5 / 0.1 nm |
| 16 um | parallel | 2.4e-3 / 0.4 nm | 3.5e-4 / 0.6 nm |

At NA 0.5 the orthogonal neighbours' interference falls by a factor of 26 (the residual comes from the
depolarisation of high-angle rays, mainly via Ez). The rows two apart (8 um) share a polarisation and then
dominate. The moving-tweezer criteria are met, and the position budget too (2.5 nm rms). The static depth
wanders by 0.4% rms with the fibre phases, at Hz to kHz rates, far below twice the trap frequencies.

Three costs follow. (i) The meta-objective must be polarisation-insensitive (propagation-phase meta-atoms with
fourfold symmetry rather than geometric phase). (ii) Where the orthogonal fields overlap, the polarisation becomes
elliptical in proportion to their relative phase, which produces a vector light shift. Our estimate, from the
far-detuned alkali formula for Cs at 1064 nm, is a vector-to-scalar ratio of about 0.09 g_F m_F per unit
ellipticity (RECALLED, verify). With an ellipticity of up to about 2e-2 at 4 um, this is about 2e-3 g_F m_F of the
trap depth: up to tens of kHz for stretched states, zero to first order for clock states (m_F = 0). (iii) Only
two orthogonal states exist, so rows two apart stay coherent.

**Passive y apodiser.** A Gaussian graded-transmission (or graded-reflectivity) filter at the plane of the y
telescope conjugate to the pupil apodises all chips at once. It lowers the 1D crosstalk at 4 um from 8.4e-5 to
3.5e-6 (a = 0.65 pupil radii, beta = 0.7 um). However, it costs a factor 1.84 in power and widens the y spot from
1.35 to 1.73 um. In 2D, the coherent sum of a full row, with the circular pupil, still leaves 4.7e-3 of static
depth change at 4 um. It is not competitive.

## 5. Recommendation

Use common phase modulation with fibre delays in place of per-chip AOMs: one modulator, no power penalty, and
the same row independence as frequency offsets. Alternating polarisation costs nothing to add, and it gives a
second, independent factor of about 26 on the nearest rows if the modulation fails or is switched off.

Next step (i) is simulated in `docs/pm_am.md` (filters narrower than about 100 pm and reflection pairs with R1 R2 >
1.3e-3 must be avoided; a representative chain gives -134.5 dB/Hz). Next steps: (i) measure the delayed self-coherence of the modulated seed after the amplifier (PM-to-AM conversion
by any spectral filtering adds intensity noise at GHz; its density at 2 nu_trap should be checked); (ii) choose
the noise drive (or a deterministic multi-tone drive with the same autocorrelation); (iii) obtain the second-order
differential light shift of the 1064.4350 nm double-magic condition over +/-5 pm.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| option table | Sec. 1 | SIM/DERIVED | results/m15_summary.json from m15_coherent and m15_prbs (review) |
| 2D coherent / incoherent, vector fields | Sec. 1, 4 | SIM | results/m15_coherent.json (review) |
| residual coherence, spectra, chromatic smear | Sec. 3 | SIM | results/m15_prbs.json (review) |
| apodiser scan | Sec. 4 | SIM | results/m15_apod.json (review) |
| heating rates | Sec. 3 | DERIVED | O'Hara et al., arXiv:physics/0003049 (used and checked in module 3), density from m15_prbs |
| 16 um/nm axial chromatic shift | Sec. 3 | DERIVED | f/lambda for a diffractive objective, f = 17.3 mm |
| vector-to-scalar ratio 0.09 g_F m_F | Sec. 4 | RECALLED | far-detuned formula (Grimm et al., arXiv:physics/9902072), not read this session; counter-rotating terms neglected |
| 3 rad rms, B 1-2 GHz, 0.3 m delay steps | Sec. 1, 3 | CHOICE | design choices |

Open placeholders: [REF?] SBS suppression by phase modulation in high-power 1064 nm fibre amplifiers;
[VALUE?] AOM insertion loss at 3 W; [VALUE?] second-order differential light shift of the double-magic condition
over a 5 pm spectrum; vector polarisability of Cs at 1064 nm from a primary source.
