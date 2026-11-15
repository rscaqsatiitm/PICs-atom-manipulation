# Addressing at 894.6 nm: the PIC router against crossed AODs and an SLM/DMD modulator

Module 31. Code: `scripts/m31_addressing.py`; results `results/m31_patterns.json`, `results/m31_leakage.json`,
`results/m31_frequency.json`, `results/m31_patterns.png`. Device figures for the PIC come from modules 4 and 5
(`docs/findings.md`, Sec. 5). Values are SIM, DERIVED or SOURCE as marked, and should be reviewed before they are
quoted.

## 1. Result

Addressing is where the PIC route has a structural lead over crossed AODs. It is not alone there: an SLM/DMD
modulator shares most of the lead.

(i) **Arbitrary patterns.** Crossed AODs light products of a row set and a column set, so an arbitrary on/off
pattern must be split into such rectangles, one round each. For random patterns this means close to one round per
row:
- 7.8 to 16 rounds on a 16 x 16 block (densities 0.05 to 0.5);
- 31 to 40 on 40 x 40.
With 1 us pulses and 4.5 to 7.5 us AOD access times this is 84 to 340 us per pattern. A router with one switch per
site, or an SLM/DMD modulator, needs one round: 1.2 us for the PIC (26 ns rise plus 0.19 us of shaped settling) and
24 us for a 43 kHz DMD.

(ii) **One optical frequency.**
- **AOD.** A crossed AOD tags each site with its own RF offset, 0.13 MHz per site for a 20 MHz, 150-spot
  deflector. For a light-shift (Z) gate this is harmless. For coherent two-photon addressing against a global beam
  the relative phase walks by 0.84 rad per site separation per microsecond, so it must be tracked per site; 40 sites
  apart, 1 us costs 34 rad.
- **PIC router and DMD.** Both keep one frequency for all sites.

(iii) **Idle-site light.**
- **PIC, measured.** The measured AlN-on-SiN switch reaches 71.4 dB mean extinction and -68 dB on-chip crosstalk.
  After free-space delivery the crosstalk is -50.8 dB, so the delivery optics, not the switch, set the floor
  (`docs/findings.md`, Sec. 2(iv)).
- **PIC, simulated.** The cascaded switch of module 4 exceeds 88 dB (5th percentile), capped at 70 dB by the assumed
  stray-light floor.
- **SLM/DMD.** The demonstrated modulator reaches 46 dB contrast and -44 dB nearest-neighbour crosstalk.
- **AOD.** With k tones on an equally spaced grid, every third-order product 2 f_i - f_j lands on a site of the same
  line. With 8 tones on a 40-site line, 17 idle sites (54%) receive such light at the intermodulation level, for
  which we found no sourced figure [VALUE?].

Hence the defensible "better than conventional" claim is about addressing, not transport:
- arbitrary patterns in one round;
- microsecond switching;
- one frequency;
- extinction above 70 dB at the chip.

Its competitor is the SLM/DMD modulator rather than the AOD. It matches the single round and the single frequency at
more than 10,000 sites, but is about 20 times slower (43 kHz). Its contrast is 46 dB against 71 dB extinction at
the PIC chip, but after free-space delivery the nearest-neighbour crosstalk differs by only 7 dB (-44 against
-51 dB).

## 2. Arbitrary patterns

The number of crossed-AOD rounds is the size of an exact partition of the pattern's ones into combinatorial
rectangles (each site lit once). We bound it below by the real rank of the 0/1 matrix and above by the best of
three: a greedy partition, the number of distinct non-empty rows, and the number of distinct non-empty columns. For
random patterns the bounds nearly coincide:

| block | density 0.05 | 0.10 | 0.25 | 0.50 |
|---|---|---|---|---|
| 4 x 4 | 0.8 | 1.1 | 2.4 | 3.2 |
| 8 x 8 | 2.5 | 3.9 | 6.3 | 7.5 to 7.8 |
| 16 x 16 | 7.8 | 11.4 | 15.3 to 15.6 | 16.0 |
| 24 x 24 | 13.4 to 13.5 | 20.0 to 20.3 | 24.0 | 24.0 |
| 40 x 40 | 30.9 to 31.3 | 38.7 to 38.8 | 40.0 | 40.0 |

Time per arbitrary pattern on a 16 x 16 block at density 0.25:

| pulse | crossed AOD | PIC router | SLM/DMD |
|---|---|---|---|
| 0.1 us | 70 to 118 us | 0.32 us | 23 us |
| 1 us | 84 to 132 us | 1.22 us | 24 us |
| 10 us | 222 to 272 us | 10.2 us | 33 us |

Structured patterns (whole rows, checkerboards) need far fewer AOD rounds; these numbers are for independent
sites.

## 3. Leakage, frequency, channels

**Third-order products** on one AOD line (fraction of idle sites that receive them):

| sites per line | 2 tones | 4 | 8 | 16 |
|---|---|---|---|---|
| 16 | 6% | 36% | 93% | - |
| 40 | 2% | 15% | 54% | 98% |

**Frequency tags** (0.133 MHz per site): the relative phase after a pulse is 0.08 / 0.84 / 8.4 rad for adjacent
sites at 0.1 / 1 / 10 us, and 3.4 / 34 / 335 rad for sites 40 apart.

**Channels and power.**
- PIC: one cascaded switch (two modulators) per site plus slow router trims; 16 sites at 20 uW need about 1 mW of
  input (module 4).
- AOD: 2 RF channels per pair.
- DMD: one device with digital control.

The PIC needs as many fast channels as sites; at large site counts that is the cost of this lead.

## 4. Limits

- **Patterns.** Random patterns with independent sites are taken as the workload.
- **AOD intermodulation.** The level is not sourced; we report where the products land, not how bright they are.
- **PIC settling.** The settling time (0.19 us) rests on an assumed mechanical mode (20 MHz, Q 10) and ZVD shaping
  (module 5).
- **SLM/DMD figures.** These come from one demonstration, read in abstract only.
- **Wavelength.** The measured PIC figures are at 795 nm, not 894.6 nm.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| rounds per pattern, times | Sec. 1, 2 | SIM | results/m31_patterns.json |
| third-order product coverage | Sec. 1, 3 | SIM | results/m31_leakage.json |
| frequency tags and phases | Sec. 1, 3 | DERIVED | results/m31_frequency.json |
| 71.4 dB, -68 dB, -50.8 dB, 26 ns, 8 channels (795 nm) | Sec. 1 | SOURCE | Zhao et al., arXiv:2508.09920, abstract read 2026-11-15 |
| 46 dB, -44 dB, 43 kHz, > 10,000 spots | Sec. 1 | SOURCE | Zhang et al., arXiv:2310.08539 (Optica 2024), abstract read 2026-11-15 |
| AOD access 4.5 / 7.5 us; 20 MHz band, 150 spots | Sec. 1, 3 | SOURCE | 45035-x-6.5DEG-1.06 (sintec.sg/ao/497.html, read 2026-11-14) |
| 0.19 us settling (ZVD, 5% mode error) | Sec. 1 | SIM | m5 (assumed f0 20 MHz, Q 10) |
| >= 88 dB cascade, 70 dB stray floor | Sec. 1 | SIM / CHOICE | m4 (docs/findings.md Sec. 5) |
| AODs limited to row or column patterns, intermodulation limits contrast | Sec. 1 | SEARCH SUMMARY | Radnaev et al., arXiv:2408.08288 (introduction not read) |

Open placeholders: [VALUE?] AOD third-order intermodulation level at 5 to 20 tones; [REF?] PIC switch performance
at 894.6 nm.
