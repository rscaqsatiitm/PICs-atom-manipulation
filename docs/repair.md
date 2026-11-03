# Repair rounds, scheduler variants and the axial handoff

Follow-up to `docs/rearrangement.md`. Code: `picsim/rearrange.py` (`repair`, `cycle_with_repair`,
selection modes), `scripts/m21_repair.py`, `scripts/m19_quantum.py axial`; results `results/m21_*.json`,
`results/m19_axial.json`. Values are SIM and valid in the regime stated. They should be reviewed before they
are quoted.

## 1. Result

After the rearrangement cycle (32 x 32 block from a 50 x 45 array at 50% loading, 17 ms), one round of imaging and
direct repair takes 10 to 13 ms: a 10 ms image and 0.1 to 1.3 repair passes. Each repair moves a single reservoir
atom through the lanes into a hole. The repair works, but it cannot get below a floor set by the image itself.
Atoms lost because of the image are not seen by it, and every atom is also exposed to vacuum loss for the round's
duration. The floor is

    defects ~ N (eps_img + t_round / tau) + p_move x (repairs),   N = 1024,

about 0.14 defects for eps_img = 1e-4 and tau = 300 s. A second round does not lower it. Hence the probability of
a defect-free 32 x 32 block is set by the imaging loss and the vacuum lifetime, not by the transport, once the
per-move loss is about 1e-4 or less:

| per-move loss | vacuum lifetime | imaging loss | P(defect-free): after cycle / after 1 repair / after 2 | defects after 1 repair |
|---|---|---|---|---|
| 1e-3 | 300 s | 1e-4 | 0.27 / **0.87** / 0.87 | 0.14 |
| 1e-4 | 300 s | 1e-4 | **0.84** / 0.87 / 0.87 | 0.14 |
| 1e-3 | 300 s | 1e-3 | 0.27 / 0.35 / 0.33 | 1.04 |
| 1e-4 | 300 s | 1e-3 | **0.84 / 0.35** / 0.33 | 1.03 |
| 1e-3 | 30 s | 1e-4 | 0.14 / 0.59 / 0.62 | 0.53 |
| 1e-4 | 30 s | 1e-4 | 0.48 / 0.62 / 0.63 | 0.48 |

(40 loadings x 25 loss draws per row.) With good transport (1e-4) and an imaging loss of 1e-3, a repair round
makes things worse: 0.84 becomes 0.35, because the image costs about one atom of the block. Imaging should then be
skipped, or used only to post-select. For a 99% defect-free block one needs roughly 1e-5 for each term:
per-move loss, imaging loss per image, and t_round / tau. These are properties of the platform, not of the photonic
chip. The chip's contribution is 1.1 moves per atom and 17 ms per cycle.

## 2. Scheduler variants

Balancing-pass selection on 40 loadings (32 x 32, 50%, one stack):

| mode | passes (all) | cycle time, mean / max |
|---|---|---|
| greedy (most atoms per pass) | 8.05 | 17.8 / 23.1 ms |
| rate (most atoms per unit pass time) | 8.05 | 17.6 / 22.8 ms |
| lookahead (best two-pass sequence among the 8 best first passes) | 7.98 | 17.6 / 22.9 ms |

The three modes agree to within 1%. The balancing stage has a floor of about four passes, set by the geometry.
In one pass the 16-row window can hold at most 9 non-target rows above or below the 32-row target band, so it can
feed at most 9 target rows. Feeding all 32 therefore takes at least four passes, plus two compaction passes. The
greedy schedule (5.8 + 2) is within about 30% of that floor. Further gains would come from shorter passes, not
fewer: ramps, lane steps and scanner settling make up 1.3 of the 2.2 ms per pass.

## 3. Axial motion during the handoff

Pending: the 1D quantum run along the optical axis (`scripts/m19_quantum.py axial`) is in progress; its results will be added here and stored in `results/m19_axial.json`.

## 4. Limits and next steps

The repair model moves one reservoir atom per hole, through the lanes; repaired atoms come from the nearest rows
under the pass rules. Detection is assumed perfect apart from the imaging-induced loss, with no false positives or
negatives. The imaging time (10 ms) and losses are CHOICE values, not measured. The axial check is 1D, with the atom
assumed to sit at the moving-trap centre laterally.

## Verification notes

| value | location | type | pointer |
|---|---|---|---|
| repair table | Sec. 1 | SIM | results/m21_repair.json (review) |
| scheduler modes | Sec. 2 | SIM | results/m21_modes.json |
| axial quantum and Wigner results | Sec. 3 | SIM | results/m19_axial.json (review) |
| eps_img 1e-3 / 1e-4, t_img 10 ms, tau 30 / 300 s, p_move 1e-3 / 1e-4 | Sec. 1 | CHOICE | assumptions |
