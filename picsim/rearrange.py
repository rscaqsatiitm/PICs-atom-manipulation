"""Rearrangement scheduler for a stacked focal-plane OPA serving a static 2D array.

Hardware model (docs/scaling_beyond.md, docs/transport3d.md):
  * the static array has R rows x C columns (4 um pitch);
  * a stack of S chips addresses S consecutive rows at a time (row select places the stack); each stack row
    carries any number of independent tweezers along x (one per site at most);
  * all tweezers of the stack share the y motion: pick-up, a 2 um step into the lane between rows, optional
    travel of j rows along mid-column lanes (all picked atoms shift by the same j), a 2 um step back, drop-off;
  * within a stack row the x moves preserve order (tweezers never cross), which also keeps them >= 4 um apart.

Algorithm (greedy, two stages):
  1. vertical balancing: passes (j, k) that move surplus atoms of rows r in [k, k+S) to rows r + j that lack atoms
     for the centred target block (R_t x C_t), choosing (j, k) that moves the most atoms; dropped atoms go
     straight to empty target columns where possible;
  2. horizontal compaction: one pass per group of S target rows; in each row the atoms are matched in order to
     the target columns, extras pushed outside the block.
Every move (pick-up to drop-off) of an atom is counted, with its time, so that heating and loss per move can be
accumulated.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Timing:
    t_ramp: float = 100e-6
    t_lane: float = 400e-6
    t_x40: float = 400e-6  # time for a 40 um x move (scaled as sqrt(d/40 um), minimum t_min)
    t_min: float = 200e-6
    t_reposition: float = 300e-6  # row-select scanner settling between passes (CHOICE)
    pitch: float = 4e-6

    def move(self, d):
        if d <= 0:
            return 0.0
        return max(self.t_min, self.t_x40 * np.sqrt(d / 40e-6))

    def pass_time(self, dx_max, dy):
        # two ramps, two lane steps, x to the mid-column, y travel, x to the destination column
        return 2 * self.t_ramp + 2 * self.t_lane + self.move(dx_max) + (self.move(dy) + self.move(self.pitch / 2)
                                                                        if dy else 0.0) + self.t_reposition


@dataclass
class Result:
    passes_vertical: int = 0
    passes_horizontal: int = 0
    time_s: float = 0.0
    moves_per_atom: list = field(default_factory=list)
    target_filled: bool = False
    deficit: int = 0


def _monotone_match(src, dst):
    """Match sorted sources to a sorted subset of destinations of equal size, preserving order, minimising the
    largest displacement (simple: choose the window of destinations closest to the sources' span)."""
    src, dst = np.sort(src), np.sort(dst)
    n = len(src)
    if n == 0 or len(dst) < n:
        return None
    best, best_cost = None, np.inf
    for a in range(len(dst) - n + 1):
        cost = np.max(np.abs(dst[a:a + n] - src))
        if cost < best_cost:
            best, best_cost = dst[a:a + n], cost
    return list(zip(src, best))


def _transfers(n, in_t, C_t, S):
    """Candidate vertical passes: T[jj, k] = atoms that a pass with row shift j = J[jj] and stack position k moves
    from surplus rows into target rows lacking atoms (counts only)."""
    R = len(n)
    sur = np.clip(np.where(in_t, n - C_t, n), 0, None)
    dfc = np.where(in_t, np.clip(C_t - n, 0, None), 0)
    J = np.array([j for j in range(-(R - 1), R) if j != 0])
    T = np.zeros((len(J), R - S + 1), dtype=int)
    for a, j in enumerate(J):
        t = np.zeros(R, dtype=int)
        r = np.arange(R)
        ok = (r + j >= 0) & (r + j < R)
        t[ok] = np.minimum(sur[ok], dfc[r[ok] + j])
        cs = np.concatenate([[0], np.cumsum(t)])
        T[a] = cs[S:] - cs[:-S]
    return J, T


def _apply_counts(n, in_t, C_t, S, j, k):
    n = n.copy()
    sur = np.clip(np.where(in_t, n - C_t, n), 0, None)
    dfc = np.where(in_t, np.clip(C_t - n, 0, None), 0)
    for r in range(k, k + S):
        rr = r + j
        if 0 <= rr < len(n):
            m = min(sur[r], dfc[rr])
            n[r] -= m
            n[rr] += m
    return n


def _choose(n, in_t, C_t, S, timing, mode, n_top=8):
    J, T = _transfers(n, in_t, C_t, S)
    if T.max() <= 0:
        return None
    t_pass = np.array([timing.pass_time(40e-6, abs(j) * timing.pitch) for j in J])[:, None]
    if mode == "greedy":  # most atoms, ties to the smallest |j|
        score = T - 1e-3 * np.abs(J)[:, None]
    elif mode == "rate":  # most atoms per unit pass time
        score = T / t_pass
    elif mode == "lookahead":  # most atoms over this pass and the best next one
        flat = np.argsort(T, axis=None)[::-1][:n_top]
        best, bsc = None, -np.inf
        for f in flat:
            a, k = np.unravel_index(f, T.shape)
            if T[a, k] <= 0:
                continue
            n2 = _apply_counts(n, in_t, C_t, S, J[a], k)
            _, T2 = _transfers(n2, in_t, C_t, S)
            sc = T[a, k] + T2.max() - 1e-3 * abs(J[a])
            if sc > bsc:
                best, bsc = (int(T[a, k]), int(J[a]), int(k)), sc
        return best
    else:
        raise ValueError(mode)
    a, k = np.unravel_index(np.argmax(score), score.shape)
    return int(T[a, k]), int(J[a]), int(k)


def schedule(occ, R_t, C_t, S, timing=Timing(), max_passes=200, mode="greedy"):
    """occ: boolean array (R, C) of loaded sites; returns Result and the final occupancy and move counts.
    mode: pass selection in the balancing stage ('greedy', 'rate' or 'lookahead')."""
    R, C = occ.shape
    occ = occ.copy()
    moves = np.zeros((R, C), dtype=int)  # move count of the atom currently at each site
    r0, c0 = (R - R_t) // 2, (C - C_t) // 2
    trows = np.arange(r0, r0 + R_t)
    tcols = np.arange(c0, c0 + C_t)
    in_t = np.zeros(R, bool)
    in_t[trows] = True
    res = Result()

    def surplus(r):
        n = occ[r].sum()
        return n - C_t if in_t[r] else n

    def deficit(r):
        return max(0, C_t - occ[r].sum()) if in_t[r] else 0

    # ---------------- stage 1: vertical balancing
    for _ in range(max_passes):
        best = _choose(occ.sum(axis=1), in_t, C_t, S, timing, mode)
        if best is None:
            break
        _, j, k = best
        dxm = 0.0
        arrived = np.zeros((R, C), bool)
        for r in range(k, k + S):
            rr = r + j
            if not (0 <= rr < R and in_t[rr]):
                continue
            m = min(max(surplus(r), 0), deficit(rr))
            if m <= 0:
                continue
            cols = np.nonzero(occ[r] & ~arrived[r])[0]
            if in_t[r]:  # a target row gives away atoms outside the target columns first
                cols = sorted(cols, key=lambda c: (c0 <= c < c0 + C_t, abs(c - (c0 + C_t / 2))))
            else:
                cols = sorted(cols, key=lambda c: abs(c - (c0 + C_t / 2)))
            src = np.array(sorted(cols[:m]))
            empt = [c for c in tcols if not occ[rr, c]]
            if len(empt) < m:
                empt = empt + [c for c in range(C) if not occ[rr, c] and c not in empt]
            pairs = _monotone_match(src, np.array(empt))
            if pairs is None:
                continue
            for a, b in pairs:
                occ[r, a] = False
                occ[rr, b] = True
                arrived[rr, b] = True
                moves[rr, b] = moves[r, a] + 1
                moves[r, a] = 0
                dxm = max(dxm, abs(b - a) * timing.pitch)
        res.passes_vertical += 1
        res.time_s += timing.pass_time(dxm, abs(j) * timing.pitch)

    # ---------------- stage 2: horizontal compaction of the target rows
    for g in range(0, R_t, S):
        dxm = 0.0
        any_move = False
        for r in trows[g:g + S]:
            cols = np.nonzero(occ[r])[0]
            n = len(cols)
            if n < C_t:
                continue
            extra = n - C_t
            best = None
            for L in range(0, extra + 1):  # L extras to the left of the block, extra - L to the right
                if L > c0 or extra - L > C - (c0 + C_t):
                    continue
                dst = list(range(c0 - L, c0)) + list(tcols) + list(range(c0 + C_t, c0 + C_t + extra - L))
                cost = np.max(np.abs(np.array(dst) - cols))
                if best is None or cost < best[0]:
                    best = (cost, dst)
            if best is None:
                continue
            new = np.zeros(C, bool)
            newm = np.zeros(C, int)
            for a, b in zip(cols, best[1]):
                new[b] = True
                newm[b] = moves[r, a] + (1 if a != b else 0)
                if a != b:
                    any_move = True
                    dxm = max(dxm, abs(b - a) * timing.pitch)
            occ[r], moves[r] = new, newm
        if any_move:
            res.passes_horizontal += 1
            res.time_s += timing.pass_time(dxm, 0.0)

    block = occ[np.ix_(trows, tcols)]
    res.target_filled = bool(block.all())
    res.deficit = int((~block).sum())
    res.moves_per_atom = moves[np.ix_(trows, tcols)][block].tolist()
    return res, occ, moves


# ------------------------------------------------------------------ losses, imaging and repair
def lose(occ, rng, p_keep):
    """Remove atoms independently with keep probability p_keep (array or scalar)."""
    return occ & (rng.random(occ.shape) < p_keep)


def repair(occ, R_t, C_t, S, timing=Timing(), max_passes=50):
    """Fill the holes of the block directly from reservoir atoms (outside the block), one move per repaired atom.

    Each pass has a common row shift j (j = 0 allowed: a reservoir atom in the same row) and stack position k; in
    every stack row the chosen reservoir atoms are matched in order to holes of row r + j (the lanes let them pass
    other atoms). Returns the new occupancy, a mask of repaired sites, passes and time."""
    R, C = occ.shape
    occ = occ.copy()
    r0, c0 = (R - R_t) // 2, (C - C_t) // 2
    blk = np.zeros((R, C), bool)
    blk[r0:r0 + R_t, c0:c0 + C_t] = True
    repaired = np.zeros((R, C), bool)
    pairs = []  # ((r, a), (rr, b)) source and destination of every repair move
    passes, t = 0, 0.0
    for _ in range(max_passes):
        res_n = (occ & ~blk).sum(axis=1)
        hole_n = (~occ & blk).sum(axis=1)
        if hole_n.sum() == 0 or res_n.sum() == 0:
            break
        best = None
        for j in range(-(R - 1), R):
            tr = np.zeros(R, dtype=int)
            r = np.arange(R)
            ok = (r + j >= 0) & (r + j < R)
            tr[ok] = np.minimum(res_n[ok], hole_n[r[ok] + j])
            cs = np.concatenate([[0], np.cumsum(tr)])
            win = cs[S:] - cs[:-S]
            k = int(np.argmax(win))
            sc = win[k] - 1e-3 * abs(j)
            if win[k] > 0 and (best is None or sc > best[0]):
                best = (sc, j, k)
        if best is None:
            break
        _, j, k = best
        dxm = 0.0
        for r in range(k, k + S):
            rr = r + j
            if not (0 <= rr < R):
                continue
            holes = np.nonzero(~occ[rr] & blk[rr])[0]
            resv = np.nonzero(occ[r] & ~blk[r])[0]
            m = min(len(holes), len(resv))
            if m == 0:
                continue
            centre = holes.mean()
            resv = np.sort(resv[np.argsort(np.abs(resv - centre))[:m]])
            for a, b in zip(resv, np.sort(holes[:m])):
                occ[r, a] = False
                occ[rr, b] = True
                repaired[rr, b] = True
                pairs.append(((r, a), (rr, b)))
                dxm = max(dxm, abs(b - a) * timing.pitch)
        passes += 1
        t += timing.pass_time(dxm, abs(j) * timing.pitch)
    return occ, repaired, passes, t, pairs


def cycle_with_repair(occ0, R_t, C_t, S, rng, p_move, tau, eps_img, t_img=10e-3, rounds=2, timing=Timing(),
                      sched=None):
    """One rearrangement cycle followed by `rounds` of (image, repair). Detection precedes the imaging-induced loss,
    so atoms lost by an image are not seen. Returns defects after the cycle and after each round, and times."""
    R, C = occ0.shape
    r0, c0 = (R - R_t) // 2, (C - C_t) // 2
    blk = np.zeros((R, C), bool)
    blk[r0:r0 + R_t, c0:c0 + C_t] = True
    res, occ, moves = sched if sched is not None else schedule(occ0, R_t, C_t, S, timing)
    occ = lose(occ, rng, (1 - p_move) ** moves * np.exp(-res.time_s / tau))
    out = {"t_cycle_ms": res.time_s * 1e3, "defects": [int((~occ & blk).sum())], "repair_passes": [], "t_round_ms": []}
    for _ in range(rounds):
        seen = occ.copy()  # image: detection
        occ = lose(occ, rng, (1 - eps_img) * np.exp(-t_img / tau))  # imaging-induced and vacuum loss, unseen
        _, _, npass, trep, pairs = repair(seen, R_t, C_t, S, timing)
        moved = np.zeros_like(occ)
        for (r, a), (rr, b) in pairs:  # a repair succeeds only if its source atom survived the image
            if occ[r, a]:
                occ[r, a] = False
                occ[rr, b] = True
                moved[rr, b] = True
        occ = lose(occ, rng, np.where(moved, 1 - p_move, 1.0) * np.exp(-trep / tau))
        out["defects"].append(int((~occ & blk).sum()))
        out["repair_passes"].append(npass)
        out["t_round_ms"].append((t_img + trep) * 1e3)
    return out
