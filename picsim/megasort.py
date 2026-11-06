"""Steady-state rearrangement of large arrays (1e4 to 1e6 sites) by a block of transport engines.

Layout: the array has R x C sites; rows are grouped in bands of B rows, the first B_c of each band computation
rows and the rest reservoir rows, so every computation row has a reservoir row within a few pitches.
Engines: a block of S rows x Mx engines per row; each engine serves one row segment of `seg` sites (its reach);
with Mx > 1 the engines of a row address every other segment (shared slit plane, m22 tiling), so a block position
covers segments g0, g0 + 2, ... The block is placed by the field-select scanner.

Each visit runs repair passes inside the block: a pass has one common row shift j (all tweezers share the y motion),
and in every (row, segment) the chosen reservoir atoms are matched in order to the known holes of row r + j in the
same segment; one move per refilled hole. Only holes and reservoir atoms seen in the last image are acted on.

Losses: vacuum (lifetime tau, every atom, Poisson-sampled per interval), imaging (eps per atom per image, after
detection), per move (p_move). Reservoir rows are refreshed every t_reload to loading fraction f (an external
supply; its feasibility is a separate question). Time per pass from picsim.rearrange.Timing.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .rearrange import Timing


@dataclass
class Layout:
    R: int
    C: int
    B: int  # band period (rows)
    B_c: int  # computation rows per band

    def masks(self):
        rows = np.arange(self.R) % self.B < self.B_c
        comp = np.repeat(rows[:, None], self.C, axis=1)
        return comp, ~comp


@dataclass
class Block:
    S: int  # rows
    Mx: int  # engines per row
    seg: int  # sites per segment

    def positions(self, R, C):
        nseg = int(np.ceil(C / self.seg))
        pos = []
        ks = list(range(0, max(R - self.S, 0) + 1, self.S))
        if ks[-1] + self.S < R:  # last block flush with the array edge, so it always spans S rows (and a reservoir)
            ks.append(R - self.S)
        for k in ks:
            if self.Mx == 1:
                for g in range(nseg):
                    pos.append((k, (g,)))
            else:
                for base in range(0, nseg, 2 * self.Mx):
                    for par in (0, 1):
                        gs = tuple(g for g in range(base + par, min(base + 2 * self.Mx, nseg), 2))
                        if gs:
                            pos.append((k, gs))
        return pos


@dataclass
class Params:
    tau: float = 300.0  # s
    eps_img: float = 1e-4
    p_move: float = 1e-4
    t_img: float = 10e-3  # s
    T_img: float = 100e-3  # s, imaging period
    t_reload: float = 1.0  # s, reservoir refresh period
    f_load: float = 0.5
    max_passes: int = 4  # per visit
    timing: Timing = field(default_factory=Timing)
    timing_fast: Timing | None = None  # if set, used while the last imaged defect fraction exceeds switch_defect
    switch_defect: float = 1e-2
    heat_slow: float = 0.3  # uK per move with `timing` (m19 quantum, m18 slow-lane 3D)
    heat_fast: float = 10.0  # uK per move with `timing_fast` (m18 classical 3D, 100 us steps, 0.3 mK)


class Sim:
    def __init__(self, layout: Layout, block: Block, prm: Params, seed=0):
        self.L, self.Bk, self.P = layout, block, prm
        self.rng = np.random.default_rng(seed)
        self.comp, self.res = layout.masks()
        self.occ = self.rng.random((layout.R, layout.C)) < prm.f_load
        self.known_holes = np.zeros_like(self.occ)
        self.known_res = np.zeros_like(self.occ)
        self.moves = np.zeros(self.occ.shape, dtype=np.int16)
        self.heat = np.zeros(self.occ.shape, dtype=np.float32)  # uK added by moves to the atom at each site
        self.t = 0.0
        self.n_moves = 0
        self.n_fills = 0
        self.consumed = 0
        self.positions = block.positions(layout.R, layout.C)
        self.seg_of_col = np.arange(layout.C) // block.seg
        self.log = []  # (t, defect fraction just after an image, just before it, time integral of the defect fraction)
        self._integ = 0.0
        self._next_img = 0.0
        self._next_reload = prm.t_reload
        self._cursor = 0  # round-robin position of the block (continues across images)

    # ------------------------------------------------------------------ losses
    def _vacuum(self, dt):
        if dt <= 0:
            return
        idx = np.flatnonzero(self.occ)
        k = self.rng.poisson(len(idx) * (1 - np.exp(-dt / self.P.tau)))
        if k:
            self.occ.flat[self.rng.choice(idx, size=min(k, len(idx)), replace=False)] = False

    def _advance(self, dt):
        self._integ += dt * self.defect_fraction()  # left Riemann sum; dt is one visit, image or idle interval
        self._vacuum(dt)
        self.t += dt

    def image(self):
        pre = self.defect_fraction()
        self.known_holes = self.comp & ~self.occ
        self.known_res = self.res & self.occ
        idx = np.flatnonzero(self.occ)
        k = self.rng.binomial(len(idx), self.P.eps_img)
        if k:
            self.occ.flat[self.rng.choice(idx, size=k, replace=False)] = False
        self._advance(self.P.t_img)
        self.log.append((self.t, self.defect_fraction(), pre, self._integ))

    def reload(self):
        empty = self.res & ~self.occ
        new = empty & (self.rng.random(self.occ.shape) < self.P.f_load)
        self.occ |= new
        self.known_res |= new
        self.heat[new] = 0.0  # fresh atoms
        self.moves[new] = 0

    def heat_stats(self):
        h = self.heat[self.comp & self.occ]
        return {"mean_uK": float(h.mean()), "p99_uK": float(np.percentile(h, 99)), "max_uK": float(h.max()),
                "moved_fraction": float((self.moves[self.comp & self.occ] > 0).mean())}

    def mean_defect(self, t0):
        """Time-averaged defect fraction from the first image at or after t0 to the last image."""
        d = np.array(self.log)
        i = np.searchsorted(d[:, 0], t0)
        return float((d[-1, 3] - d[i, 3]) / max(d[-1, 0] - d[i, 0], 1e-12))

    def defect_fraction(self):
        return float((self.comp & ~self.occ).sum() / self.comp.sum())

    # ------------------------------------------------------------------ one visit
    def _timing(self):
        if self.P.timing_fast is not None and (not self.log or self.log[-1][1] > self.P.switch_defect):
            return self.P.timing_fast, self.P.heat_fast
        return self.P.timing, self.P.heat_slow

    def visit(self, k, gs):
        S, seg = self.Bk.S, self.Bk.seg
        tm, h_move = self._timing()
        py = tm.pitch if tm.pitch_y is None else tm.pitch_y
        rows = np.arange(k, min(k + S, self.L.R))
        cols = np.concatenate([np.arange(g * seg, min((g + 1) * seg, self.L.C)) for g in gs])
        t_used = tm.t_reposition
        n_pass = 0
        for _ in range(self.P.max_passes):
            H = self.known_holes[np.ix_(rows, cols)]
            Rv = self.known_res[np.ix_(rows, cols)]
            sg = self.seg_of_col[cols]
            useg = np.unique(sg)
            hc = np.stack([H[:, sg == g].sum(axis=1) for g in useg], axis=1)  # (rows, segs)
            rc = np.stack([Rv[:, sg == g].sum(axis=1) for g in useg], axis=1)
            if hc.sum() == 0 or rc.sum() == 0:
                break
            best, bj = 0, None
            n = len(rows)
            for j in range(-(n - 1), n):
                a = np.arange(max(0, -j), min(n, n - j))
                tot = np.minimum(rc[a], hc[a + j]).sum()
                if tot > best or (tot == best and bj is not None and abs(j) < abs(bj)):
                    best, bj = tot, j
            if best == 0:
                break
            dxm = 0.0
            j = bj
            for ia in range(max(0, -j), min(n, n - j)):
                r_src, r_dst = rows[ia], rows[ia + j]
                for g in useg:
                    cg = cols[sg == g]
                    src = cg[self.known_res[r_src, cg]]
                    dst = cg[self.known_holes[r_dst, cg]]
                    m = min(len(src), len(dst))
                    if m == 0:
                        continue
                    dst = np.sort(dst[:m]) if m == len(dst) else np.sort(dst[np.argsort(np.abs(dst - dst.mean()))[:m]])
                    src = np.sort(src[np.argsort(np.abs(src - dst.mean()))[:m]])
                    for a_, b_ in zip(src, dst):
                        self.known_res[r_src, a_] = False
                        self.known_holes[r_dst, b_] = False
                        self.n_moves += 1
                        if self.occ[r_src, a_]:
                            self.occ[r_src, a_] = False
                            self.consumed += 1
                            if self.rng.random() > self.P.p_move:
                                self.occ[r_dst, b_] = True
                                self.moves[r_dst, b_] = self.moves[r_src, a_] + 1
                                self.heat[r_dst, b_] = self.heat[r_src, a_] + h_move
                                self.n_fills += 1
                        dxm = max(dxm, abs(int(b_) - int(a_)) * tm.pitch)
            t_used += tm.pass_time(dxm, abs(j) * py) - tm.t_reposition
            n_pass += 1
        return t_used, n_pass

    # ------------------------------------------------------------------ run
    def run(self, t_end, target=None):
        """Run until t_end (s). Returns the time at which the defect fraction first fell below `target`."""
        t_reach = None
        while self.t < t_end:
            if self.t >= self._next_img:
                self.image()
                self._next_img = self.t + self.P.T_img
                if target is not None and t_reach is None and self.log[-1][1] <= target:
                    t_reach = self.t
            if self.t >= self._next_reload:
                self.reload()
                self._next_reload += self.P.t_reload
            work = False
            n_pos = len(self.positions)
            for step in range(n_pos):
                if self.t >= self._next_img:
                    break
                k, gs = self.positions[(self._cursor + step) % n_pos]
                rows = slice(k, k + self.Bk.S)
                cols = np.concatenate([np.arange(g * self.Bk.seg, min((g + 1) * self.Bk.seg, self.L.C)) for g in gs])
                if not self.known_holes[rows][:, cols].any() or not self.known_res[rows][:, cols].any():
                    continue
                dt, npass = self.visit(k, gs)
                if npass:
                    work = True
                    self._advance(dt)
                    self.n_visits = getattr(self, "n_visits", 0) + 1
                    self._cursor = (self._cursor + step + 1) % n_pos
            if not work:  # idle until the next image
                self._advance(max(self._next_img - self.t, 1e-4))
        return t_reach
