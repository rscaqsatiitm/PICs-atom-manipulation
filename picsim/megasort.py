"""Steady-state rearrangement of large arrays (1e4 to 1e6 sites) by a block of transport engines.

Layout: the array has R x C sites; rows are grouped in bands of B rows, the first B_c of each band computation
rows and the rest reservoir rows, so every computation row has a reservoir row within a few pitches.
Engines: a block of S rows x Mx engines per row; each engine serves one row segment of `seg` sites (its reach);
with Mx > 1 the engines of a row address every other segment (shared slit plane, m22 tiling), so a block position
covers segments g0, g0 + 2, ... The block is placed by the field-select scanner.

Each visit runs repair passes inside the block: a pass has one common row shift j (all tweezers share the y motion),
and in every (row, segment) the chosen reservoir atoms are matched in order to the known holes of row r + j in the
same segment; one move per refilled hole. Only holes and reservoir atoms seen in the last image are acted on.

Block modes (move constraints; everything else, timing, losses, imaging, is shared):
  pic       the stacked PIC engines above (default)
  aod_row   crossed AODs over a field of S rows x seg columns (Mx = 1): per pass one row tone and many column tones,
            i.e. the reservoir atoms of one source row move by one common row shift to the holes of one destination
            row, matched in order (column tones do not cross). Multi-row passes with a shared column pattern are not
            used, so this is a lower bound on what crossed AODs can do.
  single    one moving tweezer in the field: one atom per pass, the (reservoir, hole) pair of shortest pass time
  parallel  idealised ceiling: every known hole in the field is matched in one pass to a nearby unused reservoir
            atom, with independent 2D moves (crossing allowed); not a hardware model

Losses: vacuum (lifetime tau, every atom, Poisson-sampled per interval), imaging (eps per atom per image, after
detection), per move (p_move). Reservoir rows are refreshed every t_reload to loading fraction f (an external
supply; its feasibility is a separate question). Freshly loaded reservoir atoms are taken as known at once (the
reservoir is assumed imaged as part of the reload). Time per pass from picsim.rearrange.Timing.
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
    mode: str = "pic"  # "pic", "aod_row", "single", "parallel" (see the module docstring)

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
    p_move_fast: float | None = None  # loss per move with `timing_fast`; None: p_move
    max_passes_single: int = 400  # per visit in mode "single" (also stopped by the next image)
    skip_reposition_same: bool = False  # no scanner reposition when the block visits the same position again


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
        self._last_pos = None
        self.pass_tweezers = []  # number of atoms moved in each pass (simultaneous moving tweezers)

    # ------------------------------------------------------------------ losses
    def _vacuum(self, dt):
        if dt <= 0:
            return
        idx = np.flatnonzero(self.occ)
        k = self.rng.poisson(len(idx) * (1 - np.exp(-dt / self.P.tau)))
        if k:
            self.occ.flat[self.rng.choice(idx, size=min(k, len(idx)), replace=False)] = False

    def _advance(self, dt):
        d0 = self.defect_fraction()
        self._vacuum(dt)  # losses are uniform in time over dt, so the trapezoid rule is exact in expectation
        self._integ += 0.5 * dt * (d0 + self.defect_fraction())
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
            p = self.P.p_move if self.P.p_move_fast is None else self.P.p_move_fast
            return self.P.timing_fast, self.P.heat_fast, p
        return self.P.timing, self.P.heat_slow, self.P.p_move

    def visit(self, k, gs):
        mode = self.Bk.mode
        if mode != "pic":
            return {"aod_row": self._visit_aod_row, "single": self._visit_single,
                    "parallel": self._visit_parallel}[mode](k, gs)
        S, seg = self.Bk.S, self.Bk.seg
        tm, h_move, p_move = self._timing()
        py = tm.pitch if tm.pitch_y is None else tm.pitch_y
        rows = np.arange(k, min(k + S, self.L.R))
        cols = np.concatenate([np.arange(g * seg, min((g + 1) * seg, self.L.C)) for g in gs])
        t_used = self._reposition((k, gs), tm)
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
            n_tw = 0
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
                    n_tw += m
                    for a_, b_ in zip(src, dst):
                        self.known_res[r_src, a_] = False
                        self.known_holes[r_dst, b_] = False
                        self.n_moves += 1
                        if self.occ[r_src, a_]:
                            self.occ[r_src, a_] = False
                            self.consumed += 1
                            if self.rng.random() > p_move:
                                self.occ[r_dst, b_] = True
                                self.moves[r_dst, b_] = self.moves[r_src, a_] + 1
                                self.heat[r_dst, b_] = self.heat[r_src, a_] + h_move
                                self.n_fills += 1
                        dxm = max(dxm, abs(int(b_) - int(a_)) * tm.pitch)
            t_used += tm.pass_time(dxm, abs(j) * py) - tm.t_reposition
            n_pass += 1
            self.pass_tweezers.append(n_tw)
        return t_used, n_pass

    # ------------------------------------------------------------------ other move constraints
    def _reposition(self, pos, tm):
        same = self.P.skip_reposition_same and pos == self._last_pos
        self._last_pos = pos
        return 0.0 if same else tm.t_reposition

    def _field(self, k, gs):
        seg = self.Bk.seg
        rows = np.arange(k, min(k + self.Bk.S, self.L.R))
        cols = np.concatenate([np.arange(g * seg, min((g + 1) * seg, self.L.C)) for g in gs])
        return rows, cols

    def _move(self, r_src, a_, r_dst, b_, h_move, p_move):
        """One atom from (r_src, a_) to (r_dst, b_): the bookkeeping of the pic visit."""
        self.known_res[r_src, a_] = False
        self.known_holes[r_dst, b_] = False
        self.n_moves += 1
        if self.occ[r_src, a_]:
            self.occ[r_src, a_] = False
            self.consumed += 1
            if self.rng.random() > p_move:
                self.occ[r_dst, b_] = True
                self.moves[r_dst, b_] = self.moves[r_src, a_] + 1
                self.heat[r_dst, b_] = self.heat[r_src, a_] + h_move
                self.n_fills += 1

    @staticmethod
    def _tmove(tm, d):
        d = np.asarray(d, float)
        return np.where(d > 0, np.maximum(tm.t_min, tm.t_x40 * np.sqrt(np.maximum(d, 0) / 40e-6)), 0.0)

    def _pass_times(self, tm, dx, dy):
        """Vectorised Timing.pass_time for arrays of x and y travel (m)."""
        dy = np.asarray(dy, float)
        extra = np.where(dy > 0, self._tmove(tm, dy) + self._tmove(tm, tm.pitch / 2), 0.0)
        return 2 * tm.t_ramp + 2 * tm.t_lane + self._tmove(tm, dx) + extra + tm.t_reposition

    def _visit_aod_row(self, k, gs):
        tm, h_move, p_move = self._timing()
        py = tm.pitch if tm.pitch_y is None else tm.pitch_y
        rows, cols = self._field(k, gs)
        t_used = self._reposition((k, gs), tm)
        n_pass = 0
        for _ in range(self.P.max_passes):
            H = self.known_holes[np.ix_(rows, cols)]
            Rv = self.known_res[np.ix_(rows, cols)]
            hc, rc = H.sum(axis=1), Rv.sum(axis=1)
            if hc.sum() == 0 or rc.sum() == 0:
                break
            M = np.minimum(rc[:, None], hc[None, :]).astype(float)
            dist = np.abs(np.arange(len(rows))[:, None] - np.arange(len(rows))[None, :])
            score = M - 1e-3 * dist / max(len(rows), 1)  # most refills; ties: shortest row shift
            a, b = np.unravel_index(np.argmax(score), score.shape)
            m = int(M[a, b])
            if m == 0:
                break
            r_src, r_dst = rows[a], rows[b]
            src = cols[self.known_res[r_src, cols]]
            dst = cols[self.known_holes[r_dst, cols]]
            dst = np.sort(dst[:m]) if m == len(dst) else np.sort(dst[np.argsort(np.abs(dst - dst.mean()))[:m]])
            src = np.sort(src[np.argsort(np.abs(src - dst.mean()))[:m]])
            for a_, b_ in zip(src, dst):
                self._move(r_src, a_, r_dst, b_, h_move, p_move)
            dxm = float(np.max(np.abs(dst - src))) * tm.pitch
            t_used += tm.pass_time(dxm, abs(int(r_dst) - int(r_src)) * py) - tm.t_reposition
            n_pass += 1
            self.pass_tweezers.append(m)
        return t_used, n_pass

    def _visit_single(self, k, gs):
        tm, h_move, p_move = self._timing()
        py = tm.pitch if tm.pitch_y is None else tm.pitch_y
        rows, cols = self._field(k, gs)
        t_used = self._reposition((k, gs), tm)
        n_pass = 0
        for _ in range(self.P.max_passes_single):
            if n_pass and self.t + t_used >= self._next_img:
                break
            sub = np.ix_(rows, cols)
            hr, hcol = np.nonzero(self.known_holes[sub])
            rr, rcol = np.nonzero(self.known_res[sub])
            if len(hr) == 0 or len(rr) == 0:
                break
            dx = np.abs(cols[hcol][:, None] - cols[rcol][None, :]) * tm.pitch
            dy = np.abs(rows[hr][:, None] - rows[rr][None, :]) * py
            T = self._pass_times(tm, dx, dy)
            i, j = np.unravel_index(np.argmin(T), T.shape)
            self._move(rows[rr[j]], cols[rcol[j]], rows[hr[i]], cols[hcol[i]], h_move, p_move)
            t_used += float(T[i, j]) - tm.t_reposition
            n_pass += 1
            self.pass_tweezers.append(1)
        return t_used, n_pass

    def _visit_parallel(self, k, gs, k_near=8):
        from scipy.spatial import cKDTree

        tm, h_move, p_move = self._timing()
        py = tm.pitch if tm.pitch_y is None else tm.pitch_y
        rows, cols = self._field(k, gs)
        t_used = self._reposition((k, gs), tm)
        n_pass = 0
        for _ in range(self.P.max_passes):
            sub = np.ix_(rows, cols)
            hr, hcol = np.nonzero(self.known_holes[sub])
            rr, rcol = np.nonzero(self.known_res[sub])
            if len(hr) == 0 or len(rr) == 0:
                break
            hxy = np.c_[cols[hcol] * tm.pitch, rows[hr] * py]
            rxy = np.c_[cols[rcol] * tm.pitch, rows[rr] * py]
            kk = min(k_near, len(rr))
            dist, idx = cKDTree(rxy).query(hxy, k=kk)
            dist, idx = dist.reshape(len(hr), kk), idx.reshape(len(hr), kk)
            used = np.zeros(len(rr), bool)
            pairs = []
            for i in np.argsort(dist[:, 0], kind="stable"):  # closest holes first
                for q in range(kk):
                    j = idx[i, q]
                    if not used[j]:
                        used[j] = True
                        pairs.append((i, j))
                        break
            if not pairs:
                break
            dxm = dym = 0.0
            for i, j in pairs:
                self._move(rows[rr[j]], cols[rcol[j]], rows[hr[i]], cols[hcol[i]], h_move, p_move)
                dxm = max(dxm, abs(hxy[i, 0] - rxy[j, 0]))
                dym = max(dym, abs(hxy[i, 1] - rxy[j, 1]))
            t_used += tm.pass_time(dxm, dym) - tm.t_reposition
            n_pass += 1
            self.pass_tweezers.append(len(pairs))
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
            start = self._cursor  # one sweep in order from where the last one stopped (no position skipped)
            for step in range(n_pos):
                if self.t >= self._next_img:
                    break
                idx = (start + step) % n_pos
                k, gs = self.positions[idx]
                rows = slice(k, k + self.Bk.S)
                cols = np.concatenate([np.arange(g * self.Bk.seg, min((g + 1) * self.Bk.seg, self.L.C)) for g in gs])
                if not self.known_holes[rows][:, cols].any() or not self.known_res[rows][:, cols].any():
                    continue
                dt, npass = self.visit(k, gs)
                if npass:
                    work = True
                    self._advance(dt)
                    self.n_visits = getattr(self, "n_visits", 0) + 1
                    self._cursor = (idx + 1) % n_pos
            if not work:  # idle until the next image
                self._advance(max(self._next_img - self.t, 1e-4))
        return t_reach
