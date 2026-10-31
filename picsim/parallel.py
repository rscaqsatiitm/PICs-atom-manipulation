"""Process-parallel map for independent mode solves (gmsh and scipy are not thread-safe)."""
from __future__ import annotations

import multiprocessing as mp
import os


def pmap(func, items, processes=None):
    processes = processes or min(4, os.cpu_count() or 1)
    if processes == 1 or len(items) == 1:
        return [func(i) for i in items]
    ctx = mp.get_context("fork")
    with ctx.Pool(processes, maxtasksperchild=8) as pool:
        return pool.map(func, items, chunksize=1)
