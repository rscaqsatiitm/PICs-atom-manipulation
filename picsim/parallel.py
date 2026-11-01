"""Process-parallel map for independent mode solves (gmsh and scipy are not thread-safe)."""
from __future__ import annotations

import multiprocessing as mp
import os


def pmap(func, items, processes=None):
    processes = processes or min(4, os.cpu_count() or 1)
    if processes == 1 or len(items) == 1:
        return [func(i) for i in items]
    # fork all workers up front (no maxtasksperchild): re-forking a parent whose pool threads are running
    # can copy a held lock into the child and deadlock it
    ctx = mp.get_context("fork")
    with ctx.Pool(processes) as pool:
        return pool.map(func, items, chunksize=1)
