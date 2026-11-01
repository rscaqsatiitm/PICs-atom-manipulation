"""Process-parallel map for independent mode solves (gmsh and scipy are not thread-safe)."""
from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

# one BLAS/OpenMP thread per worker process; set before numpy/scipy are imported by scripts
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")


def pmap(func, items, processes=None):
    """Ordered parallel map. A worker that dies (e.g. out of memory) raises BrokenProcessPool instead of
    hanging, which multiprocessing.Pool does when it re-forks replacement workers."""
    processes = processes or min(4, os.cpu_count() or 1)
    items = list(items)
    if processes == 1 or len(items) == 1:
        return [func(i) for i in items]
    with ProcessPoolExecutor(max_workers=processes, mp_context=mp.get_context("fork")) as ex:
        return list(ex.map(func, items, chunksize=1))
