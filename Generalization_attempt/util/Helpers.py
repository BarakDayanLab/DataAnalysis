import numpy as np


def _merge(a, b):
    if len(a) < len(b):
        b, a = a, b
    c = np.empty(len(a) + len(b), dtype=a.dtype)
    b_indices = np.arange(len(b)) + np.searchsorted(a, b)
    a_indices = np.ones(len(c), dtype=bool)
    a_indices[b_indices] = False
    c[b_indices] = b
    c[a_indices] = a
    return c


def merge(*x):
    if len(x) == 1:
        return np.array(x[0])
    return _merge(np.array(x[0]), merge(*x[1:]))


