from .vectors import Vector
from .partitions import Partition


def sh_horizontal_strip(a, b):
    assert 0 <= a <= b
    def op(mu):
        if mu is None: return None
        i = 1
        j = 1 + a
        while Partition.isin_shifted(mu, i, j):
            i += 1
            j += 1
        if (i == j or Partition.isin_shifted(mu, i, j - 1)) and (i == 1 or Partition.isin_shifted(mu, i - 1, j + b - a)):
            return mu[:i - 1] + (Partition.get(mu, i) + b - a + 1,) + mu[i:]
        return None
    return op


def sh_vertical_strip(a, b):
    assert 0 <= a <= b
    def op(mu):
        if mu is None: return None
        i = 1
        j = 1 + a
        while Partition.isin_shifted(mu, i, j):
            i += 1
            j += 1
        h = i + a - b
        if Partition.isin_shifted(mu, h, j):
            return None
        if (i == j or Partition.isin_shifted(mu, i, j - 1)) and (h == 1 or Partition.isin_shifted(mu, h - 1, j)):
            return mu[:h - 1] + tuple(Partition.get(mu, k) + 1 for k in range(h, i + 1)) + mu[i:]
        return None
    return op

