from stable.kromatic import (
    kromatic,
    oriented_kromatic,
    oriented_expander
)
from stable.utils import mn_G_expansion, schur_expansion
import itertools


def natural_unit_iterval_order_incompability_graph(n, r):
    v = list(range(1, n + 1))
    e = [(i, j) for i in v for j in v if 0 < j - i < r]
    return v, e


def graphs(n):
    v = list(range(1, n + 1))
    e = list(itertools.combinations(v, 2))
    for k in range(0, len(e) + 1):
        for edges in itertools.combinations(e, k):
            yield v, edges


def is_claw_free(vertices, edges):
    e = {v: set() for v in vertices}
    for a, b in edges:
        e[a].add(b)
        e[b].add(a)
    for v in vertices:
        for a in e[v]:
            for b in e[v]:
                if a == b:
                    continue
                for c in e[v]:
                    if a == c or b == c:
                        continue
                    if a not in e[b] and a not in e[c] and b not in e[c]:
                        return False
    return True


def test_G_expansion(nvars=3, nverts=3):
    for v, e in graphs(nverts):
        f = kromatic(nvars, v, e)
        exp = mn_G_expansion(f)
        icf = is_claw_free(v, e)
        pos = exp.is_nonnegative()
        print(v, e, icf, pos, exp)
        print()
        assert (not icf) or pos


def test_G_oriented_expansion(nvars=3, n=5, r=5):
    seen = set()
    for nn in range(1, n + 1):
        for rr in range(1, r + 1):
            v, e = natural_unit_iterval_order_incompability_graph(nn, rr)
            
            v, e = tuple(sorted(v)), tuple(sorted(e))
            if (v, e) not in seen:
                seen.add((v, e))
            else:
                continue
            
            icf = is_claw_free(v, e)
            try:
                f = oriented_kromatic(nvars, v, e)
                exp = oriented_expander(mn_G_expansion, f)
                pos = exp.is_nonnegative()
                print(v, e, icf, pos, exp)
                print()
                assert (not icf) or pos
            except:
                print(v, e, icf)
                print('pass\n')
                input('')
