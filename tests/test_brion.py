from clans import Clan
from permutations import Permutation


def test(slow=False):
    _test_AI(1)
    _test_AI(2)
    _test_AI(3)
    _test_AI(4)
    _test_AI(5)
    _test_AI(6)

    _test_AII(1)
    _test_AII(3)
    _test_AII(5)
    _test_AII(7)

    if slow:
        _test_AI(7)
        _test_AI(8)

        _test_AII(9)
        _test_AII(11)


def precsim(k, n):
    def span(w):
        o = [w(i) for i in range(1, n + 1)]
        for i in range(k, len(o) - 2):
            b, c, a = o[i:i + 3]
            if a < b < c:
                new_o = o[:i] + [c, a, b] + o[i + 3:]
                yield w.__class__(*new_o)
    return span


def precapprox(k, n):
    def span(w):
        o = [w(i) for i in range(1, n + 1)]
        for i in range(k, len(o) - 3, 2):
            b, c, a, d = o[i:i + 4]
            if a < b < c < d:
                new_o = o[:i] + [a, d, b, c] + o[i + 4:]
                yield w.__class__(*new_o)
    return span


def transitive_closure(*args):
    def span(w):
        for spanner in args:
            for v in spanner(w):
                yield v
    return span


def _test_AI(rank=3):
    print('testing AI, rank =', rank)
    gamma_set = set(Permutation.involutions(rank + 1))
    invol_set = gamma_set
    rs_fn = lambda x: x
    brion_fn = lambda z: {w.inverse() for w in z.get_atoms()}
    extended_brion_fn = brion_fn
    matchings_fn = lambda z: {()}
    shape_fn = lambda w: ()
    is_aligned_fn = lambda gamma, m: True
    generator_fn = lambda z, m: z.get_max_atom().inverse()
    span_fn = precsim(0, rank + 1)
    _generic_test(gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, matchings_fn, shape_fn, is_aligned_fn, generator_fn, span_fn)


def _test_AII(rank=3):
    print('testing AII, rank =', rank)
    assert rank % 2 != 0
    gamma_set = set(Permutation.fpf_involutions(rank + 1))
    invol_set = gamma_set
    rs_fn = lambda x: x
    brion_fn = lambda z: {w.inverse() for w in z.get_fpf_atoms()}
    extended_brion_fn = brion_fn
    matchings_fn = lambda z: {()}
    shape_fn = lambda w: ()
    is_aligned_fn = lambda gamma, m: True
    generator_fn = lambda z, m: z.get_max_fpf_atom().inverse()
    span_fn = precapprox(0, rank + 1)
    _generic_test(gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, matchings_fn, shape_fn, is_aligned_fn, generator_fn, span_fn)


def span(generator, preorder_fn):
    ans = set()
    add = {generator}
    while add:
        newadd = set()
        for g in add:
            if g not in ans:
                ans.add(g)
                for h in preorder_fn(g):
                    if h not in ans or add:
                        newadd.add(h)
        add = newadd
    return ans


def _generic_test(
        gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, 
        matchings_fn, shape_fn, is_aligned_fn, 
        generator_fn, preorder_fn):
    expected_invol = {rs_fn(gamma) for gamma in gamma_set}
    assert expected_invol == invol_set

    matchings = {}
    for z in invol_set:
        extended_brion = extended_brion_fn(z)

        expected_matchings = set(matchings_fn(z))
        actual_matchings = {shape_fn(w) for w in extended_brion}
        assert expected_matchings == actual_matchings

        matchings[z] = {m: set() for m in actual_matchings}
        for w in extended_brion:
            m = shape_fn(w)
            matchings[z][m].add(w)

        for m in matchings[z]:
            gen = generator_fn(z, m)
            expected_block = span(gen, preorder_fn)
            actual_block = matchings[z][m]
            assert expected_block == actual_block

    for gamma in gamma_set:
        z = rs_fn(gamma)
        expected_brion = set()
        for m in matchings[z]:
            if is_aligned_fn(gamma, m):
                expected_brion |= matchings[z][m]
        actual_brion = brion_fn(gamma)
        assert expected_brion == actual_brion

