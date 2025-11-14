from clans import Clan
from permutations import Permutation
from signed import SignedPermutation


def test(slow=False):
    _test_CII(1)
    _test_CII(2)
    _test_CII(3)
    _test_CII(4)
    _test_CII(5)
    _test_CII(6)

    _test_CI(1)
    _test_CI(2)
    _test_CI(3)
    _test_CI(4)
    _test_CI(5)

    _test_BI(1)
    _test_BI(2)
    _test_BI(3)
    _test_BI(4)
    _test_BI(5)

    _test_AIII(1)
    _test_AIII(2)
    _test_AIII(3)
    _test_AIII(4)
    _test_AIII(5)

    _test_AII(1)
    _test_AII(3)
    _test_AII(5)
    _test_AII(7)

    _test_AI(1)
    _test_AI(2)
    _test_AI(3)
    _test_AI(4)
    _test_AI(5)
    _test_AI(6)

    if slow:
        _test_AI(7)
        _test_AI(8)

        _test_AII(9)
        _test_AII(11)
        
        _test_AIII(6)
        _test_AIII(7)

        _test_BI(6)
        _test_BI(7)

        _test_CI(6)
        _test_CI(7)

        _test_CII(6)
        _test_CII(7)


def precsim(k, n):
    def span(w):
        o = [w(i) for i in range(1, n + 1)]
        for i in range(k, len(o) - 2):
            b, c, a = o[i:i + 3]
            if a < b < c:
                new_o = o[:i] + [c, a, b] + o[i + 3:]
                yield w.__class__(*new_o)
    return span


def precsim_AIII(k, n):
    def span(w):
        o = [w(i) for i in range(1, n + 1)]
        for i in range(len(o) - 1):
            j = len(o) - 1 - i
            if (j - 2) - (i + 1) >= k:
                a1, b1, b2, a2 = o[i], o[i + 1], o[j - 1], o[j]
                if a1 < b1 and a2 < b2:
                    new_o = o[:i] + [b1, a1] + o[i + 2:j - 1] + [a2, b2] + o[j + 1:]
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


def _test_AI(rank):
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


def _test_AII(rank):
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


def _test_AIII(rank):
    n = rank + 1
    twisted = set(Permutation.twisted_involutions(n))

    print('testing AIII, rank =', rank)
    for p in range(n + 1):
        q = n - p
        k = abs(p - q)
        print('  ', 'p =', p, 'q =', q)
        
        gamma_set = set(Clan.all_a(p, q))
        invol_set = {z for z in twisted if len(z.twisted_fixed_points(n)) >= k}
        rs_fn = lambda gamma: gamma.richardson_springer_involution()
        brion_fn = lambda gamma: {w.inverse() for w in gamma.get_atoms()}
        extended_brion_fn = lambda z: {w.inverse() for w in z.get_twisted_atoms(n, offset=k)}
        
        base = lambda z: z.twisted_fixed_points(n)
        triv = lambda m: len([(a, b) for (a, b) in m if a + b == 0])
        matchings_fn = lambda z: {m for m in Permutation.ncsp_matchings(base(z)) if triv(m) == k}
        
        shape_fn = lambda w: tuple(sorted(w.inverse().twisted_shape(n, k)))
        is_aligned_fn = lambda gamma, m: gamma.is_aligned(m)
        generator_fn = lambda z, m: z.get_min_twisted_atom(n, m).inverse()
        span_fn = precsim_AIII(k, n)
        _generic_test(gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, matchings_fn, shape_fn, is_aligned_fn, generator_fn, span_fn)


def nest(o):
    ndes = []
    while True:
        i = [i for i in range(len(o) - 1) if o[i] > o[i + 1]]
        if len(i) == 0:
            break
        i = i[0]
        ndes.append((o[i], o[i + 1]))
        o = o[:i] + o[i + 2:]
    nres = o
    return ndes, nres


def cyc_pm(z, m, n):
    ans = [(a, b) for a in range(-n, n + 1) for b in range(1, n + 1) if a != 0 and (abs(a) < z(a) == b or a == z(a) == b)]
    ans += [(-b, a) for (a, b) in m if 0 < a < b]
    return sorted(ans, key=lambda pair: pair[1])


def desword(p):
    ans = []
    for a, b in p:
        if a == b:
            ans.append(a)
        else:
            ans += [b, a]
    return ans


def ascword(p):
    ans = []
    for a, b in p:
        if a == b:
            ans.append(a)
        else:
            ans += [a, b]
    return ans


def _test_BI(rank):
    n = rank
    invol = set(SignedPermutation.involutions(n))

    print('testing BI, rank =', rank)
    for p in range(2 * n + 2):
        q = 2 * n + 1 - p
        k = (abs(p - q) - 1) // 2
        print('  ', 'p =', p, 'q =', q)
        
        gamma_set = {Clan(oneline, Clan.TYPE_B) for oneline in Clan.symmetric_clans(p, q)}
        invol_set = {z for z in invol if len(z.neg()) >= k}
        rs_fn = lambda gamma: gamma.richardson_springer_involution()
        brion_fn = lambda gamma: {w.inverse() for w in gamma.get_atoms()}
        extended_brion_fn = lambda z: {w.inverse() for w in z.get_atoms(offset=k)}
        
        base = lambda z: [i for i in range(1, n + 1) if z(i) == -i]
        triv = lambda m: len([(a, b) for (a, b) in m if a + b == 0])
        matchings_fn = lambda z: {m for m in Permutation.ncsp_matchings(base(z)) if triv(m) >= k}
        
        def shape_fn(w):
            ans = []
            o = [w(i) for i in range(1, n + 1)]
            ndes, nres = nest(o[k:])
            for (a, b) in ndes:
                if 0 < a < -b:
                    ans.append((a, -b))
                    ans.append((b, -a))
            for a in nres:
                if a < 0:
                    ans.append((a, -a))
            for a in o[:k]:
                b = abs(a)
                ans.append((-b, b))
            return tuple(sorted(ans))

        is_aligned_fn = lambda gamma, m: gamma.is_aligned(m)
        
        def generator_fn(z, m):
            c = sorted([b for (a, b) in m if a + b == 0], reverse=True)
            u = list(reversed(c[:k])) + [-b for b in c[k:]]
            v = desword(cyc_pm(z, m, n))
            return SignedPermutation(*(u + v))

        span_fn = precsim(k, n)
        _generic_test(gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, matchings_fn, shape_fn, is_aligned_fn, generator_fn, span_fn)


def _test_CI(rank):
    n = rank
    print('testing CI, rank =', rank)

    gamma_set = set(Clan.all_c1(n))
    invol_set = set(SignedPermutation.involutions(n))
    rs_fn = lambda gamma: gamma.richardson_springer_involution()
    brion_fn = lambda gamma: {w.inverse() for w in gamma.get_atoms()}
    extended_brion_fn = lambda z: {w.inverse() for w in z.get_atoms()}
    
    base = lambda z: [i for i in range(1, n + 1) if z(i) == -i]
    matchings_fn = lambda z: set(Permutation.ncsp_matchings(base(z)))
    
    def shape_fn(w):
        ans = []
        o = [w(i) for i in range(1, n + 1)]
        ndes, nres = nest(o)
        for (a, b) in ndes:
            if 0 < a < -b:
                ans.append((a, -b))
                ans.append((b, -a))
        for a in nres:
            if a < 0:
                ans.append((a, -a))
        return tuple(sorted(ans))

    is_aligned_fn = lambda gamma, m: gamma.is_aligned(m)
    
    def generator_fn(z, m):
        c = sorted([b for (a, b) in m if a + b == 0], reverse=True)
        u = [-b for b in c]
        v = desword(cyc_pm(z, m, n))
        return SignedPermutation(*(u + v))

    span_fn = precsim(0, n)
    _generic_test(gamma_set, invol_set, rs_fn, brion_fn, extended_brion_fn, matchings_fn, shape_fn, is_aligned_fn, generator_fn, span_fn)


def _test_CII(rank):
    n = rank
    invol = set(SignedPermutation.fpf_involutions(n))

    print('testing CII, rank =', rank)
    for p in range(0, 2 * n + 1, 2):
        q = 2 * n - p
        k = abs(p - q) // 2
        print('  ', 'p =', p, 'q =', q)
        
        gamma_set = {Clan(oneline, Clan.TYPE_C2) for oneline in Clan.symmetric_clans(p, q, disallow_negations=True)}
        invol_set = {z for z in invol if len(z.neg()) >= k}
        rs_fn = lambda gamma: gamma.richardson_springer_involution()
        brion_fn = lambda gamma: {w.inverse() for w in gamma.get_atoms()}
        extended_brion_fn = lambda z: {w.inverse() for w in z.get_fpf_atoms(offset=k)}
        
        base = lambda z: [i for i in range(1, n + 1) if z(i) == -i]
        triv = lambda m: len([(a, b) for (a, b) in m if a + b == 0])
        matchings_fn = lambda z: {m for m in Permutation.ncsp_matchings(base(z)) if triv(m) == k}
        
        def shape_fn(w):
            ans = []
            o = [w(i) for i in range(1, n + 1)]
            for a in o[:k]:
                b = abs(a)
                ans.append((-b, b))
            o = o[k:]
            while o:
                b, c = o[:2]
                if 0 < c < -b:
                    ans.append((c, -b))
                    ans.append((b, -c))
                o = o[2:]
            return tuple(sorted(ans))

        is_aligned_fn = lambda gamma, m: gamma.is_aligned(m)
        
        def generator_fn(z, m):
            c = sorted([b for (a, b) in m if a + b == 0])
            u = c
            v = ascword(cyc_pm(z, m, n))
            return SignedPermutation(*(u + v))

        span_fn = precapprox(k, n)
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

        print('z =', z)
        expected_matchings = set(matchings_fn(z))
        actual_matchings = {shape_fn(w) for w in extended_brion}
        assert expected_matchings == actual_matchings

        matchings[z] = {m: set() for m in actual_matchings}
        for w in extended_brion:
            m = shape_fn(w)
            matchings[z][m].add(w)

        for m in matchings[z]:
            print('  ', 'm =', m)
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

