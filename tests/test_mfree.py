from clans import Clan
from permutations import Permutation
from signed import SignedPermutation
from qp_utils import (
    is_quasiparabolic, 
    is_zero_hecke_module,
    is_q_hecke_module
)
from tests.test_crystals import draw_graph
from tests.test_kn import subsets


def cover_group(c):
    g = {i for i in c.generators() if c.weak_order_action(i)[1] is not None}
    bns = [(c.simple_generator(i), c.weak_order_action(i)[0]) for i in g]

    add = {c.weyl_group_identity()}
    ans = set()
    while add:
        newadd = set()
        for w in add:
            ans.add(w)
            for i in g:
                ws = w * c.simple_generator(i)
                newadd.add(ws)
        add = newadd - ans
    return ans, bns


def test_clan_hecke_atoms(n):
    for genfn in [Clan.all_a, Clan.all_b, Clan.all_c1, Clan.all_c2, Clan.all_d1, Clan.all_d2, Clan.all_d3]:
        a = list(genfn(n))
        for c in a:
            atoms = set(c.get_atoms())
            hecke = set(c.get_hecke_atoms())
            if atoms != hecke:
                print(c)
                print()
                for w in atoms:
                    print('  ', w, '=', w.get_reduced_word())
                print()
                for w in hecke - atoms:
                    print('  ', w, '=', w.get_reduced_word())
                print()
                print()
                
                expected = set()
                group, covers = cover_group(c)
                for (s, d) in covers:
                    for w in d.get_hecke_atoms():
                        assert c.weyl_group_length(w * s) == c.weyl_group_length(w) + 1
                        for v in group:
                            h = w * s * v
                            if c.weyl_group_length(h) == c.weyl_group_length(w) + c.weyl_group_length(v) + 1:
                                expected.add(h)
                for w in expected - hecke:
                    print('  extra  ', w, '=', w.get_reduced_word())
                for w in hecke - expected:
                    print('missing  ', w, '=', w.get_reduced_word())
                      
                if expected != hecke:
                    c.draw()
                assert expected == hecke


def get_digraph(poset, length, typename, w0=None):
    if w0 is None:
        w0 = [w for w in poset if length(w) == max(map(length, poset))][0]
    e = w0.get_reduced_word(True) if typename == 'D' else w0.get_reduced_word()
    vertices = {w0}
    edges = set()
    q = [(e, w0, None)]
    seen = set()
    while q:
        e, w, j = q[0]
        q = q[1:]
        for i in range(0, len(e)):
            f = e[:i] + (None,) + e[i + 1:]
            ee = [a for a in e if a is not None]
            ff = [a for a in f if a is not None]
            v = w0.from_word(*ff) if typename == 'A' else w0.from_word(w0.rank, *ff) if typename == 'BC' else w0.from_dword(w0.rank, *ff)
            if length(v) == len(ee) - 1 and v in poset and (v, w) not in seen:
                vertices.add(v)
                edges.add((i, w, v))
                q.append((f, v, i))
                seen.add((v, w))
    return vertices, edges


def test_atom_subsets(n=4):
    types = ['A']#, 'BC', 'D']
    for t in types:
        print('check type', t, n)
        
        longest_element = lambda x: Permutation.longest_element(x) if t == 'A' else SignedPermutation.longest_element(x) if t == 'BC' else SignedPermutation.dtype_longest_element(x)
        length = lambda w: w.dlength() if t == 'D' else w.length()
        group = lambda x: Permutation.all(x) if t == 'A' else SignedPermutation.all(x, dtype=(t=='D'))
        invol = lambda x: Permutation.involutions(x) if t == 'A' else SignedPermutation.involutions(x, dtype=(t=='D'))
        leq = lambda x,y: x.dbruhat_less_equal(y) if t == 'D' else x.strong_bruhat_less_equal(y)

        for z in invol(n):
            atoms = set(z.get_atoms_d() if t == 'D' else z.get_atoms())
            for s in subsets(atoms):
                if len(s) == 0:
                    continue
                upper_poset = get_upper_poset(None, lambda x: s, lambda x: group(n), leq)
                mobius = get_mobius(upper_poset)
                print('  ', z, s, set(mobius.values()))
                try:
                    if False and len(s) == len(atoms):
                        vertices, edges = get_digraph(upper_poset, length, t)
                        draw_graph(vertices, edges, printer=lambda s: ''.join(['-' if a is None else str(a) for a in s.inverse()]), edge_labels=False)
                        input('?')
                    assert set(mobius.values()).issubset({-1, 0, 1})
                except:
                    for x in mobius:
                        if mobius[x] not in {-1, 0, 1}:
                            vertices, edges = get_digraph(upper_poset, length, t, x)
                            draw_graph(vertices, edges, printer=lambda s: ''.join(['-' if a is None else str(a) for a in s.inverse()]), edge_labels=False)
                            input('')


def test_subsets(n=4):
    types = ['A', 'BC', 'D']
    for t in types:
        print('check type', t, n)
        
        length = lambda w: w.dlength() if t == 'D' else w.length()
        group = lambda x: Permutation.all(x) if t == 'A' else SignedPermutation.all(x, dtype=(t=='D'))
        leq = lambda x,y: x.dbruhat_less_equal(y) if t == 'D' else x.strong_bruhat_less_equal(y)

        atoms = lambda n: [w for w in group(n) if (w.is_atom_d(True) if t == 'D' else w.is_atom())]
        
        bylength = {}
        for w in group(n):
            bylength[length(w)] = bylength.get(length(w), set()) | {w}
        mapping = {}
        for ell in sorted(bylength):
            for s in subsets(bylength[ell]):
                upper_poset = get_upper_poset(None, lambda x: s, lambda x: group(n), leq)
                mobius = get_mobius(upper_poset)
                mapping[tuple(sorted(s))] = mobius
                print('  ', ell, s, set(mobius.values()))
                if not set(mobius.values()).issubset({-1, 0, 1}):
                    print('')
                else:
                    for r in subsets(s):
                        mobius = mapping[tuple(sorted(r))]
                        try:
                            assert set(mobius.values()).issubset({-1, 0, 1})
                        except:
                            print()
                            print(r, s)
                            return
                #assert set(mobius.values()).issubset({-1, 0, 1})
                



def get_pseudo_hecke_atoms(m, simple, length, conjugate, translate=None):
    return Clan.get_pseudo_hecke_atoms(m, simple, length, conjugate, translate)
         

def _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn, leq=None):
    count = 0
    subcount = 0
    subsubcount = 0
    of = 0
    seen = set()
    for y in elements:
        if y in seen:
            continue
        of += 1

        xset = construct(y)
        k = min([length(x) for x in xset])
        seen |= xset
        
        height = lambda x: (length(x) - k) // 2

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_q_hecke_module(xset, height, right_action, simple, False)
        b3 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        subsubcount += int(b3)

        #if b3:
        #    assert b2        
        #if b2:
        #    assert b1
        if b1 or b2 or b3:
            minimum = [x for x in xset if height(x) == 0]
            assert test_all or len(minimum) == 1
            m = minimum[0]
            print(printer(m), len(printer(m)**2) == 0, ':', len(minimum), 'of', len(xset), b1, b2, b3)
    
            if test_mobius:
                hecke = get_pseudo_hecke_atoms(m, simple, length, right_action)
                atoms = {z: {w for w in hecke[z] if length(w) == min(map(length, hecke[z]))} for z in hecke if z is not None}
                for z in atoms:
                    upper_poset = get_upper_poset(z, lambda x: atoms[x], allfn, leq)
                    mobius = get_mobius(upper_poset)
                    actual = {x: (-1)**(length(x) - length(next(iter(atoms[z])))) for x in hecke[z]}
                    expected = {}
                    for x in upper_poset:
                        if mobius[x] != 0:
                            expected[x] = mobius[x]
                    assert actual == expected

    print()
    print('count:', count, '>=', subcount, '==', subsubcount, 'of', of)
    print()
    # assert subcount == subsubcount


def test_a_conjugation_modules(n=4, test_mobius=True, test_all=False):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    reflections = [Permutation.t_ij(i, j) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda w: w.length()
    elements = Permutation.involutions(n) if not test_all else Permutation.all(n)
    construct = lambda y: y.conjugacy_class(n)
    right_action = lambda x,s: s * x * s
    left_action = lambda s,x: s * x * s
    printer = lambda m: m
    allfn = lambda x: Permutation.all(n)
    _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn)


def test_a_twisted_conjugation_modules(n=4, test_mobius=True, test_all=False):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    reflections = [Permutation.t_ij(i, j) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda w: w.length()
    t = Permutation.longest_element(n)
    elements = Permutation.twisted_involutions(n) if not test_all else Permutation.all(n)
    construct = lambda y: y.twisted_conjugacy_class(n)
    right_action = lambda x,s: t * s * t * x * s
    left_action = lambda s,x: t * s * t * x * s
    printer = lambda m: t * m
    allfn = lambda x: Permutation.all(n)
    _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn)


def test_b_conjugation_modules(n=3, test_mobius=True, test_all=False):
    simple = [SignedPermutation.s_i(i, n) for i in range(0, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n + 1) for j in range(i, n + 1)]
    length = lambda w: w.length()
    elements = SignedPermutation.involutions(n) if not test_all else SignedPermutation.all(n)
    construct = lambda y: y.conjugacy_class()
    right_action = lambda x,s: s * x * s
    left_action = lambda s,x: s * x * s
    printer = lambda m: m
    allfn = lambda x: SignedPermutation.all(n)
    _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn)


def test_d_conjugation_modules(n=4, test_mobius=False, test_all=False):
    assert n >= 2
    simple = [SignedPermutation.ds_i(-1, n)] + [SignedPermutation.s_i(i, n) for i in range(1, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda x: x.dlength()
    leq = lambda x,y: x.dbruhat_less_equal(y)
    elements = SignedPermutation.involutions(n, dtype=True) if not test_all else SignedPermutation.all(n, dtype=True)
    construct = lambda y: y.conjugacy_class(dtype=True)
    right_action = lambda x,s: s * x * s
    left_action = lambda s,x: s * x * s
    printer = lambda m: m
    allfn = lambda x: SignedPermutation.all(n, dtype=True)
    _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn, leq)


def test_d_twisted_conjugation_modules(n=4, test_mobius=True, test_all=False):
    assert n >= 2
    simple = [SignedPermutation.ds_i(-1, n)] + [SignedPermutation.s_i(i, n) for i in range(1, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda x: x.dlength()
    leq = lambda x,y: x.dbruhat_less_equal(y)
    t = SignedPermutation.s_i(0, n)
    elements = SignedPermutation.involutions(n, dtype=True, twisted=True) if not test_all else SignedPermutation.all(n, dtype=True)
    construct = lambda y: y.twisted_conjugacy_class()
    right_action = lambda x,s: t * s * t * x * s
    left_action = lambda s,x: t * s * t * x * s
    printer = lambda m: t * m
    allfn = lambda x: SignedPermutation.all(n, dtype=True)
    _test_conjugation_modules(test_mobius, test_all, simple, reflections, length, elements, construct, right_action, left_action, printer, allfn, leq)


def test_d_fpf_module(n=4):
    assert n >= 2
    xset = SignedPermutation.fpf_class(n, dtype=True)
    k = (n + 1) // 2
    height = lambda x: (x.dlength() - k) // 2
    t = SignedPermutation.s_i(0, n) if n % 2 != 0 else SignedPermutation.identity(n)
    action = lambda x,s: t * s * t * x * s
    ds = SignedPermutation.ds_i
    simple = [ds(-1, n)] + [ds(i, n) for i in range(1, n)]
    assert is_zero_hecke_module(xset, height, action, simple)


def test_d_fpf_quasiparabolic(n=4):
    # fails if n is odd
    assert n >= 2
    xset = SignedPermutation.fpf_class(n, dtype=True)
    k = n // 2
    height = lambda x: (x.dlength() - k) // 2
    t = SignedPermutation.s_i(0, n) if n % 2 != 0 else SignedPermutation.identity(n)
    action = lambda s,x: t * s * t * x * s
    ds = SignedPermutation.ds_i
    simple = [ds(-1, n)] + [ds(i, n) for i in range(1, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    assert is_quasiparabolic(xset, height, action, simple, reflections)


def test_b_fpf_module(n=4):
    xset = SignedPermutation.fpf_class(n)
    k = (n + 1) // 2
    height = lambda x: (x.length() - k) // 2
    action = lambda x,s: s * x * s
    simple = [SignedPermutation.s_i(i, n) for i in range(n)]
    assert is_zero_hecke_module(xset, height, action, simple)


def test_b_fpf_quasiparabolic(n=4):
    # fails if n is odd
    xset = SignedPermutation.fpf_class(n)
    k = (n + 1) // 2
    height = lambda x: (x.length() - k) // 2
    action = lambda s,x: s * x * s
    simple = [SignedPermutation.s_i(i, n) for i in range(n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n + 1) for j in range(i, n + 1)]
    assert is_quasiparabolic(xset, height, action, simple, reflections)


def get_upper_poset(z, atomsfn, allfn, leq=None):
    if leq is None:
        leq = lambda a,w: a.strong_bruhat_less_equal(w)

    atoms = atomsfn(z)
    upper_poset = {}
    for w in allfn(z):
        for a in atoms:
            if leq(a, w):
                upper_poset[w] = set()
    for x in upper_poset:
        for y in upper_poset:
            if leq(x, y):
                upper_poset[y].add(x)
    return upper_poset


def get_dinv_poset(z, twisted):
    n = z.rank
    atomsfn = lambda x: x.get_atoms_d(twisted)
    allfn = lambda x: x.all(n, dtype=True)
    leq = lambda x,y: x.dbruhat_less_equal(y)
    return get_upper_poset(z, atomsfn, allfn, leq)


def get_inv_poset(z):
    n = z.rank
    return get_upper_poset(z, lambda x: x.get_atoms(), lambda x: x.all(n))


def get_fpf_poset(z):
    n = z.rank
    return get_upper_poset(z, lambda x: x.get_fpf_atoms(), lambda x: x.all(n))


def get_twisted_poset(z, n):
    return get_upper_poset(z, lambda x: x.get_twisted_atoms(n), lambda x: x.all(n))
                

def get_mobius(upper_poset):
    ans = {}

    def mobius(x):
        if x not in ans:
            mu = 0
            for z in upper_poset[x] - {x}:
                mu += mobius(z)
            ans[x] = 1 - mu
        return ans[x]

    for x in upper_poset:
        mobius(x)
    return ans


def test_twisted_hecke_mobius(n=4):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    length = lambda x: x.length()
    t = Permutation.longest_element(n)
    conjugate = lambda x,s: t * s * t * x * s
    translate = lambda x,s: x * s

    for y in Permutation.twisted_involutions(n):
        yhecke = get_pseudo_hecke_atoms(y, simple, length, conjugate, translate)

        for z in Permutation.twisted_involutions(n):
            if z not in yhecke:
                continue

            hecke = yhecke[z] 
            minlen = min(map(len, hecke))
            atomsfn = lambda z: {x for x in hecke if x.length() == minlen}
            allfn = lambda x: x.all(n)

            upper_poset = get_upper_poset(z, atomsfn, allfn, leq=None)
            mobius = get_mobius(upper_poset)
            actual = {x: (-1)**(x.length() - z.twisted_involution_length(n) + y.twisted_involution_length(n)) for x in hecke}
            expected = {}
            for x in upper_poset:
                if mobius[x] != 0:
                    expected[x] = mobius[x]
                    # print('  ', x, expected[x], 'vs', actual.get(x, '*'))
            print(set(actual) == set(expected), 'y =', y, 'z =', z)
            if len(y) == 0 or set(actual) == set(expected):
                assert actual == expected


def test_fpf_hecke_mobius(n=4):
    # tests formula in Theorem 3 of https://arxiv.org/pdf/0902.1930
    for z in Permutation.fpf_involutions(n):
        print('z =', z)
        upper_poset = get_fpf_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_symplectic_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.fpf_involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_inv_hecke_mobius(n=4):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    length = lambda x: x.length()
    conjugate = lambda x,s: s * x * s
    translate = lambda x,s: x * s

    for y in Permutation.involutions(n):
        yhecke = get_pseudo_hecke_atoms(y, simple, length, conjugate, translate)

        for z in Permutation.involutions(n):
            if z not in yhecke:
                continue

            hecke = yhecke[z] 
            minlen = min(map(len, hecke))
            atomsfn = lambda z: {x for x in hecke if x.length() == minlen}
            allfn = lambda x: x.all(n)

            upper_poset = get_upper_poset(z, atomsfn, allfn, leq=None)
            mobius = get_mobius(upper_poset)
            actual = {x: (-1)**(x.length() - z.involution_length() + y.involution_length()) for x in hecke}
            expected = {}
            for x in upper_poset:
                if mobius[x] != 0:
                    expected[x] = mobius[x]
                    # print('  ', x, expected[x], 'vs', actual.get(x, '*'))
            print(set(actual) == set(expected), 'y =', y, 'z =', z)
            if len(y) == 0 or set(actual) == set(expected):
                assert actual == expected


def test_binv_hecke_mobius(n=4):
    for z in SignedPermutation.involutions(n):
        print('z =', z)
        upper_poset = get_inv_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_dinv_hecke_mobius(n=4):
    for z in SignedPermutation.involutions(n, dtype=True):
        print('z =', z)
        upper_poset = get_dinv_poset(z, False)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms(dtype=True))
        actual = {x: (-1)**(x.length() - z.involution_length(dtype=True)) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        if actual != expected:
            for x in set(actual) | set(expected):
                print('  ', x, actual.get(x, 0), expected.get(x, 0))
        assert actual == expected


def test_dtwisted_hecke_mobius(n=4):
    for z in SignedPermutation.involutions(n, dtype=True, twisted=True):
        print('z =', z)
        upper_poset = get_dinv_poset(z, True)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms(dtype=True, twisted=True))
        actual = {x: (-1)**(x.length() - z.involution_length(dtype=True, twisted=True)) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        if actual != expected:
            for x in set(actual) | set(expected):
                print('  ', x, actual.get(x, 0), expected.get(x, 0))
        assert actual == expected