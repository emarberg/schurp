from permutations import Permutation
from signed import SignedPermutation
from qp_utils import is_quasiparabolic, is_zero_hecke_module


def get_pseudo_hecke_atoms(m, simple, length, act):
    one = m * m.inverse()
    action = {one: m}
    level = {one}
    while level:
        newlevel = set()
        for w in level:
            z = action[w]
            for s in simple:
                ws = w * s
                if length(ws) > length(w):
                    newlevel.add(ws)
                    szs = None if z is None else act(z, s)
                    if szs is None or length(szs) == length(z):
                        action[ws] = None
                    elif length(szs) < length(z):
                        action[ws] = z
                    else:
                        action[ws] = szs
            level = newlevel
    ans = {}
    for (w, z) in action.items():
        if z not in ans:
            ans[z] = set()
        ans[z].add(w)
    return ans
         

def test_a_conjugation_modules(n):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    reflections = [Permutation.t_ij(i, j) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda w: w.length()

    count = 0
    subcount = 0
    of = 0
    seen = set()
    for y in Permutation.involutions(n):
        if y in seen:
            continue
        of += 1

        xset = y.conjugacy_class(n)
        k = min([length(x) for x in xset])
        seen |= xset
        
        height = lambda x: (length(x) - k) // 2
        right_action = lambda x,s: s * x * s
        left_action = lambda s,x: s * x * s

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        
        if b2:
            assert b1
        if b1:
            minimum = [x for x in xset if height(x) == 0]
            assert len(minimum) == 1
            m = minimum[0]
            print(m.cycle_repr(), b1, b2)
    
            hecke = get_pseudo_hecke_atoms(m, simple, length, right_action)
            atoms = {z: {w for w in hecke[z] if length(w) == min(map(length, hecke[z]))} for z in hecke if z is not None}
            for z in atoms:
                upper_poset = get_upper_poset(z, lambda x: atoms[x], lambda x: z.all(n))
                mobius = get_mobius(upper_poset)
                actual = {x: (-1)**(length(x) - length(next(iter(atoms[z])))) for x in hecke[z]}
                expected = {}
                for x in upper_poset:
                    if mobius[x] != 0:
                        expected[x] = mobius[x]
                assert actual == expected

    print()
    print('count:', count, '>', subcount, 'of', of)
    print()


def test_a_twisted_conjugation_modules(n):
    simple = [Permutation.s_i(i) for i in range(1, n)]
    reflections = [Permutation.t_ij(i, j) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda w: w.length()
    t = Permutation.longest_element(n)

    count = 0
    subcount = 0
    of = 0
    seen = set()
    for y in Permutation.twisted_involutions(n):
        if y in seen:
            continue
        of += 1

        xset = y.twisted_conjugacy_class(n)
        k = min([length(x) for x in xset])
        seen |= xset
        
        height = lambda x: (length(x) - k) // 2

        right_action = lambda x,s: t * s * t * x * s
        left_action = lambda s,x: t * s * t * x * s

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        
        if b2:
            assert b1
        if b1:
            minimum = [x for x in xset if height(x) == 0]
            assert len(minimum) == 1
            m = minimum[0]
            print((t * m).cycle_repr(), b1, b2)
            
            hecke = get_pseudo_hecke_atoms(m, simple, length, right_action)
            atoms = {z: {w for w in hecke[z] if length(w) == min(map(length, hecke[z]))} for z in hecke if z is not None}
            for z in atoms:
                upper_poset = get_upper_poset(z, lambda x: atoms[x], lambda x: z.all(n))
                mobius = get_mobius(upper_poset)
                actual = {x: (-1)**(length(x) - length(next(iter(atoms[z])))) for x in hecke[z]}
                expected = {}
                for x in upper_poset:
                    if mobius[x] != 0:
                        expected[x] = mobius[x]
                assert actual == expected

    print()
    print('count:', count, '>', subcount, 'of', of)
    print()


def test_b_conjugation_modules(n):
    simple = [SignedPermutation.s_i(i, n) for i in range(0, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n + 1) for j in range(i, n + 1)]
    length = lambda w: w.length()

    count = 0
    subcount = 0
    of = 0
    seen = set()
    for y in SignedPermutation.involutions(n):
        if y in seen:
            continue
        of += 1

        xset = y.conjugacy_class()
        k = min([length(x) for x in xset])
        seen |= xset
        
        height = lambda x: (length(x) - k) // 2
        right_action = lambda x,s: s * x * s
        left_action = lambda s,x: s * x * s

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        
        if b2:
            assert b1
        if b1:
            minimum = [x for x in xset if height(x) == 0]
            assert len(minimum) == 1
            m = minimum[0]
            print(m, b1, b2)
            
            hecke = get_pseudo_hecke_atoms(m, simple, length, right_action)
            atoms = {z: {w for w in hecke[z] if length(w) == min(map(length, hecke[z]))} for z in hecke if z is not None}
            for z in atoms:
                upper_poset = get_upper_poset(z, lambda x: atoms[x], lambda x: z.all(n))
                mobius = get_mobius(upper_poset)
                actual = {x: (-1)**(length(x) - length(next(iter(atoms[z])))) for x in hecke[z]}
                expected = {}
                for x in upper_poset:
                    if mobius[x] != 0:
                        expected[x] = mobius[x]
                assert actual == expected

    print()
    print('count:', count, '>', subcount, 'of', of)
    print()


def test_d_conjugation_modules(n):
    assert n >= 2
    simple = [SignedPermutation.ds_i(-1, n)] + [SignedPermutation.s_i(i, n) for i in range(1, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    length = lambda x: x.dlength()
    leq = lambda x,y: x.dbruhat_less_or_equal(y)

    count = 0
    subcount = 0
    of = 0
    seen = set()
    for y in SignedPermutation.involutions(n, dtype=True):
        if y in seen:
            continue
        of += 1

        xset = y.conjugacy_class(dtype=True)
        k = min([length(x) for x in xset])
        seen |= xset
        
        height = lambda x: (length(x) - k) // 2
        right_action = lambda x,s: s * x * s
        left_action = lambda s,x: s * x * s

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        
        if b2:
            assert b1
        if b1:
            minimum = [x for x in xset if height(x) == 0]
            assert len(minimum) == 1
            m = minimum[0]
            print(m, b1, b2)
            
            hecke = get_pseudo_hecke_atoms(m, simple, length, right_action)
            atoms = {z: {w for w in hecke[z] if length(w) == min(map(length, hecke[z]))} for z in hecke if z is not None}
            for z in atoms:
                upper_poset = get_upper_poset(z, lambda x: atoms[x], lambda x: z.all(n))
                mobius = get_mobius(upper_poset)
                actual = {x: (-1)**(length(x) - length(next(iter(atoms[z])))) for x in hecke[z]}
                expected = {}
                for x in upper_poset:
                    if mobius[x] != 0:
                        expected[x] = mobius[x]
                assert actual == expected

    print()
    print('count:', count, '>', subcount, 'of', of)
    print()


def test_d_twisted_conjugation_modules(n):
    assert n >= 2
    simple = [SignedPermutation.ds_i(-1, n)] + [SignedPermutation.s_i(i, n) for i in range(1, n)]
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    t = SignedPermutation.s_i(0, n)

    count = 0
    subcount = 0
    of = 0
    seen = set()
    for y in SignedPermutation.involutions(n, dtype=True, twisted=True):
        if y in seen:
            continue
        of += 1

        xset = y.twisted_conjugacy_class()
        k = min([x.dlength() for x in xset])
        seen |= xset
        
        height = lambda x: (x.dlength() - k) // 2

        right_action = lambda x,s: t * s * t * x * s
        left_action = lambda s,x: t * s * t * x * s

        b1 = is_zero_hecke_module(xset, height, right_action, simple, False)
        b2 = is_quasiparabolic(xset, height, left_action, simple, reflections, False)
        
        count += int(b1)
        subcount += int(b2)
        
        if b1 or b2:
            print([t * x for x in xset if x.dlength() == k], b1, b2)
        if b2:
            assert b1

    print()
    print('count:', count, '>', subcount, 'of', of)
    print()


def test_d_fpf_module(n):
    assert n >= 2
    xset = SignedPermutation.fpf_class(n, dtype=True)
    k = (n + 1) // 2
    height = lambda x: (x.dlength() - k) // 2
    t = SignedPermutation.s_i(0, n) if n % 2 != 0 else SignedPermutation.identity(n)
    action = lambda x,s: t * s * t * x * s
    ds = SignedPermutation.ds_i
    simple = [ds(-1, n)] + [ds(i, n) for i in range(1, n)]
    assert is_zero_hecke_module(xset, height, action, simple)


def test_d_fpf_quasiparabolic(n):
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


def test_b_fpf_module(n):
    xset = SignedPermutation.fpf_class(n)
    k = (n + 1) // 2
    height = lambda x: (x.length() - k) // 2
    action = lambda x,s: s * x * s
    simple = [SignedPermutation.s_i(i, n) for i in range(n)]
    assert is_zero_hecke_module(xset, height, action, simple)


def test_b_fpf_quasiparabolic(n):
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


def test_twisted(n=4):
    for z in Permutation.twisted_involutions(n):
        print('z =', z)
        upper_poset = get_twisted_poset(z, n)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_twisted_hecke_atoms(n))
        actual = {x: (-1)**(x.length() - z.twisted_involution_length(n)) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_fpf(n=4):
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


def test_inv(n=4):
    for z in Permutation.involutions(n):
        print('z =', z)
        upper_poset = get_inv_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
                # print('  ', x, expected[x], 'vs', actual.get(x, '*'))
        assert actual == expected


def test_fpf(n=4):
    for z in SignedPermutation.fpf_involutions(n):
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


def test_binv(n=4):
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