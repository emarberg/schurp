from clans import Clan
from permutations import Permutation
from signed import SignedPermutation


def _test_hecke_atoms(cl):
    print(cl)
    atoms = cl.get_atoms()
    for w in cl.get_hecke_atoms():
        print('  ', w, 'atom' if w in atoms else '')
    print()


def test_hecke_atoms_a(n=4):
    for p in range(n + 1):
        q = n - p
        for cl in Clan.all_a(p, q):
            _test_hecke_atoms(cl)

def test_hecke_atoms_b(n=4):
    for p in range(n + 1):
        q = n - p
        for cl in Clan.all_b(p, q):
            _test_hecke_atoms(cl)


def test_hecke_atoms_c2(n=4):
    for p in range(n + 1):
        q = n - p
        for cl in Clan.all_c2(p, q):
            _test_hecke_atoms(cl)


def test_hecke_atoms_c1(n=4):
    for cl in Clan.all_c1(n):
        _test_hecke_atoms(cl)


def test_init():
    w = Clan([1, True, 1, False])
    assert w.cycles() == [(1, 3)]
    assert w(1) == 3
    assert w(2) is True
    assert w(3) == 1
    assert w(4) is False


def test_multiply():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])

    assert 1 * w == w
    assert 3 * w == v

    assert w.richardson_springer_map() == Permutation(3, 2, 1, 4)
    assert v.richardson_springer_map() == Permutation(4, 2, 3, 1)


def test_multiply_b():
    w = Clan([True, True, False, False, False, True, True], Clan.TYPE_B)
    v = Clan([True, 1, 1, False, 2, 2, True], Clan.TYPE_B)
    u2 = Clan([1, True, 1, False, 2, True, 2], Clan.TYPE_B)
    u = Clan([True, 1, 2, False, 1, 2, True], Clan.TYPE_B)
    t = Clan([1, True, 2, False, 1, True, 2], Clan.TYPE_B)
    s = Clan([True, 1, 2, False, 2, 1, True], Clan.TYPE_B)
    r = Clan([1, True, 2, False, 2, True, 1], Clan.TYPE_B)
    q = Clan([1, 2, True, False, True, 2, 1], Clan.TYPE_B)
    p = Clan([1, 2, 3, True, 3, 2, 1], Clan.TYPE_B)

    assert w.richardson_springer_map() == SignedPermutation(1, 2, 3)
    assert v.richardson_springer_map() == SignedPermutation(2, 1, 3)
    assert u.richardson_springer_map() == SignedPermutation(-2, -1, 3)

    assert 0 * w == w
    assert 1 * w == v
    assert 2 * w == w

    assert 0 * v == u
    assert 1 * v == v
    assert 2 * v == u2

    assert 0 * u == u
    assert 1 * u == s
    assert 2 * u == t

    assert 0 * s == s
    assert 1 * s == s
    assert 2 * s == r

    assert 0 * r == r
    assert 1 * r == q
    assert 2 * r == r

    assert 0 * q == p
    assert 1 * q == q
    assert 2 * q == q


def test_multiply_c1():
    w = Clan([True, False, True, False], Clan.TYPE_C1)
    v = Clan([True, 1, 1, False], Clan.TYPE_C1)
    u = Clan([1, 1, 2, 2], Clan.TYPE_C1)

    assert 0 * w == v
    assert 1 * w == u

    assert w.richardson_springer_map() == SignedPermutation(1, 2)
    assert v.richardson_springer_map() == SignedPermutation(-1, 2)
    assert u.richardson_springer_map() == SignedPermutation(2, 1)

    w = u
    u = Clan([1, 2, 1, 2], Clan.TYPE_C1)

    assert 0 * w == u
    assert 1 * w == w

    w = u
    u = Clan([1, 2, 2, 1], Clan.TYPE_C1)

    assert 0 * w == w
    assert 1 * w == u

    w = Clan([True, True, False, False], Clan.TYPE_C1)

    assert 0 * w == v
    assert 1 * w == w

    w = v
    v = Clan([1, True, False, 1], Clan.TYPE_C1)

    assert 0 * w == w
    assert 1 * w == v

    w = v
    v = Clan([1, 2, 2, 1], Clan.TYPE_C1)

    assert 0 * w == v
    assert 1 * w == w


def test_atoms():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])
    assert set(w.get_clan_words()) == {(2, 3)}
    assert set(v.get_clan_words()) == {(2,)}
    assert w.get_atoms() == {Permutation(1, 3, 4, 2)}
    assert v.get_atoms() == {Permutation(1, 3, 2)}


def test_all_a():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])
    a = set(Clan.all_a(2, 2))
    assert v in a and w in a
    assert len(a) == 6 + 12 + 3

    for p in range(4):
        for q in range(4):
            for cl in Clan.all_a(p, q):
                assert cl.clan_type() == p - q


def test_all_b():
    assert len(list(Clan.all_b(2, 1))) == 25

    for p in range(4):
        for q in range(4):
            for cl in Clan.all_b(p, q):
                assert cl.clan_type() == 2 * p - 2 * q - 1


def test_all_c2():
    assert len(list(Clan.all_c2(2, 1))) == 9

    for p in range(4):
        for q in range(4):
            for cl in Clan.all_c2(p, q):
                assert cl.clan_type() == 2 * p - 2 * q


def test_all_c1():
    assert len(list(Clan.all_c1(2))) == 11

    for n in range(4):
        for cl in Clan.all_c1(n):
            assert cl.clan_type() == 0


def test_atoms_a3(n=4):
    w0 = Permutation.longest_element(n)
    for p in range(1, n):
        q = n - p
        offset = (n - abs(p - q)) // 2
        for clan in Clan.all_a(p, q):
            z = w0 * clan.richardson_springer_map()
            atoms = set(clan.get_atoms())
            print(clan)
            for a in atoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            print(z, offset)
            btoms = set(z.get_twisted_atoms(n, offset))
            for a in btoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            assert atoms.issubset(btoms)


def test_atoms_b(n=4):
    for p in range(1, n):
        q = n - p
        offset = abs(2 * p - 2 * q - 1) // 2
        for clan in Clan.all_b(p, q):
            z = -clan.richardson_springer_map()
            atoms = set(clan.get_atoms())
            print(clan)
            for a in atoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            print(z, offset)
            btoms = set(z.get_atoms(offset))
            for a in btoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            assert atoms.issubset(btoms)


def test_atoms_c1(n=4):
    for m in range(n + 1):
        for clan in Clan.all_c1(m):
            z = -clan.richardson_springer_map()
            assert set(clan.get_atoms()).issubset(set(z.get_atoms()))


def test_atoms_c2(n=4):
    for p in range(1, n):
        q = n - p
        offset = abs(p - q)
        for clan in Clan.all_c2(p, q):
            z = -clan.richardson_springer_map()
            atoms = clan.get_atoms()
            print(clan)
            for a in atoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            print(z, offset)
            btoms = set(z.get_fpf_atoms(offset))
            for a in btoms:
                print('  ', a.inverse(), a.inverse().get_reduced_word())
            assert atoms.issubset(btoms)


def relatom_shape_test(aword, y):
    n = y.rank
    yfixed = {i for i in range(1, n + 1) if y(i) == i}
    v = SignedPermutation.identity(n)
    sh = set()
    for a in aword:
        if a > 0 and y(a) == a and y(a + 1) == a + 1:
            e, f = tuple(sorted([v(a), v(a + 1)]))
            sh |= {(e, f), (-f, -e)}
        if a == 0 and y(1) == 1:
            e, f = tuple(sorted([v(-1), v(1)]))
            sh |= {(e, f)}
        s = SignedPermutation.s_i(a, n)
        v *= s
        y = s % y % s
    f = {i for p in sh for i in p}
    return sh | {(-i, i) for i in yfixed - f}


def relatom_shape(w, y):
    return relatom_shape_test(w.get_reduced_word(), y)


def test_relatom_shape(n=3):
    for y in SignedPermutation.involutions(n):
        for w in y.get_atoms():
            sh = w.shape()
            for a in w.inverse().get_reduced_words():
                rsh = relatom_shape_test(a, -y)
                # print(y, w, sh, '?=', rsh)
                assert rsh == sh


def test_relatoms(n=3):
    for k in range(0, n + 1):
        y = SignedPermutation.longest_element(n, k)
        g = SignedPermutation.grassmannian_element(n, k)
        for z in SignedPermutation.involutions(n):
            for w in SignedPermutation.relative_atoms(y, z):
                sh = relatom_shape(w.inverse(), -z)
                gw = g * w
                print(k, y, g, '-z =', -z, 'w =', w, sh, gw.shape(), w.inverse().get_reduced_word())
                assert gw in z.get_atoms()
                assert gw.shape() == sh


def test_atoms_a_refined(n=4, verbose=False):
    w0 = Permutation.longest_element(n)
    for p in range(n + 1):
        q = n - p
        k = (n - abs(p - q)) // 2
        for clan in Clan.all_a(p, q):
            z = clan.richardson_springer_map()
            base = z.fixed(n)
            excluded_guess = {
                m for m in Permutation.nc_matchings(base)
                if not clan.is_aligned(m, verbose=verbose)
            }
            if verbose:
                print('\nexcluded_guess:', excluded_guess)

            atoms_by_shape = {}
            for w in Permutation.get_twisted_atoms(w0 * z, n, k):
                sh = w.twisted_shape(n)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            if verbose:
                print('atoms_by_shape:', atoms_by_shape)

            _test_refinement(clan, atoms_by_shape, excluded_guess)


def test_atoms_b_refined(n=4):
    for p in range(n + 1):
        q = n - p
        k = abs(2 * p - 2 * q - 1) // 2
        y = SignedPermutation.longest_element(n, k)
        for clan in Clan.all_b(p, q):
            z = -clan.richardson_springer_map()
            base = z.negated_points()
            excluded_guess = {
                m for m in SignedPermutation.ncsp_matchings(base)
                if not clan.is_aligned(m, True)
            }
            print('\nexcluded_guess:', excluded_guess)

            atoms_by_shape = {}
            for w in SignedPermutation.relative_atoms(y, z):
                sh = relatom_shape(w.inverse(), -z)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            print('atoms_by_shape:', atoms_by_shape)

            _test_refinement(clan, atoms_by_shape, excluded_guess)


def test_atoms_c2_refined(n=4):
    for p in range(n + 1):
        q = n - p
        k = abs(p - q)
        for clan in Clan.all_c2(p, q):
            z = -clan.richardson_springer_map()
            base = z.negated_points()
            excluded_guess = {
                m for m in SignedPermutation.ncsp_matchings(base)
                if not clan.is_aligned(m)
            }
            print('\nexcluded_guess:', excluded_guess)

            atoms_by_shape = {}
            for w in z.get_fpf_atoms(offset=k):
                sh = w.fpf_shape(offset=k)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            print('atoms_by_shape:', atoms_by_shape)

            _test_refinement(clan, atoms_by_shape, excluded_guess)


def test_atoms_c1_refined(n=4):
    for clan in Clan.all_c1(n):
        z = -clan.richardson_springer_map()
        base = z.negated_points()
        excluded_guess = {
            m for m in SignedPermutation.ncsp_matchings(base)
            if not clan.is_aligned(m)
        }
        z_atoms_by_shape = z.get_atoms_by_shape()

        _test_refinement(clan, z_atoms_by_shape, excluded_guess)


def _test_refinement(clan, atoms_by_shape, excluded_guess):
    clan_atoms = set(clan.get_atoms())

    lines = []
    lines += [str(clan)]
    lines += ['']

    excluded = list(set(excluded_guess) - set(atoms_by_shape))
    otherset = set()
    for sh, subset in atoms_by_shape.items():
        if subset.issubset(clan_atoms):
            lines += ['   ' + str(set(sh))]
            otherset |= subset
        else:
            assert len(subset & clan_atoms) == 0
            excluded.append(sh)
    assert otherset == clan_atoms

    lines += ['']
    if excluded:
        lines += ['EXCLUDED:']
        lines += ['']
        for sh in excluded:
            lines += ['   ' + str(set(sh))]
        lines += ['']
    lines += ['']

    extra = set(excluded_guess) - set(excluded)
    if extra:
        lines += ['NOT EXCLUDED BUT PREDICTED TO BE EXCLUDED:']
        lines += ['']
        for sh in extra:
            lines += ['   ' + str(set(sh))]
        lines += ['']
    lines += ['']

    if set(excluded_guess) != set(excluded):
        print('\n', 'clan atoms:', clan_atoms, '\n')
        print('\n'.join(lines))

    assert set(excluded_guess) == set(excluded)
