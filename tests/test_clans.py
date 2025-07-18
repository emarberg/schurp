from clans import Clan
from permutations import Permutation
from signed import SignedPermutation
from even import EvenSignedPermutation


def test_action_a(n=3):
    for a in Clan.all_a(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)
            u = a.richardson_springer_map()
            v = b.richardson_springer_map()
            signs_a = a.signs()
            signs_b = b.signs()

            if type(a(i)) == bool and type(a(i + 1)) == bool and a(i) == a(i + 1):
                assert a == b
            else:
                s = a.simple_generator(i)
                assert s % u % s == v
                assert len(signs_a) - len(signs_b) in [0, 2]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def test_action_b(n=3):
    for a in Clan.all_b(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()

            amap = {i - n: a.oneline[i] for i in range(len(a))}
            bmap = {i - n: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            if type(amap[i]) == bool and type(amap[i + 1]) == bool and amap[i] == amap[i + 1]:
                assert a == b
            else:
                s = a.simple_generator(i)
                assert s % u % s == v
                assert len(signs_a) - len(signs_b) in [0, 2, 4]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                elif len(signs_a) == len(signs_b) + 4:
                    assert b.covers(a)
                else:
                    assert b.covers(a, n)


def test_action_c1(n=3):
    for a in Clan.all_c1(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()

            amap = {i - n + 1: a.oneline[i] for i in range(len(a))}
            bmap = {i - n + 1: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            if type(amap[i]) == bool and type(amap[i + 1]) == bool and amap[i] == amap[i + 1]:
                assert a == b
            else:
                s = a.simple_generator(i)
                assert s % u % s == v
                assert len(signs_a) - len(signs_b) in [0, 2, 4]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def test_action_c2(n=3):
    for a in Clan.all_c2(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()

            amap = {i - n + 1: a.oneline[i] for i in range(len(a))}
            bmap = {i - n + 1: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            if type(amap[i]) == bool and type(amap[i + 1]) == bool and amap[i] == amap[i + 1]:
                assert a == b
            else:
                s = a.simple_generator(i)
                w = s % u % s
                if w.negated_points():
                    assert a == b
                else:
                    assert v == w
                assert len(signs_a) - len(signs_b) in [0, 4]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def test_action_d1(n=3):
    for a in Clan.all_d1(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()
            
            u = EvenSignedPermutation(*u)
            v = EvenSignedPermutation(*v)

            amap = {i - n + 1: a.oneline[i] for i in range(len(a))}
            bmap = {i - n + 1: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            if type(amap[i]) == bool and type(amap[i + 1]) == bool and amap[i] == amap[i + 1]:
                assert a == b
            else:
                s = EvenSignedPermutation.s_i(max(i, 0), n)
                assert v == s % u % s
                assert len(signs_a) - len(signs_b) in [0, 4]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def test_action_d2(n=3):
    for a in Clan.all_d2(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()
            
            u = EvenSignedPermutation(*u)
            v = EvenSignedPermutation(*v)

            amap = {i - n + 1: a.oneline[i] for i in range(len(a))}
            bmap = {i - n + 1: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            if type(amap[i]) == bool and type(amap[i + 1]) == bool and amap[i] == amap[i + 1]:
                assert a == b
            else:
                s = EvenSignedPermutation.s_i(max(i, 0), n)
                assert v == s % u % s.star()
                assert len(signs_a) - len(signs_b) in [0, 4]
                if len(signs_a) == len(signs_b):
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def test_action_d3(n=3):
    for a in Clan.all_d3(n):
        for i in a.generators():
            b, _ = a.weak_order_action(i)

            u = a.richardson_springer_map()
            v = b.richardson_springer_map()
            
            u = EvenSignedPermutation(*u)
            v = EvenSignedPermutation(*v)

            amap = {i - n + 1: a.oneline[i] for i in range(len(a))}
            bmap = {i - n + 1: b.oneline[i] for i in range(len(b))}

            signs_a = a.signs()
            signs_b = b.signs()

            j = i + (1 if i > 0 else 2)
            if type(amap[i]) == bool and type(amap[j]) == bool and amap[i] == amap[j]:
                assert a == b
            else:
                s = EvenSignedPermutation.s_i(max(i, 0), n)
                w = s % u % s
                if w.negated_points():
                    assert a == b
                else:
                    assert v == w
                assert len(signs_a) - len(signs_b) in [0, 4]
                if len(signs_a) == len(signs_b):
                    if i == -1 and u != v and (type(amap[1]) == bool or type(amap[2]) == bool):
                        k = len(signs_a) // 2
                        signs_a = list(signs_a)
                        signs_a[k - 1], signs_a[k] = signs_a[k], signs_a[k - 1]
                        signs_a = tuple(signs_a)
                    assert signs_a == signs_b
                else:
                    assert b.covers(a)


def _test_hecke_atoms(cl, dtype=False, verbose=False):
    length = lambda x: cl.weyl_group_length(x)
    get_reduced_word = lambda w: w.get_reduced_word(dtype=dtype) if dtype else w.get_reduced_word()

    atoms = set(cl.get_atoms())
    hecke = set(cl.get_hecke_atoms())
    extended = set(cl.get_hecke_atoms_extended())

    print('clan =', cl)
    print()

    for w in sorted(extended, key=length):
        shapes = set()
        print('   b =', cl.richardson_springer_base())
        print('   z =', cl.richardson_springer_involution())
        print('   w =', w.inverse(), 'atom' if w in atoms else 'hecke atom' if w in hecke else 'EXTRA')
        print()
        for v in atoms:
            if cl.weyl_group_bruhat_leq(v, w):
                sh = cl.weyl_group_shape(v)
                shapes.add(sh)
                print('  ', 'v =', v.inverse(), 'sh =', sh)
        print()
        print('  possible shapes:', len(shapes))
        print()
        print()

    expected = {w for w in cl.get_hecke_atoms_extended() if any(cl.weyl_group_bruhat_leq(v, w) for v in atoms)}
    if verbose:
        print(' extended:', {get_reduced_word(w) for w in extended})
        print(' computed:', {get_reduced_word(w) for w in hecke})
        print('predicted:', {get_reduced_word(w) for w in expected})
        print()
    print(hecke == expected, ': is strongly alternating?', cl.is_strongly_alternating())
    print()

    if not hecke.issubset(expected):
        print({get_reduced_word(w) for w in hecke - expected})
    assert hecke.issubset(expected)

    if cl.is_strongly_alternating():
        assert hecke == expected


def test_hecke_atoms_a(n=3, verbose=False):
    for cl in Clan.all_a(n):
        _test_hecke_atoms(cl, verbose=verbose)


def test_hecke_atoms_b(n=3, verbose=False):
    for cl in Clan.all_b(n):
        _test_hecke_atoms(cl, verbose=verbose)


def test_hecke_atoms_c1(n=3, verbose=False):
    for cl in Clan.all_c1(n):
        _test_hecke_atoms(cl, verbose=verbose)


def test_hecke_atoms_c2(n=3, verbose=False):
    for cl in Clan.all_c2(n):
        _test_hecke_atoms(cl, verbose=verbose)


def test_hecke_atoms_d1(n=3, verbose=False):
    for cl in Clan.all_d1(n):
        _test_hecke_atoms(cl, dtype=True, verbose=verbose)


def test_hecke_atoms_d2(n=3, verbose=False):
    for cl in Clan.all_d2(n):
        _test_hecke_atoms(cl, dtype=True, verbose=verbose)


def test_hecke_atoms_d3(n=3, verbose=False):
    for cl in Clan.all_d3(n):
        _test_hecke_atoms(cl, dtype=True, verbose=verbose)


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


def test_atoms_a(n=4):
    for p in range(1, n):
        q = n - p
        atoms_by_inv = {}
        extended_atoms = {}
        for clan in Clan.all_a(p, q):
            atoms = set(clan.get_atoms())
            for a in atoms:
                assert clan.weyl_group_weight(a) == 0
            btoms = set(clan.get_atoms_extended())
            assert atoms.issubset(btoms)
            assert (atoms == btoms) == clan.is_alternating()

            z = clan.richardson_springer_involution()
            atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
            if z not in extended_atoms:
                extended_atoms[z] = btoms
            else:
                assert extended_atoms[z] == btoms

        unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
        print('p =', p, 'q =', q, len(atoms_by_inv), len(unions_match))
        assert abs(p - q) > 1 or len(atoms_by_inv) == len(unions_match)

def test_atoms_b(n=4, verbose=False):
    for p in range(1, n):
        q = n - p
        k = abs(2 * p - 2 * q - 1) // 2
        g = SignedPermutation.bbase_atom(n, k)

        atoms_by_inv = {}
        extended_atoms = {}

        for clan in Clan.all_b(p, q):
            z = clan.richardson_springer_involution()
            d = len([i for i in range(1, n + 1) if i > z(i) > 0])
            atoms = set(clan.get_atoms())
            lengths_a = set()
            weights = set()
            for a in atoms:
                word = a.inverse().get_reduced_word()
                lengths_a.add(len(word))
                assert clan.weyl_group_weight(a) == (g * a).brion_length_b() - g.brion_length_b()
                assert clan.weyl_group_weight(a) == d + a.ell_zero()
                weights.add(d + a.ell_zero())
            #if weights == {0}:
            #    print(z)

            btoms = set(clan.get_atoms_extended())
            lengths_b = set()
            for a in btoms:
                word = a.inverse().get_reduced_word()
                lengths_b.add(len(word))

            assert atoms.issubset(btoms)
            assert lengths_a == lengths_b
            assert len(lengths_a) == 1
            assert (atoms == btoms) == clan.is_alternating()

            z = clan.richardson_springer_involution()
            atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
            if z not in extended_atoms:
                extended_atoms[z] = btoms
            else:
                assert extended_atoms[z] == btoms

        unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
        print('p =', p, 'q =', q, len(atoms_by_inv), len(unions_match))
        assert abs(p - q) > 1 or len(atoms_by_inv) == len(unions_match)


def test_atoms_c1(n=4):
    for m in range(n + 1):
        atoms_by_inv = {}
        extended_atoms = {}
        for clan in Clan.all_c1(m):
            z = clan.richardson_springer_involution()
            d = len([i for i in range(1, n + 1) if i > z(i)])
            atoms = set(clan.get_atoms())
            btoms = set(clan.get_atoms_extended())
            assert atoms.issubset(btoms)
            for a in atoms:
                assert clan.weyl_group_weight(a) == d - a.ell_zero()
            assert (atoms == btoms) == clan.is_alternating()

            atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
            if z not in extended_atoms:
                extended_atoms[z] = btoms
            else:
                assert extended_atoms[z] == btoms

        unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
        print('n =', m, len(atoms_by_inv), len(unions_match))
        assert len(atoms_by_inv) == len(unions_match)


def test_atoms_c2(n=4):
    for p in range(1, n):
        q = n - p
        
        atoms_by_inv = {}
        extended_atoms = {}

        for clan in Clan.all_c2(p, q):
            atoms = clan.get_atoms()
            for a in atoms:
                assert clan.weyl_group_weight(a) == 0
            btoms = set(clan.get_atoms_extended())
            assert atoms.issubset(btoms)
            assert (atoms == btoms) == clan.is_alternating()

            z = clan.richardson_springer_involution()
            atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
            if z not in extended_atoms:
                extended_atoms[z] = btoms
            else:
                assert extended_atoms[z] == btoms

        unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
        print('p =', p, 'q =', q, len(atoms_by_inv), len(unions_match))
        assert abs(p - q) > 1 or len(atoms_by_inv) == len(unions_match)


def test_atoms_d1(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        for p in range(1, n):
            q = n - p
            offset = abs(p - q)
            twisted = offset % 2 != 0
            g = SignedPermutation.dbase_atom(n, twisted, offset)

            print('n =', n, '(p, q) =', (p, q))
            atoms_by_inv = {}
            extended_atoms = {}

            for clan in Clan.all_d1(p, q):
                z = clan.richardson_springer_involution()
                if twisted:
                    z = SignedPermutation.s_i(0, n) * z
                d = len([i for i in range(-n, n + 1) if i > z(i)]) - offset
                assert d % 2 == 0
                    
                atoms = set(clan.get_atoms())
                lengths_a = set()
                weights = set()
                for a in atoms:
                    word = a.inverse().get_reduced_word(dtype=True)
                    lengths_a.add(len(word))
                    assert g.brion_length_d(twisted) == offset // 2
                    assert clan.weyl_group_weight(a) == (g * a).brion_length_d(twisted) - g.brion_length_d(twisted)
                    assert clan.weyl_group_weight(a) == d // 2
                    weights.add(d // 2)
                #if weights == {0}:
                #    print(z)
                
                btoms = set(clan.get_atoms_extended())
                lengths_b = set()
                for a in btoms:
                    word = a.inverse().get_reduced_word(dtype=True)
                    lengths_b.add(len(word))

                assert atoms.issubset(btoms)
                assert lengths_a == lengths_b
                assert len(lengths_a) == 1
                assert (atoms == btoms) == clan.is_alternating()

                z = clan.richardson_springer_involution()
                atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
                if z not in extended_atoms:
                    extended_atoms[z] = btoms
                else:
                    assert extended_atoms[z] == btoms

            unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
            print('p =', p, 'q =', q, len(atoms_by_inv), len(unions_match))
            assert abs(p - q) > 1 or len(atoms_by_inv) == len(unions_match)


def test_atoms_d2(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        for p in range(1, n + 1):
            q = n + 1 - p
            offset = abs(p - q)
            twisted = offset % 2 != 0
            g = SignedPermutation.dbase_atom(n, twisted, offset)

            print('n =', n, '(p, q) =', (p, q))

            atoms_by_inv = {}
            extended_atoms = {}

            equal = set()
            every = set(Clan.all_d2(p, q))
            extra = set()
            for clan in every:
                z = clan.richardson_springer_involution()
                if twisted:
                    z = SignedPermutation.s_i(0, n) * z
                d = len([i for i in range(-n, n + 1) if i > z(i)]) - offset
                assert d % 2 == 0

                atoms = set(clan.get_atoms())
                for a in atoms:
                    assert clan.weyl_group_weight(a) == (g * a).brion_length_d(twisted) - g.brion_length_d(twisted)
                    assert clan.weyl_group_weight(a) == d // 2
                btoms = set(clan.get_atoms_extended())
                assert atoms.issubset(btoms)
                assert (atoms == btoms) == clan.is_alternating()

                z = clan.richardson_springer_involution()
                if atoms == btoms:
                    equal.add(z)
                    # print(clan, clan.is_alternating())
                else:
                    extra.add(clan)

                atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
                if z not in extended_atoms:
                    extended_atoms[z] = btoms
                else:
                    assert extended_atoms[z] == btoms

            #print()
            #print('extra:')
            #for clan in extra:
            #    print(clan, clan.is_alternating())

            unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
            print()
            print('p =', p, 'q =', q, len(atoms_by_inv), len(unions_match), len(equal))
            print()
            assert abs(p - q) > 1 or len(atoms_by_inv) == len(unions_match)


def test_atoms_d3(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        atoms_by_inv = {}
        extended_atoms = {}
        equal = set()
        every = set(Clan.all_d3(n))
        extra = set()
        for clan in every:
            phi = clan.richardson_springer_map()
            z = phi.dtype_longest_element(n) * phi.inverse()

            atoms = set(clan.get_atoms())
            lengths_a = set()
            for a in atoms:
                word = a.inverse().get_reduced_word(dtype=True)
                lengths_a.add(len(word))
                assert clan.weyl_group_weight(a) == 0

            btoms = set(z.get_fpf_atoms_d())
            lengths_b = set()
            for a in btoms:
                word = a.inverse().get_reduced_word(dtype=True)
                lengths_b.add(len(word))

            assert atoms.issubset(btoms)
            assert lengths_a == lengths_b
            assert len(lengths_a) == 1
            assert (atoms == btoms) == clan.is_alternating()

            z = clan.richardson_springer_involution()
            if atoms == btoms:
                equal.add(z)
                # print(clan, clan.is_alternating())
            else:
                extra.add(clan)

            atoms_by_inv[z] = atoms_by_inv.get(z, set()) | atoms
            if z not in extended_atoms:
                extended_atoms[z] = btoms
            else:
                assert extended_atoms[z] == btoms

        # print()
        # print('extra:')
        # for clan in extra:
        #     print(clan, clan.is_alternating())

        unions_match = [z for z in atoms_by_inv if atoms_by_inv[z] == extended_atoms[z]]
        print()
        print('n =', n, len(atoms_by_inv), len(unions_match), len(equal))
        print()
        assert len(atoms_by_inv) == len(unions_match)


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
        k = abs(p - q)
        print('n =', n, '(p, q) =', (p, q), 'k =', k)
        for clan in Clan.all_a(p, q):
            z = clan.richardson_springer_map()
            base = z.fixed(n)
            z = w0 * z
            expected_shapes = {
                m for m in Permutation.ncsp_matchings(base)
                if clan.is_aligned(m, verbose=verbose)
            }
            atoms_by_shape = {}
            for w in Permutation.get_twisted_atoms(z, n, k):
                sh = w.twisted_shape(n, k)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            _test_atype_atoms_by_shape(z, n, k, atoms_by_shape)
            _test_refinement(clan, atoms_by_shape, expected_shapes)


def _test_atype_atoms_by_shape(z, n, k, atoms_by_shape):
    def span(v):
        v = tuple(v(i) for i in range(1, n + 1))
        level = {(None, v)}
        while level:
            nextlevel = set()
            for u, v in level:
                yield Permutation(*v).inverse()

                for i in range(n):
                    if (n - i - 3) - (i + 1) < k:
                        break
                    b1, a1, a2, b2 = v[i], v[i + 1], v[n - i - 2], v[n - i - 1]
                    if a1 < b1 and a2 < b2:
                        w = v[:i] + (a1, b1,) + v[i + 2:n - i - 2] + (b2, a2) + v[n - i:]
                        nextlevel.add((v, w))

            level = nextlevel

    for sh, atoms in atoms_by_shape.items():
        a = z.get_max_twisted_atom(n, sh)
        start = a.inverse()
        btoms = set(span(start))
        if atoms != btoms:
            print('  z =', z, 'a =', a.inverse(), 'k =', k, sh)
            print('  ', {w.inverse() for w in atoms})
            print('  ', {w.inverse() for w in btoms})
        assert atoms == btoms


def test_atoms_b_refined(n=3, verbose=False):
    for p in range(n + 1):
        q = n - p
        k = abs(2 * p - 2 * q - 1) // 2
        y = SignedPermutation.longest_element(n, k)
        g = SignedPermutation.bbase_atom(n, k)

        print('n =', n, '(p, q) =', (p, q), 'k =', k)
        for clan in Clan.all_b(p, q):
            z = clan.richardson_springer_involution()
            base = z.negated_points()
            expected_shapes = {
                m for m in SignedPermutation.ncsp_matchings(base)
                if clan.is_aligned(m)
            }
            atoms_by_shape = {}
            for w in SignedPermutation.relative_atoms(y, z):
                sh = (g * w).shape(k)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            _test_refinement(clan, atoms_by_shape, expected_shapes)


def test_atoms_c2_refined(n=3):
    for p in range(n + 1):
        q = n - p
        k = abs(p - q)

        print('n =', n, '(p, q) =', (p, q), 'k =', k)
        for clan in Clan.all_c2(p, q):
            z = -clan.richardson_springer_map()
            base = z.negated_points()
            expected_shapes = {
                m for m in SignedPermutation.ncsp_matchings(base)
                if clan.is_aligned(m)
            }
            atoms_by_shape = {}
            for w in z.get_fpf_atoms(offset=k):
                sh = w.fpf_shape(offset=k)
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}
            _test_refinement(clan, atoms_by_shape, expected_shapes)


def test_atoms_c1_refined(n=4):
    for clan in Clan.all_c1(n):
        z = -clan.richardson_springer_map()
        base = z.negated_points()
        expected_shapes = {
            m for m in SignedPermutation.ncsp_matchings(base)
            if clan.is_aligned(m)
        }
        z_atoms_by_shape = z.get_atoms_by_shape()
        _test_refinement(clan, z_atoms_by_shape, expected_shapes)


def test_atoms_d1_refined(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        for p in range(0, n + 1):
            q = n - p
            k = abs(p - q)
            twisted = k % 2 != 0
            g = SignedPermutation.dbase_atom(n, twisted, k)
            
            print('n =', n, '(p, q) =', (p, q), 'k =', k)
            
            for clan in Clan.all_d1(p, q):
                z = clan.richardson_springer_involution()
                t = SignedPermutation.s_i(0, n) if twisted else SignedPermutation.identity(n)
                base = (t * z).negated_points()
                
                expected_shapes = {
                    m for m in SignedPermutation.ncsp_matchings(base)
                    if clan.is_aligned(m)
                }
                
                atoms_by_shape = {}
                for w in z.get_atoms_d(twisted=twisted, offset=k):
                    assert (g * w).dlength() == g.dlength() + w.dlength()
                    sh = (g * w).dshape(k)
                    sh = tuple(sorted(sh))
                    atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}

                _test_dtype_atoms_by_shape(z, g, k, atoms_by_shape)
                _test_refinement(clan, atoms_by_shape, expected_shapes)


def test_atoms_d2_refined(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        for p in range(1, n + 1):
            q = n + 1 - p
            k = abs(p - q)
            twisted = k % 2 != 0
            g = SignedPermutation.dbase_atom(n, twisted, k)

            print('n =', n, '(p, q) =', (p, q), 'k =', k)

            for clan in Clan.all_d2(p, q):
                z = clan.richardson_springer_involution()
                t = SignedPermutation.s_i(0, n) if twisted else SignedPermutation.identity(n)
                base = (t * z).negated_points()
                
                expected_shapes = {
                    m for m in SignedPermutation.ncsp_matchings(base)
                    if clan.is_aligned(m)
                }
                
                atoms_by_shape = {}
                for w in z.get_atoms_d(twisted=twisted, offset=k):
                    assert (g * w).dlength() == g.dlength() + w.dlength()
                    sh = (g * w).dshape(k)
                    sh = tuple(sorted(sh))
                    atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}

                _test_dtype_atoms_by_shape(z, g, k, atoms_by_shape)
                _test_refinement(clan, atoms_by_shape, expected_shapes)


def test_atoms_d3_refined(nn=4, verbose=False):
    for n in [nn, nn + 1]:
        g = SignedPermutation.one_fpf_d(n)
        print('n =', n)
        for clan in Clan.all_d3(n):
            z = clan.richardson_springer_involution()
            t = SignedPermutation.s_i(0, n) if (n % 2 != 0) else SignedPermutation.identity(n)
            base = (t * z).negated_points()

            expected_shapes = {
                m for m in SignedPermutation.ncsp_matchings(base)
                if clan.is_aligned(m)
            }

            atoms_by_shape = {}
            for w in z.get_fpf_atoms_d():
                assert (g * w).dlength() == g.dlength() + w.dlength()
                sh = w.fpf_dshape()
                sh = tuple(sorted(sh))
                atoms_by_shape[sh] = atoms_by_shape.get(sh, set()) | {w}

            _test_dtype_fpf_atoms_by_shape(z, atoms_by_shape)
            _test_refinement(clan, atoms_by_shape, expected_shapes)


def _test_dtype_atoms_by_shape(z, g, k, atoms_by_shape):
    def forms(v):
        yield v
        if k >= 1 and len(v) > k:
             b, a = v[0], v[k]
             if abs(a) < abs(b):
                yield (-b,) + v[1:k] + (-a,) + v[k + 1:]

    def span(v):
        v = v.oneline
        level = {(None, x) for x in forms(v)}
        while level:
            nextlevel = set()
            for u, v in level:
                yield SignedPermutation(*v).inverse()

                for i in range(k, len(v) - 2):
                    b, c, a = v[i: i + 3]
                    if a < b < c:
                        for w in forms(v[:i] + (c, a, b) + v[i + 3:]):
                            nextlevel.add((v, w))

                if k == 0 and len(v) >= 3:
                    b, c, a = v[:3]
                    if a < -b < c:
                        w = (-c, a, -b) + v[3:]
                        nextlevel.add((v, w))

            level = nextlevel

    z = EvenSignedPermutation(*z)
    g = EvenSignedPermutation(*g)
    for sh, atoms in atoms_by_shape.items():
        a = z.get_max_atom(sh) if k % 2 ==0 else z.get_max_twisted_atom(sh)
        
        assert (g.inverse() * a).length() == a.length() - g.length()
        a = a.inverse() * g
        btoms = set(span(a))
        assert atoms == btoms


def _test_dtype_fpf_atoms_by_shape(z, atoms_by_shape):
    n = z.rank

    def span(v):
        v = v.oneline
        level = {(None, v)}
        while level:
            nextlevel = set()
            for u, v in level:
                yield SignedPermutation(*v).inverse()

                for i in range(0 if n % 2 == 0 else 1, len(v) - 3, 2):
                    b, c, a, d = v[i: i + 4]
                    if a < b < c < d:
                        w = v[:i] + (a, d, b, c) + v[i + 4:]
                        nextlevel.add((v, w))

            level = nextlevel

    z = EvenSignedPermutation(*z)
    for sh, atoms in atoms_by_shape.items():
        a = z.get_max_fpf_atom(sh)
        btoms = set(span(a.inverse()))
        assert atoms == btoms


def _test_refinement(clan, atoms_by_shape, expected_shapes):
    clan_atoms = set(clan.get_atoms())
    
    actual_shapes = set()
    otherset = set()
    for sh, subset in atoms_by_shape.items():
        if subset.issubset(clan_atoms):
            otherset |= subset
            actual_shapes.add(sh)
        else:
            assert len(subset & clan_atoms) == 0
    assert otherset == clan_atoms
    if expected_shapes != actual_shapes:
        print(clan)
        print('actual shapes')
        for a in actual_shapes:
            print('  ', a)
        print('expected shapes')
        for a in expected_shapes:
            print('  ', a)
        print('extra shapes')
        for a in expected_shapes - actual_shapes:
            print('  ', a)
        assert actual_shapes.issubset(expected_shapes)
    assert expected_shapes == actual_shapes
