from even import EvenSignedPermutation


def test_get_minimal_fpf_involution():
    cls = EvenSignedPermutation
    assert cls.get_minimal_fpf_involution(1) == cls(1)
    assert cls.get_minimal_fpf_involution(2) == cls(2, 1)
    assert cls.get_minimal_fpf_involution(3) == cls(1, 3, 2)
    assert cls.get_minimal_fpf_involution(4) == cls(2, 1, 4, 3)


def test_get_involution_word():
    s = EvenSignedPermutation.s_i(0, 3)
    t = EvenSignedPermutation.s_i(2, 3)
    assert (t * s * t).get_involution_word() in {(0, 2), (2, 0)}
    assert (s * t * s).get_involution_word() in {(0, 2), (2, 0)}
    assert set((t * s * t).get_involution_words()) == {(0, 2), (2, 0)}
    assert set((s * t * s).get_involution_words()) == {(0, 2), (2, 0)}


def test_fpf_involution_words():
    for n in [2, 3]:
        for w in EvenSignedPermutation.fpf_involutions(n):
            words = set(w.get_fpf_involution_words())
            if words:
                print(w, '->', len(words))


def test_get_atoms():
    for n in [3, 4]:
        i = list(EvenSignedPermutation.involutions(n))
        for w in i:
            words = set()
            for a in w.get_atoms():
                assert len(a) == w.involution_length()
                assert a.inverse() % a == w
                words |= set(a.get_reduced_words())
            assert words == set(w.get_involution_words())


def test_get_twisted_atoms():
    for n in [3, 4]:
        i = list(EvenSignedPermutation.twisted_involutions(n))
        for w in i:
            words = set()
            for a in w.get_twisted_atoms():
                assert len(a) == w.twisted_involution_length()
                assert a.inverse().star() % a == w
                words |= set(a.get_reduced_words())
            assert words == set(w.get_twisted_involution_words())


def test_twisted_involutions(n=5):
    assert n % 2 != 0
    w0 = EvenSignedPermutation.longest_element(n)
    assert set(EvenSignedPermutation.twisted_involutions(n)) == {
        w0 * w for w in EvenSignedPermutation.involutions(n)
    }


def test_twisted(n=4):
    for w in EvenSignedPermutation.twisted_involutions(n):
        for a in w.get_twisted_atoms():
            v = a.star().inverse() % a
            print(w, a, a.get_reduced_word(), v)
            assert v == w
            assert a.length() == w.twisted_involution_length()


def _nneg(w):
    des = 0
    for i in range(len(w) - 1):
        a, b = w[i:i + 2]
        if a > b:
            des += 1
            v = w[:i] + w[i + 2:]
            for nn in _nneg(v):
                yield {(abs(a), b)} | nn
    if des == 0:
        yield {w}


def test_nneg(n=4):
    for w in EvenSignedPermutation.involutions(n):
        for atom in w.get_atoms():
            o = atom.inverse().oneline
            nn = {tuple(sorted(x)) for x in _nneg(o)}

            ndes, fix = atom.ndes(), atom.nfix()
            neg = []
            for a, b in ndes:
                if a < -b:
                    neg += [a, -b]
            neg = tuple(sorted(neg))

            sh = {(a, -b) for a, b in ndes if 0 < a < -b}
            sh |= {(b, -a) for a, b in ndes if 0 < a < -b}
            sh = tuple(sorted(sh))

            pair = sorted((b, a) for a, b in ndes if not (0 < a < -b))

            cfix = tuple(i for i in w.fixed_points() if i > 0)
            cneg = tuple(i for i in w.negated_points() if i > 0)
            csh = tuple(sorted(atom.shape()))
            cpair = w.pair()

            print(w, ':', list(o))
            print('         shape:', csh, '==', sh)
            print('         pairs:', cpair, '==', pair)
            print('  fixed points:', cfix, '==', fix)
            print('negated points:', cneg, '==', neg)
            print()
            assert fix == cfix
            assert neg == cneg
            assert len(nn) == 1
            assert csh == sh
            assert cpair == pair


def _nneg(w):
    des = 0
    for i in range(len(w) - 1):
        a, b = w[i:i + 2]
        if a > b:
            des += 1
            v = w[:i] + w[i + 2:]
            for nn in _nneg(v):
                yield {(abs(a), b)} | nn
    if des == 0:
        yield {w}


def test_twisted_nneg(n=4):
    for w in EvenSignedPermutation.twisted_involutions(n):
        for atom in w.get_twisted_atoms():
            o = atom.inverse().oneline[1:]
            nn = {tuple(sorted(x)) for x in _nneg(o)}

            ndes, fix = atom.twisted_ndes(), atom.twisted_nfix()
            neg = []
            for a, b in ndes:
                if a < -b:
                    neg += [a, -b]
                elif a == -b:
                    neg += [a]
            neg = tuple(sorted(neg))

            sh = {(-a, a) for a, b in ndes if a + b == 0}
            sh |= {(a, -b) for a, b in ndes if 0 < a < -b}
            sh |= {(b, -a) for a, b in ndes if 0 < a < -b}
            sh = tuple(sorted(sh))

            pair = sorted((b, a) for a, b in ndes if not (0 < a <= -b))

            cfix = tuple(i for i in w.twisted_fixed_points() if i > 0)
            cneg = tuple(i for i in w.twisted_negated_points() if i > 0)
            csh = tuple(sorted(atom.twisted_shape()))
            cpair = w.pair()

            print(w, ':', list(o))
            print('         shape:', csh, '==', sh)
            print('         pairs:', cpair, '==', pair)
            print('  fixed points:', cfix, '==', fix)
            print('negated points:', cneg, '==', neg)
            print()
            assert fix == cfix
            assert neg == cneg
            assert len(nn) == 1
            assert csh == sh
            assert cpair == pair


def test_shape(n=4):
    cls = EvenSignedPermutation
    w0 = cls.longest_element(n)
    for w in cls.involutions(n):
        shapes = {}
        for a in w.get_atoms():
            sh = tuple(sorted(a.shape()))
            shapes[sh] = shapes.get(sh, []) + [a]
        print(w, ' :: ', w * w0, ' :: ', (w * w0).involution_fixed_points(n % 2 != 0))
        print()
        for sh, atoms in shapes.items():
            print(' ', set(sh), '->', atoms)
            print(' ', len(str(set(sh))) * ' ', '  ', [a.inverse() for a in atoms])
            print()
            assert not any(i == j for i, j in sh)
        print()
        print()


def test_twisted_shape(n=4):
    cls = EvenSignedPermutation
    w0 = cls.longest_element(n)
    for w in cls.twisted_involutions(n):
        shapes = {}
        for a in w.get_twisted_atoms():
            sh = tuple(sorted(a.twisted_shape()))
            shapes[sh] = shapes.get(sh, []) + [a]
        print(w, ' :: ', w0 * w, ' :: ', (w0 * w).involution_fixed_points(n % 2 == 0))
        print()
        for sh, atoms in shapes.items():
            print(' ', set(sh), '->', atoms)
            print(' ', len(str(set(sh))) * ' ', '  ', [a.inverse() for a in atoms])
            print()
            assert len([i for i, j in sh if i == -j]) == 1
        print()
        print()
