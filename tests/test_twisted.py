from permutations import Permutation


def test_twisted_involutions(n=5):
    w0 = Permutation.longest_element(n)
    assert set(Permutation.twisted_involutions(n)) == {
        w0 * w for w in Permutation.involutions(n)
    }


def test_twisted(n=5):
    for w in Permutation.twisted_involutions(n):
        for a in w.get_twisted_atoms(n):
            v = a.star(n).inverse() % a
            print(w, a, a.get_reduced_word(), v)
            assert v == w
            assert a.length() == w.twisted_involution_length(n)


def _test_shape(w, n):
    y = w.star(n).inverse() % w
    assert y.twisted_involution_length(n) == w.length()

    w0 = Permutation.longest_element(n)
    y = w0 * y
    aword = reversed(w.get_reduced_word())

    yfixed = {i for i in range(1, n + 1) if y(i) == i}
    v = Permutation()
    sh = set()
    for a in aword:
        if a > 0 and y(a) == a and y(a + 1) == a + 1:
            e, f = tuple(sorted([v(a), v(a + 1)]))
            sh |= {(e, f)}
        s = Permutation.s_i(a)
        v *= s
        z = s % y % s
        assert z.involution_length() == y.involution_length() + 1
        y = z
    f = {i for p in sh for i in p}
    return sh | {(i, i) for i in yfixed - f}


def test_twisted_shape(n=5):
    w0 = Permutation.longest_element(n)
    for w in Permutation.twisted_involutions(n):
        shapes = {}
        for a in w.get_twisted_atoms(n):
            sh = tuple(sorted(a.twisted_shape(n)))
            shapes[sh] = shapes.get(sh, []) + [a]
        print(w, ' :: ', (w0 * w).fixed(n))
        print()
        for sh, atoms in shapes.items():
            print(' ', set(sh), '->', atoms)
            print(' ', len(str(set(sh))) * ' ', '  ', [a.inverse() for a in atoms])
            print()
        print()
        print()
        for sh, atoms in shapes.items():
            sh = set(sh)
            for a in atoms:
                if _test_shape(a, n) != sh:
                    print(' *', a, ':', sh, '!=', _test_shape(a, n))
                    raise Exception
