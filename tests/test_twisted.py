from permutations import Permutation


def test_twisted(n=5):
    for w in Permutation.twisted_involutions(n):
        for a in w.get_twisted_atoms(n):
            v = a.star(n).inverse() % a
            print(w, a, a.get_reduced_word(), v)
            assert v == w
            assert a.length() == w.twisted_involution_length(n)


def test_twisted_shape(n=7):
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
    assert False