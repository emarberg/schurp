from permutations import Permutation



def test_inversions_commutations(n=5):
    for w in Permutation.involutions(n):
        for word in w.get_involution_words():
            inv = Permutation.inversions(word)
            com = Permutation.find_commutations(word)
            print(word, inv, com)
            for i in range(len(inv)):
                a, b = inv[i]
                if i in com:
                    assert b == w(a)
                else:
                    assert b != w(a)


def test_inversions(n=5):
    for w in Permutation.all(n):
        word = w.get_reduced_word()
        inv = Permutation.inversions(word)
        exp = {(a, b) for a in range(1, n) for b in range(a + 1, n + 1) if w(a) > w(b)}
        # print(w, word, inv, exp)
        # print()
        assert exp == set(inv)


def test_init():
    w = Permutation(1)
    v = Permutation([1])
    u = Permutation(*[1])
    assert u == v == w


def test_lt():
    assert Permutation(1, 3, 2) < Permutation(2, 1)
    assert Permutation() < Permutation(2, 1)


def test_fpf_atoms():
    z = Permutation.cycles([[1, 8], [2, 6], [3, 10], [4, 7], [5, 9]])
    assert tuple(z.get_min_fpf_atom().inverse().oneline) == (1, 8, 2, 6, 3, 10, 4, 7, 5, 9)
    assert len(z.get_fpf_atoms()) == 8
    assert Permutation(2, 6, 4, 7, 1, 8, 5, 9, 3, 10).inverse() in z.get_fpf_atoms()


def trim(t):
    a = list(t)
    while a and a[-1] == 0:
        a = a[:-1]
    return tuple(a)


def test_codes():
    for i in range(4):
        for j in range(4):
            for k in range(4):
                a = [0] + [i, j, k]
                b = [i] + [0] + [j, k]
                c = [i, j] + [0] + [k]
                d = [i, j, k] + [0]

                for x in [a, b, c, d]:
                    w = Permutation.from_code(x)
                    assert w.code() == trim(x)

    assert Permutation([1, 2, 3, 4, 5]).code() == ()


def test_involution_codes():
    for w in Permutation.involutions(5):
        assert sum(w.involution_code()) == w.involution_length()

        v = Permutation.from_involution_code(w.involution_code())
        print('w =', w)
        print('v =', v)
        print(w.involution_code())
        print()

        assert v == w
        assert w.involution_code() == w.get_min_atom().code()


def test_fpf_involution_codes():
    for w in Permutation.fpf_involutions(6):
        assert all(w.fpf_involution_length() == a.length() for a in w.get_fpf_atoms())
        assert sum(w.fpf_involution_code()) == w.fpf_involution_length()

        v = Permutation.from_fpf_involution_code(w.fpf_involution_code())
        print('w =', w, type(w))
        print('v =', v, type(v))
        print(w.fpf_involution_code())
        print()

        assert v == w
        assert trim(w.fpf_involution_code()) == w.get_min_fpf_atom().code()


def test_from_fpf_involution_codes():
    w = Permutation(3, 5, 1, 6, 2, 4)
    assert Permutation.from_fpf_involution_code(w.fpf_involution_code()) == w


def test_is_perfect(m=5):
    def is_perfect(w, n):
        for t in Permutation.reflections(n):
            wtw = w * t.star(n) * w.star(n)
            if t * wtw != wtw * t:
                return False
        return True

    for n in range(m + 1):
        seen = {}
        for w in Permutation.twisted_involutions(n):
            if is_perfect(w, n):
                seen[w] = 1

        for w in seen:
            print(n, ':', (w * Permutation.longest_element(n)).cycle_repr())
        # assert len(seen) == 1 and list(seen)[0].is_identity()

