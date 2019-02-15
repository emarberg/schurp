from schubert import (
    Schubert,
    Grothendieck,
    InvSchubert,
    FPFSchubert,
    FPFGrothendieck,
    MPolynomial
)
import schubert
from permutations import Permutation


def test_divided_differences():
    x = MPolynomial.monomial(3)
    y = MPolynomial.monomial(4)
    f = x**3
    assert f.divided_difference(3) == x**2 + x * y + y**2
    assert (x * y).isobaric_divided_difference(3) == x * y


def test_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation.s_i(3)
    assert Schubert.get(w) == x + y + z

    w = Permutation(3, 4, 1, 2)
    assert Schubert.get(w) == x * x * y * y

    w = Permutation(1, 4, 3, 2)
    assert Schubert.get(w) == x**2 * y + x**2 * z + x * y**2 + x * y * z + y**2 * z


def test_inv_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation()
    assert InvSchubert.get(w) == 1

    w = Permutation(1, 4, 3, 2)
    assert InvSchubert.get(w) == x**2 + 2 * x * y + x * z + y**2 + y * z

    w = Permutation(3, 4, 1, 2)
    assert InvSchubert.get(w) == x**2 * y + x * y**2


def test_fpf_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation(2, 1, 4, 3)
    assert FPFSchubert.get(w) == 1

    w = Permutation(3, 4, 1, 2)
    assert FPFSchubert.get(w) == x + y

    w = Permutation(4, 3, 2, 1)
    assert FPFSchubert.get(w) == x * x + x * y + x * z + y * z


def test_fpf_grothendieck():
    n = 6
    g = list(Permutation.fpf_involutions(n))

    for w in g:
        a = FPFGrothendieck.get(w) == FPFGrothendieck.product(w)
        b = w.is_fpf_dominant()
        assert a == b
        f = FPFGrothendieck.get(w)
        g = sum([
            (-1) ** (u.length() - w.fpf_involution_length()) * Grothendieck.get(u)
            for u in w.get_symplectic_hecke_atoms()
        ])
        print(w)
        print(f)
        print(g)
        print()
        print(list(w.get_symplectic_hecke_atoms()))
        print()
        assert f == g


def test_min():
    w = Permutation(3, 5, 1, 6, 2, 4)
    f = FPFSchubert.get(w)
    assert FPFSchubert.least_term(f) == {2: 1, 4: 1}


def test_fpf_transitions():
    n = 4

    def terms(w, j):
        queue = [(w, n + 1)]
        while queue:
            y, k = queue[0]
            queue = queue[1:]

            if k <= j:
                continue

            s = Permutation.transposition(j, k)
            z = s * y * s
            if z.fpf_involution_length() == y.fpf_involution_length() + 1:
                yield z
                queue.append((z, k - 1))
            queue.append((y, k - 1))

    # g = list(Permutation.fpf_involutions(n))
    # for w in g:
    #     cyc = [
    #         (i, j) for i, j in w.cycles
    #         if not any(k < i and l < j for k, l in w. cycles)
    #     ]
    #     w = w * Permutation.s_i(n + 1)
    #     for i, j in cyc:
    #         v = schubert.x(i) + schubert.x(j) - schubert.x(i) * schubert.x(j)
    #         f = FPFGrothendieck.get(w) * v
    #         a = 0
    #         print(list(terms(w, j)))
    #         for z in terms(w, j):
    #             if (z.fpf_involution_length() - w.fpf_involution_length()) % 2 == 0:
    #                 sgn = -1
    #             else:
    #                 sgn = 1
    #             a += FPFGrothendieck.get(z) * sgn
    #         print('G_%s * (%s)' % (w, v))
    #         print(f)
    #         print()
    #         print(a)
    #         print()
    #         print()
    #         assert f == a

    g = list(Permutation.fpf_involutions(6))
    for v in g[1:]:
        c = max(v.get_fpf_visible_inversions())
        t = Permutation.cycle(c)
        u = t * v * t
        s = set(v.fpf_rothe_diagram()) - set(u.fpf_rothe_diagram())

        print(v, c)
        v.print_fpf_rothe_diagram()
        u.print_fpf_rothe_diagram()
        print()
        print(s)
        print()

        assert len(s) == 1
        i, j = s.pop()

        f = FPFGrothendieck.get(v)
        g = FPFGrothendieck.get(u)
        d = g - f

        # print(f)
        # print()
        # print(g)
        # print()
        # print(d)
        # print()

        e = d.divide_linear(i, -1).divide_linear(j, -1)
        #print(e)
        #print()
        print('G_u - G_v = (1 - x_%i) (1 - x_%i) [' % (i, j), FPFGrothendieck.decompose(e), ']')
        print()
        print()
        print()
        print()
    assert False
