from schubert import (
    Schubert,
    Grothendieck,
    InvSchubert,
    FPFSchubert,
    FPFGrothendieck,
    InvGrothendieck,
)
from polynomials import (
    MPolynomial,
    Operator
)
from permutations import Permutation
import schubert
import pytest


def test_inv_vex_tab_formula(n=5):
    for w in Permutation.involutions(n):
        if w.is_vexillary():
            a = InvSchubert.vexillary_tableau_formula(w)
            b = InvSchubert.get(w) * 2**w.number_two_cycles()
            assert a == b


def test_formal_operators():
    D = Operator.create(1) # noqa
    x1 = MPolynomial.monomial(1)
    x2 = MPolynomial.monomial(2)
    x3 = MPolynomial.monomial(3)

    (x1 * x2).toggle(2) == x1 * x3
    (x1 * x2).toggle(1) == x1 * x2
    x1.toggle(1) == x2
    x2.toggle(1) == x1

    E = D * x1 # noqa
    F = x1 * D # noqa

    print('1:', D, type(D))
    print('2:', E, type(E))
    print('3:', F, type(F))
    print()

    D + D
    D - D
    D * D
    assert D * x1 != x1 * D
    assert D * x2 != x2 * D
    assert x3 * D == D * x3

    assert D * x1 == Operator.create() + x2 * D

    A = Operator.create(1) * x1 # noqa
    B = Operator.create(2) * x2 # noqa

    print(A)
    print(Operator.create(1) * x2 * Operator.create(1))
    print(A * Operator.create(1))
    print(A * Operator.create(1) - Operator.create(1))
    print()

    assert A * Operator.create() == A
    assert Operator.create() * A == A
    assert A * Operator.create(1) == Operator.create(1)
    assert Operator.create(1) * A == 0
    assert A * A == A
    assert B * B == B
    assert A * B * A == B * A * B

    a = Operator.create(1) * (1 - x2)
    b = Operator.create(2) * (1 - x3)
    assert a * a == a
    assert b * b == b
    assert a * b * a == b * a * b

    A = Operator.create(1) * x1 * (1 - x2) # noqa
    B = Operator.create(2) * x2 * (1 - x3) # noqa
    print(A * B * A)
    print(B * A * B)
    print()

    assert A * A == A
    assert B * B == B
    assert A * B * A == B * A * B


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


def test_inv_grothendieck(n=4):
    i = list(Permutation.involutions(n))
    s = {w: InvGrothendieck.get(w) for w in i}
    for w in s:
        if not w.is_vexillary():
            assert s[w] == 0
        else:
            assert s[w] != 0


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


@pytest.mark.slow
def test_fpf_transitions():
    beta = FPFGrothendieck.beta

    for n in [2, 4, 6, 8]:
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

        g = list(Permutation.fpf_involutions(n))
        for w in g:
            cyc = [
                (i, j) for i, j in w.cycles
            ]
            w = w * Permutation.s_i(n + 1)
            for i, j in cyc:
                var = (1 + beta * schubert.x(i)) * (1 + beta * schubert.x(j))
                ts = []
                for k in range(1, i):
                    t = Permutation.cycle([k, i])
                    v = t * w * t
                    if v.fpf_involution_length() == w.fpf_involution_length() + 1:
                        ts.append(k)
                ttt = [(w, 1)]
                for k in ts:
                    t = Permutation.cycle([k, i])
                    ttt += [(t * v * t, beta * a) for v, a in ttt]
                f = 0
                for v, a in ttt:
                    f += FPFGrothendieck.get(v) * a
                f = f * var
                sp = ''.join(['(1 + beta t_{%s,%s})' % (k, i) for k in ts]) if ts else '1'
                print('G_%s * %s * (%s) = ' % (w, sp, var))
                print()

                # try:
                #     print('    ', f)
                #     dec = FPFGrothendieck.decompose(f)
                # except:
                #     print('     halted computation')
                #     assert False
                # print()
                # print('    ', dec)
                # print()

                a = FPFGrothendieck.get(w)
                for z in terms(w, j):
                    len_diff = z.fpf_involution_length() - w.fpf_involution_length()
                    a += FPFGrothendieck.get(z) * beta**len_diff
                assert f == a

                print()
                print()


@pytest.mark.slow
def test_lenart_grothendieck_transitions():
    for n in [1, 2, 3, 4, 5, 6]:
        def terms(w, j):
            queue = [(w, j - 1, 0)]
            while queue:
                y, i, q = queue[0]
                queue = queue[1:]
                if i <= 0:
                    queue.append((y, n + 1, q))
                    continue
                if i == j:
                    continue
                s = Permutation.transposition(i, j)
                z = y * s
                if z.length() == y.length() + 1:
                    b = q + 1 if i < j else q
                    yield z, (-1)**b * Grothendieck.beta**(z.length() - w.length() - 1)
                    queue.append((z, i - 1, b))
                queue.append((y, i - 1, q))

        g = list(Permutation.all(n))
        for w in g:
            for i in range(1, n + 1):
                f = Grothendieck.get(w) * schubert.x(i)

                print('w =', w, 'i =', i)
                print([(z, c, Grothendieck.get(z)) for z, c in terms(w, i)])
                print(f)
                g = 0
                for z, c in terms(w, i):
                    g += Grothendieck.get(z) * c
                print(g)
                print()
                assert f == g


@pytest.mark.slow
def test_grothendieck_transitions():
    for n in [1, 2, 3, 4, 5, 6]:
        def terms(w, j):
            queue = [(w, n + 1)]
            while queue:
                y, k = queue[0]
                queue = queue[1:]

                if k <= j:
                    continue

                s = Permutation.transposition(j, k)
                z = y * s
                if z.length() == y.length() + 1:
                    yield z
                    queue.append((z, k - 1))
                queue.append((y, k - 1))

        g = list(Permutation.all(n))
        for w in g:
            for i in range(1, n + 1):
                var = 1 + Grothendieck.beta * schubert.x(i)

                ts = []
                for k in range(1, i):
                    t = Permutation.cycle([k, i])
                    v = w * t
                    if v.length() == w.length() + 1:
                        ts.append(k)

                ttt = [(w, 1)]
                for k in ts:
                    t = Permutation.cycle([k, i])
                    ttt += [(v * t, Grothendieck.beta * a) for v, a in ttt]

                f = 0
                for v, a in ttt:
                    f += Grothendieck.get(v) * a
                f = f * var

                sp = ''.join(['(1 + beta t_{%s,%s})' % (k, i) for k in ts]) if ts else '1'
                print('G_%s * %s * (%s) = ' % (w, sp, var))
                print()

                try:
                    dec = Grothendieck.decompose(f)
                except:
                    print('     halted computation')
                    assert False

                print()
                print('    ', dec)
                print()

                a = Grothendieck.get(w)
                for z in terms(w, i):
                    len_diff = z.length() - w.length()
                    a += Grothendieck.beta**len_diff * Grothendieck.get(z)
                assert f == a

                print()
                print()


@pytest.mark.slow
def test_product_grothendieck():
    n = 6
    g = list(Permutation.all(n))

    for w in g:
        a = Grothendieck.get(w) == Grothendieck.product(w)
        b = w.is_dominant()
        c = Grothendieck.get(w) == Schubert.get(w)
        assert a == b
        assert b == c


@pytest.mark.slow
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
