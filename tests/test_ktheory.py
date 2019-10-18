from schubert import X, Vector, Permutation


"""

X(3) * (1,3,4,5,2) == \
(1,3,5,4,2) - (1,4,3,5,2) + \
(1,4,5,3,2) - (3,4,1,5,2) + \
(3,4,5,1,2) + (3,4,2,5,1) - \
(3,4,5,2,1)

\[
(1 + \beta x_3) \fkG_{13452} =
\fkG_{13452} +
\beta\fkG_{13542} - \beta\fkG_{14352} -
\beta^2\fkG_{14532} + \beta^2\fkG_{34152} +
\beta^3\fkG_{34512} + \beta^3\fkG_{34251} +
\beta^4\fkG_{34521}
\]


f = (1 + X(0) * X(3)) * Grothendieck.get(Permutation(1,3,4,5,2))
g = Grothendieck.get(Permutation(1,3,4,5,2)) + \
X(0) * Grothendieck.get(Permutation(1,3,5,4,2)) - \
X(0) * Grothendieck.get(Permutation(1,4,3,5,2)) - \
X(0)**2 * Grothendieck.get(Permutation(1,4,5,3,2)) + \
X(0)**2 * Grothendieck.get(Permutation(3,4,1,5,2)) + \
X(0)**3 * Grothendieck.get(Permutation(3,4,5,1,2)) + \
X(0)**3 * Grothendieck.get(Permutation(3,4,2,5,1)) + \
X(0)**4 * Grothendieck.get(Permutation(3,4,5,2,1))
"""


def grothendieck_transition_terms(w, j):
    """
    Yields pairs (z, sign * beta**n) where

        z = w(a_1,k)(a_2,k)···(a_p,k)(k,b_1)(k,b_2)···(k,b_q)
        n = len(z) - len(w) - 1
        sign = (-1)**p

    and a_p < ... < a_2<a_1 < k < b_q < ... < b_2 < b_1, and the length
    increases by exactly one upon multiplication by each transposition.
    """
    beta = X(0)
    n = max(j, len(w.oneline))
    yield w, beta**0
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
            yield z, (-1)**b * beta**(z.length() - w.length())
            queue.append((z, i - 1, b))
        queue.append((y, i - 1, q))


def test_grothendieck_transition_terms():
    w = Permutation(1, 3, 4, 5, 2)
    j = 3
    assert set(grothendieck_transition_terms(w, j)) == {
        (Permutation(1, 3, 4, 5, 2), X(0)**0),
        (Permutation(1, 3, 5, 4, 2), X(0)),
        (Permutation(1, 4, 3, 5, 2), -X(0)),
        (Permutation(1, 4, 5, 3, 2), -X(0)**2),
        (Permutation(3, 4, 1, 5, 2), X(0)**2),
        (Permutation(3, 4, 5, 1, 2), X(0)**3),
        (Permutation(3, 4, 2, 5, 1), X(0)**3),
        (Permutation(3, 4, 5, 2, 1), X(0)**4)
    }


def fpf_transition_upper_terms(w, j):
    assert w(j) < j

    def trim(z):
        m = len(z.oneline)
        if m % 2 == 0 and z(m) == m - 1:
            return trim(z * Permutation.s_i(m - 1))
        return z

    beta = X(0)
    n = len(w.oneline)
    yield trim(w), beta**0

    w *= Permutation.s_i(n + 1)
    queue = [(w, n + 1)]
    while queue:
        y, k = queue[0]
        queue = queue[1:]

        if k <= j:
            continue

        s = Permutation.transposition(j, k)
        z = s * y * s
        if z.fpf_involution_length() == y.fpf_involution_length() + 1:
            yield trim(z), beta ** (z.fpf_involution_length() - w.fpf_involution_length())
            queue.append((z, k - 1))
        queue.append((y, k - 1))


def fpf_transition_lower_terms(w, j):
    assert j < w(j)

    beta = X(0)
    yield w, beta**0

    queue = [(w, 1)]
    while queue:
        y, k = queue[0]
        queue = queue[1:]

        if k >= j:
            continue

        s = Permutation.transposition(j, k)
        z = s * y * s
        if z.fpf_involution_length() == y.fpf_involution_length() + 1:
            yield z, beta ** (z.fpf_involution_length() - w.fpf_involution_length())
            queue.append((z, k + 1))
        queue.append((y, k + 1))


def multiply_via_grothendieck_transitions(w, j):
    if type(w) == Permutation:
        w = Vector({w: X(0)**0})
    assert type(w) == Vector
    ans = Vector()
    for y, b in w.items():
        for z, c in grothendieck_transition_terms(y, j):
            ans += Vector({z: b * c})
    return ans


def test_fpf_transitions(n=6):
    beta = X(0)
    for z in Permutation.fpf_involutions(n):
        for j, k in z.cycles:
            print('n =', n)
            print('z =', z)
            print('j =', j)
            print('k =', k)
            print()

            f = Vector()
            for y, c in fpf_transition_lower_terms(z, j):
                for w in y.get_symplectic_hecke_atoms():
                    f += Vector({w: c * beta ** (len(w) - y.fpf_involution_length())})
            f = multiply_via_grothendieck_transitions(f, j)
            f = multiply_via_grothendieck_transitions(f, k)
            print('LHS =', f)

            g = Vector()
            for y, c in fpf_transition_upper_terms(z, k):
                for w in y.get_symplectic_hecke_atoms():
                    g += Vector({w: c * beta ** (len(w) - y.fpf_involution_length())})
            print('RHS =', g)
            print()

            print()
            assert f == g

