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

beta = X(0)


def grothendieck_transition_terms(w, j):
    """
    Yields pairs (z, sign * beta**n) where

        z = w(a_1,k)(a_2,k)···(a_p,k)(k,b_1)(k,b_2)···(k,b_q)
        n = len(z) - len(w) - 1
        sign = (-1)**p

    and a_p < ... < a_2<a_1 < k < b_q < ... < b_2 < b_1, and the length
    increases by exactly one upon multiplication by each transposition.
    """
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


def inv_transition_upper_terms(w, j, n=None, forbidden=()):
    if n is None:
        n = max(j, len(w.oneline)) + 1
        yield w, beta**0

    queue = [(w, n)]
    while queue:
        y, k = queue[0]
        queue = queue[1:]

        if k <= j:
            continue

        z = y.tau_ij(j, k)
        repeats = (y(j) == j < y(k) < k and z == y.tau_ij(j, y(k)))
        
        if k not in forbidden and not repeats and z.involution_length() == y.involution_length() + 1:    
            yield z, beta ** (z.involution_length() - w.involution_length())
            queue.append((z, k - 1))

            if z.number_two_cycles() < y.number_two_cycles():
                extra = tuple(a for a in range(y(k) + 1, k) if y(a) > k)
                if extra:
                    print('RHS *', repeats, y, j, k, z, 'forbid =', forbidden + extra)
                for u, _ in inv_transition_upper_terms(z, y(k), k, forbidden + extra):
                    if extra:
                        print('RHS new term:', u)
                    yield u, beta ** (u.involution_length() - w.involution_length())
                    queue.append((u, k - 1))
            
        queue.append((y, k - 1))


def inv_transition_lower_terms(w, j, n=None, forbidden=()):
    if n is None:
        yield w, beta**0
        n = 1
    
    queue = [(w, n)]
    while queue:
        y, k = queue[0]
        queue = queue[1:]

        if k >= j:
            continue

        z = y.tau_ij(k, j)
        repeats = (y(j) == j > y(k) > k and z == y.tau_ij(y(k), j))

        if k not in forbidden and not repeats and z.involution_length() == y.involution_length() + 1:
            yield z, beta ** (z.involution_length() - w.involution_length())
            queue.append((z, k + 1))

            if z.number_two_cycles() < y.number_two_cycles():
                extra = tuple(a for a in range(k + 1, y(k)) if y(a) < k)
                if extra:
                    print('LHS *', y, k, j, z, 'forbid =', forbidden + extra)
                for u, _ in inv_transition_lower_terms(z, y(k), k, forbidden + extra):
                    if extra:
                        print('LHS new term:', u)
                    yield u, beta ** (u.involution_length() - w.involution_length())
                    queue.append((u, k + 1))

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


def fpfgroth(y):
    f = Vector()
    for w in y.get_symplectic_hecke_atoms():
        f += Vector({w: beta ** (len(w) - y.fpf_involution_length())})
    return f


def invgroth(y):
    f = Vector()
    for w in y.get_involution_hecke_atoms():
        f += Vector({w: beta ** (len(w) - y.involution_length())})
    return f


def test_fpf_transitions(n=4):
    for z in Permutation.fpf_involutions(n):
        for j, k in z.cycles:
            print('n =', n)
            print('z =', z)
            print('j =', j)
            print('k =', k)
            print()

            f = Vector()
            for y, c in fpf_transition_lower_terms(z, j):
                f += fpfgroth(y) * c
            f = multiply_via_grothendieck_transitions(f, j)
            f = multiply_via_grothendieck_transitions(f, k)
            print('LHS =', f)

            g = Vector()
            for y, c in fpf_transition_upper_terms(z, k):
                g += fpfgroth(y) * c
            print('RHS =', g)
            print()

            print()
            assert f == g


def decomposeinv(vec):
    ans = Vector()
    if vec:
        wl = [(w.involution_length(), w) for w in vec]
        w = sorted(wl)[0][1]
        z = w.inverse() % w
        ans = Vector({z: vec[w]})
        vec -= invgroth(z) * vec[w]
        ans += decomposeinv(vec)
    return ans


def test_inv_transitions(n=4):
    for z in Permutation.involutions(n): #[Permutation(2,1,6,8,7,3,5,4)]: #
        cyc = [(j, k) for j, k in z.get_two_cycles()] + [(j, j) for j in z.fixed(n + 1)]
        for j, k in cyc:
            print('n =', n, 'z =', z, 'j =', j, 'k =', k)

            f = Vector()
            for y, c in inv_transition_lower_terms(z, j):
                f += invgroth(y) * c
            f = multiply_via_grothendieck_transitions(f, j)
            if j < k:
                f = multiply_via_grothendieck_transitions(f, k)
            try:
                decomposeinv(f)
            except:
                raise Exception

            g = Vector()
            for y, c in inv_transition_upper_terms(z, k):
                g += invgroth(y) * c
            
            if f != g:
                print()
                # print('LHS =', f)
                # print()
                # print('RHS =', g)
                # print()
                try:
                    dec = decomposeinv(f - g)
                    for w, coeff in sorted(dec.items(), key=lambda a: a[0].involution_length()):
                        print('  ', w, coeff)
                except:
                    print('diff =', f - g)
                    return f - g
                print()
                input('?')
            print()

            #assert f == g


def lessdot(z, i, j):
    return z(i) < z(j) and not any(z(i) < z(e) < z(j) for e in range(i + 1, j))


def kappa(z):
    return z.number_two_cycles()


def gamma(z, i, j):
    if lessdot(z, i, j):
        t = z.transposition(i, j)
        if z(i) == i and z(j) == j:
            return z * t

        if z(i) <= i < j <= z(j) or i < j < z(i) < z(j) or z(i) < z(j) < i < j:
            return t * z * t

        if i < z(i) < j < z(j):
            return t * z * t * z.transposition(z(i), j)

        if z(i) < i < z(j) < j:
            return t * z * t * z.transposition(i, z(j))
    return z
        

def uop(v, i, j):
    t = v.transposition(i, j)
    if v * t == t * v:
        return v * t
    else:
        return t * v * t


def aop(v, j, S):
    ilist = sorted([i for i in S if lessdot(v, i, j)])
    queue = [(v, 1, 0)]
    while queue:
        v, c, t = queue[0]
        queue = queue[1:]

        if t == len(ilist):
            yield (v, c)
            continue

        queue.append((v, c, t + 1))
        i = ilist[t]
        g = uop(v, i, j)
        queue.append((g, c * beta, t + 1))


def bop(v, j, S):
    klist = sorted([k for k in S if lessdot(v, j, k)], reverse=True)
    queue = [(v, 1, 0)]
    while queue:
        v, c, t = queue[0]
        queue = queue[1:]

        if t == len(klist):
            yield (v, c)
            continue

        queue.append((v, c, t + 1))
        k = klist[t]
        g = uop(v, j, k)
        queue.append((g, c * beta, t + 1))


def A(z, a, b):
    queue = [(z, 1, a, b)]
    while queue:
        v, c, i, j = queue[0]
        queue = queue[1:]

        if i == j:
            yield (v, c)
            continue

        g = gamma(v, i, j)
        R = [r for r in range(i + 1, v(i)) if v(r) > i]

        if i < j and v == g:
            queue.append((v, c, i + 1, j))
        elif i < j and kappa(v) > kappa(g):
            queue.append((v, c, i + 1, j))
            for (w, d) in aop(g, v(i), R):
                queue.append((w, c * d * beta, i + 1, j))
        else:
            queue.append((v, c, i + 1, j))
            queue.append((g, c * beta, i + 1, j))


def B(z, a, b):
    queue = [(z, 1, a, b)]
    while queue:
        v, c, j, k = queue[0]
        queue = queue[1:]

        if j == k:
            yield (v, c)
            continue

        g = gamma(v, j, k)
        S = [s for s in range(v(k) + 1, k) if v(s) < k]

        if j < k and v == g:
            queue.append((v, c, j, k - 1))
        elif j < k and kappa(v) > kappa(g):
            queue.append((v, c, j, k - 1))
            for (w, d) in bop(g, v(k), S):
                queue.append((w, c * d * beta, j, k - 1))
        else:
            queue.append((v, c, j, k - 1))
            queue.append((g, c * beta, j, k - 1))


def test_reformulation(n=4):
    for z in [Permutation(2,1,6,5,4,3)]: #Permutation.involutions(n):
        cyc = [(j, k) for j, k in z.get_two_cycles()] + [(j, j) for j in z.fixed(n + 1)]
        for j, k in cyc:
            print('n =', n, 'z =', z, 'j =', j, 'k =', k)

            f = Vector()
            for y, c in A(z, 1, j):
                f += invgroth(y) * c
            f = multiply_via_grothendieck_transitions(f, j)
            if j < k:
                f = multiply_via_grothendieck_transitions(f, k)
            try:
                decomposeinv(f)
            except:
                raise Exception

            g = Vector()
            for y, c in B(z, k, n + 2):
                g += invgroth(y) * c
            
            if f != g:
                try:
                    dec = decomposeinv(f - g)
                    for w, coeff in sorted(dec.items(), key=lambda a: a[0].involution_length()):
                        print('  ', w, coeff)
                except:
                    print('diff =', f - g)
                    return f - g
                print()
                input('?')
            print()