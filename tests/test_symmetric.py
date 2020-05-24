from symmetric import SchurQ, Vector
from signed import SignedPermutation


def test_schurq_multiplication():
    assert SchurQ(1) * SchurQ() == Vector({SchurQ(1): 1})
    assert SchurQ(1) * SchurQ(1) == 2 * SchurQ(2)
    assert SchurQ(1) * SchurQ(2) == SchurQ(2, 1) + 2 * SchurQ(3)
    assert SchurQ(1) * SchurQ(3) == SchurQ(3, 1) + 2 * SchurQ(4)
    assert SchurQ(1) * SchurQ(4) == SchurQ(4, 1) + 2 * SchurQ(5)
    assert SchurQ(1) * SchurQ(5) == SchurQ(5, 1) + 2 * SchurQ(6)

    assert SchurQ(2) * SchurQ(2) == 2 * SchurQ(3, 1) + 2 * SchurQ(4)
    assert SchurQ(2) * SchurQ(3) == SchurQ(3, 2) + 2 * SchurQ(4, 1) + 2 * SchurQ(5)
    assert SchurQ(2) * SchurQ(4) == SchurQ(4, 2) + 2 * SchurQ(5, 1) + 2 * SchurQ(6)
    assert SchurQ(2) * SchurQ(5) == SchurQ(5, 2) + 2 * SchurQ(6, 1) + 2 * SchurQ(7)

    assert SchurQ(3) * SchurQ(3) == 2 * SchurQ(4, 2) + 2 * SchurQ(5, 1) + 2 * SchurQ(6)
    assert SchurQ(3) * SchurQ(4) == SchurQ(4, 3) + 2 * SchurQ(5, 2) + 2 * SchurQ(6, 1) + 2 * SchurQ(7)
    assert SchurQ(3) * SchurQ(5) == SchurQ(5, 3) + 2 * SchurQ(6, 2) + 2 * SchurQ(7, 1) + 2 * SchurQ(8)

    assert SchurQ(4) * SchurQ(4) == 2 * SchurQ(5, 3) + 2 * SchurQ(6, 2) + 2 * SchurQ(7, 1) + 2 * SchurQ(8)
    assert SchurQ(4) * SchurQ(5) == SchurQ(5, 4) + 2 * SchurQ(6, 3) + 2 * SchurQ(7, 2) + 2 * SchurQ(8, 1) + 2 * SchurQ(9)

    assert SchurQ(5) * SchurQ(5) == 2 * SchurQ(6, 4) + 2 * SchurQ(7, 3) + 2 * SchurQ(8, 2) + 2 * SchurQ(9, 1) + 2 * SchurQ(10)


def conjectural_stanley_schur_decomposition(w, starting_rank=None, verbose=True, step=0):
    w = w.reduce()
    n = w.rank

    if starting_rank is None:
        starting_rank = n

    space = ((4 + starting_rank) * step) * ' '
    if w in SignedPermutation.longest_element(starting_rank - 1).get_atoms():
        print(space, '*', w, 'is atom')
        print()
        verbose = False
        yield (w, 1)
        return

    sh = w.increasing_shape()
    if sh is not None:
        return

    def get_shape(oneline):
        while oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        oneline = oneline[1:]
        ans = []
        while oneline:
            for i in range(len(oneline)):
                a = oneline[i]
                if i == 0 and a < 0:
                    ans += [(a, -a)]
                    oneline = oneline[1:]
                    break
                if i + 1 >= len(oneline):
                    continue
                b = oneline[i + 1]
                if 0 < a < -b:
                    ans += [(a, -b)]
                    oneline = oneline[:i] + oneline[i + 2:]
                    break
        return ans

    if verbose:
        r = 0
        s = n
        while r < s - 1 and (w(r + 1) == n or w(r + 1) == w(r) - 1):
            r += 1
        try:
            z = SignedPermutation(*w.oneline[r:]).reduce()
            a = SignedPermutation.longest_element(n - r - 1).get_atoms()
            assert z in a
        except:
            try:
                r += 1
                shape = get_shape(tuple(w(i) for i in range(r, n + 1)))
                s = w.inverse()(
                    max([a for a, b in shape if 0 < a < w(r) < b])
                )
                print('(r,s) = (%s,%s)' % (r, s))
            except:
                print(space, 'no inversions:', w, '\n')
                w, r, s = eval(input('w, r, s = '))
                print('(r,s) = (%s,%s)' % (r, s))
    else:
        r = w.last_descent()
        s = w.last_inversion(r)

    v = w * SignedPermutation.reflection_t(r, s, n)

    v_len = len(v)
    assert v_len + 1 == len(w)

    indices = [v * SignedPermutation.reflection_s(r, r, n)]
    indices += [v * SignedPermutation.reflection_t(i, r, n) for i in range(1, r)]
    newline = v.oneline + (n + 1,)
    v = SignedPermutation(*newline)
    indices += [v * SignedPermutation.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
    indices = [x.reduce() for x in indices if v_len + 1 == len(x)]

    subindices = [v * SignedPermutation.reflection_t(r, i, n + 1) for i in range(r + 1, n + 2) if i != s]
    subindices = [x.reduce() for x in subindices if v_len + 1 == len(x)]

    if verbose:
        print(space, w, '->', v.reduce(), '->', indices, ('- ' + str(subindices)) if subindices else '')
        print()

    for x in indices:
        for a, coeff in conjectural_stanley_schur_decomposition(starting_rank, verbose, step + 1):
            yield (a, coeff)

    for x in subindices:
        yield (x, -1)
