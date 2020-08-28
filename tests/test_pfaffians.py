from polynomials import MPolynomial
from permutations import Permutation
from cached import cached_value
import itertools

OFFSET_CACHE = {}
PFAFFIAN_CACHE = {}
DETERMINANT_CACHE = {}
PFSGN_CACHE = {}
BLOCK_CACHE = {}


class PfPol(MPolynomial):

    @cached_value(OFFSET_CACHE)
    def to_pair(cls, i): # noqa
        assert i > 0
        j, t = 1, 1
        while t < i:
            j += 1
            t += j
        return (j - (t - i), 1 + (t - i))

    @cached_value(OFFSET_CACHE)
    def from_pair(cls, i, j): # noqa
        d = i + j - 2
        return (d * (d + 1)) // 2 + i

    @classmethod
    def monomial(cls, i, j, power=1):
        assert i > 0 and j > 0
        if i == j:
            return cls.zero()
        c = 1 if i > j else -1
        (i, j) = (i, j) if i > j else (j, i)
        return super(PfPol, cls).monomial(cls.from_pair(i, j), power) * c

    @classmethod
    def letters(cls, i):
        assert i > 0
        return "x_{%i,%i}" % cls.to_pair(i)

    @classmethod
    def sorter(cls, index):
        big = 10
        assert not any(max(cls.to_pair(p)) > big for p in index)
        return tuple(
            # # lexicographic term order
            # -index.get(cls.from_pair(i, j), 0)
            # for i in range(big + 1)
            # for j in range(i - 1, 0, -1)
            # reverse lexicographic term order
            index.get(cls.from_pair(i, j), 0)
            for i in range(big, 0, -1)
            for j in range(i - 1, 0, -1)
        )

    @cached_value(PFAFFIAN_CACHE)
    def pfminor(cls, a, b): # noqa
        verbose = False
        a = tuple(range(1, a + 1)) if type(a) == int else a
        b = tuple(range(1, b + 1)) if type(b) == int else b
        assert len(a) <= len(b)
        dot = tuple(zip(a, reversed(b)))
        a_tild = a  # tuple(s for (s, t) in dot if (t, s) not in dot)
        a_less = tuple(s for s in a_tild if s not in b)
        ans = cls.zero()
        if verbose:
            print('A =', a)
            print('B =', b)
            print('dot =', dot)
            print('A_tilde =', a_tild)
            print('A_tilde \\ B =', a_less)
            print()
        for k in range(len(a_less) + 1):
            for sub in itertools.combinations(a_less, k):
                i = tuple(s for s in a_tild if s not in sub)
                j = tuple(sorted(b + sub))
                if verbose:
                    print()
                    print('* S =', sub)
                    print('* A_tilde \\ S =', i)
                    print('* B cup S =', j)
                if len(i) % 2 != 0 or len(j) % 2 != 0:
                    continue
                c = (-1) ** (len(i) // 2) * cls._sgn_helper(a_tild, sub) * cls._sgn_helper(j, sub)
                v = cls.pfaffian(i) * cls.pfaffian(j)
                if verbose:
                    print('* c =', (-1) ** (len(i) // 2), '*', cls._sgn_helper(a_tild, sub), '*', cls._sgn_helper(j, b))
                    print('* v =', v)
                ans += c * v
        return ans

    @cached_value(PFAFFIAN_CACHE)
    def pfaffian(cls, i): # noqa
        i = list(range(1, i + 1)) if type(i) == int else i
        ans = cls.zero()
        for z in Permutation.fpf_involutions(len(i)):
            c = (-1) ** (z.fpf_involution_length() + len(i) // 2)
            v = cls.one()
            for t in range(1, len(i) + 1):
                if z(t) < t:
                    v *= cls.monomial(i[z(t) - 1], i[t - 1])
            ans += c * v
        return ans

    @cached_value(BLOCK_CACHE)
    def blockpf(cls, a, b): # noqa
        assert len(a) <= len(b)
        ans = cls.zero()
        n = len(a) + len(b)
        for z in Permutation.fpf_involutions(n):
            if any(len(b) < z(t) for t in range(len(b) + 1, n + 1)):
                continue
            c = (-1) ** (z.fpf_involution_length() + n // 2)
            v = cls.one()
            for t in range(1, n + 1):
                if z(t) < t:
                    i = b[z(t) - 1] if z(t) <= len(b) else a[z(t) - 1 - len(b)]
                    j = b[t - 1] if t <= len(b) else a[t - 1 - len(b)]
                    v *= cls.monomial(i, j)
            ans += c * v
        return ans

    @cached_value(PFSGN_CACHE)
    def _sgn_helper(cls, i, j): # noqa
        k = tuple(sorted(t for t in i if t not in j))

        mapping = {}
        for t in range(len(i)):
            mapping[i[t]] = k[t] if t < len(k) else j[t - len(k)]

        ans, base = 1, set(i)
        while base:
            c = [base.pop()]
            while mapping[c[-1]] != c[0]:
                c.append(mapping[c[-1]])
            base = base - set(c)
            if len(c) % 2 == 0:
                ans *= -1
        return ans

    @cached_value(DETERMINANT_CACHE)
    def _minor(cls, i, j): # noqa
        j = i if j is None else j
        i = list(range(1, i + 1)) if type(i) == int else i
        j = list(range(1, j + 1)) if type(j) == int else j
        assert len(i) == len(j)
        ans = cls.zero()
        for z in Permutation.all(len(i)):
            c = z.sgn()
            v = cls.one()
            for t in range(1, len(i) + 1):
                v *= cls.monomial(i[t - 1], j[z(t) - 1])
            ans += c * v
        return ans

    @classmethod
    def minor(cls, i, j=None):
        return cls._minor(i, j)

    @classmethod
    def determinant(cls, i, j=None):
        return cls._minor(i, j)


def transform(x, y):
    if len(x) < 2 or x[1] > y[0]:
        return x, y
    assert x[0] < y[0]
    i = [_ for _ in range(1, len(x)) if y[0] >= x[_]][-1]
    if y[i] <= x[0]:
        x, y = tuple(reversed(y[:i + 1])) + x[i + 1:], tuple(reversed(x[:i + 1])) + y[i + 1:]
    newx, newy = transform(x[1:i + 1], y[1:i + 1])
    return (x[0],) + newx + x[i + 1:], (y[0],) + newy + y[i + 1:]


def twists(x, y):
    ans = 0
    for i in range(len(x)):
        for j in range(i + 1, len(y)):
            b, c = x[i], x[j]
            d, a = y[i], y[j]
            if a < b < c < d:
                ans += 1
    return ans


def test_relation(q=7):
    t = tuple(range(1, 2 * q + 1))
    for x in itertools.combinations(t, q):
        for y in itertools.combinations(tuple(reversed(t)), q):
            m = twists(x, y)
            tx, ty = transform(x, y)
            if m > 0:
                print('A =', x)
                print('B =', y)
                print()
                print('twists =', m)
                print()
                print('transformed A =', tx)
                print('transformed B =', ty)
                print()
                print('twists =', twists(tx, ty))
                print()
                print()
                print()
            assert all(tx[i] < tx[i + 1] for i in range(q - 1))
            assert all(ty[i] > ty[i + 1] for i in range(q - 1))
            assert twists(tx, ty) == 0


def test_leading(m=2):
    r = tuple(range(1, 1 + 2 * m))
    for k in range(m, 2 * m):
        for b in itertools.combinations(r, k):
            b = tuple(reversed(b))
            for a in itertools.combinations(r, k):
                dot = tuple(zip(a, b))
                if any(x == y for (x, y) in dot):
                    continue
                pairs = {(x, y) if x > y else (y, x) for (x, y) in dot}
                a_tilde = tuple(x for (x, y) in dot if (y, x) not in dot)

                t = tuple(i if i in a_tilde else 0 for i in a)

                print('A:', a)
                print('B:', b)
                print()
                print('dot =', dot)
                print('pairs =', pairs)
                print()
                print('A:', t)
                print('B:', b)

                u = PfPol.one()
                for (x, y) in pairs:
                    u *= PfPol.monomial(x, y)

                block = PfPol.blockpf(a_tilde, b)
                f = PfPol({list(block)[0]: 1})

                print()
                print('f_AB =', block)
                print()
                print('u =', u)
                print('f =', f)
                print()
                print('twists =', twists(a, b))
                assert u == f or twists(a, b) > 0
                print()
                print()
                print()


def test_blockfp(m=7):
    r = tuple(range(1, 1 + m))
    for k in range(m):
        for l in range(k + 1):
            for b in itertools.combinations(r, k):
                for a in itertools.combinations(r, l):
                    dot = tuple(zip(a, reversed(b)))
                    if any((y, x) in dot for (x, y) in dot):
                        continue
                    g = PfPol.pfminor(a, b)
                    block = PfPol.blockpf(a, b)
                    print(a)
                    print(tuple(reversed(b)))
                    # print(g)
                    print(block)
                    print(MPolynomial(block.coeffs))
                    # print(g == block)
                    print()
                    assert (g == block)


def test_pfminor(m=7):
    for i in range(m):
        a = tuple(t + 1 for t in range(i))
        b = tuple(i + t for t in a)
        print(a, b)
        # print(PfPol.pfminor(a, b))
        # print(PfPol.minor(a, b))

        f = PfPol.pfminor(a, b)
        g = PfPol.minor(a, b)
        assert f == g or f == -g
        print('det:', f == g)

        dot = tuple(zip(a, reversed(a)))
        c = tuple(x for (x, y) in dot if (y, x) not in dot)
        f = PfPol.pfminor(c, a)
        g = PfPol.pfaffian(a)
        assert f == g or f == -g
        print(' pf:', f == g)
        print()


def test_pfaffian(m=7):
    for i in range(m):
        print(i)
        assert PfPol.pfaffian(i) ** 2 == PfPol.minor(i)


def test_offsets(m=100):
    for i in range(1, m):
        print(i, PfPol.to_pair(i), PfPol.from_pair(*PfPol.to_pair(i)))
        assert PfPol.from_pair(*PfPol.to_pair(i)) == i
