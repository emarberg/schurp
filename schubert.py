from vectors import Vector
from permutations import Permutation
from signed import SignedPermutation
from polynomials import (
    MPolynomial,
    one,
    x,
    y,
    X,
)
from tableaux import Tableau
from stable.utils import GQ_expansion, GP_expansion, SymmetricPolynomial


SCHUBERT_CACHE = {}
DOUBLE_SCHUBERT_CACHE = {}

GROTHENDIECK_CACHE = {}
DOUBLE_GROTHENDIECK_CACHE = {}

FPF_SCHUBERT_CACHE = {}
FPF_GROTHENDIECK_CACHE = {}

INV_SCHUBERT_CACHE = {}
INV_GROTHENDIECK_CACHE = {}
ALT_INV_GROTHENDIECK_CACHE = {}

C_GROTHENDIECK_CACHE = {}
C_DOUBLE_GROTHENDIECK_CACHE = {}

B_GROTHENDIECK_CACHE = {}
B_DOUBLE_GROTHENDIECK_CACHE = {}

D_GROTHENDIECK_CACHE = {}
D_DOUBLE_GROTHENDIECK_CACHE = {}

C_SYMMETRIC_GROTHENDIECK_SIMPLE = {}
B_SYMMETRIC_GROTHENDIECK_SIMPLE = {}
D_SYMMETRIC_GROTHENDIECK_SIMPLE = {}

C_SYMMETRIC_GROTHENDIECK_CACHE = {}
B_SYMMETRIC_GROTHENDIECK_CACHE = {}
D_SYMMETRIC_GROTHENDIECK_CACHE = {}


class AbstractSchubert(object):

    @classmethod
    def double_lettering(cls):
        def ans(i):
            if i > 0 and i % 2 != 0:
                return "x_" + str((i + 1) // 2)
            elif i > 0 and i % 2 == 0:
                return "y_" + str(i // 2)
            elif i == 0:
                return "\u03B2"
            else:
                return "z_" + str(-i)
        return ans

    @classmethod
    def cache(cls):
        raise NotImplementedError

    @classmethod
    def top(cls, w):
        raise NotImplementedError

    @classmethod
    def invalid_case(cls):
        raise Exception

    @classmethod
    def is_valid(cls, w):
        return True

    @classmethod
    def get_ascent(cls, w):
        raise NotImplementedError

    @classmethod
    def divided_difference(cls, f, i):
        return f.divided_difference(i)

    @classmethod
    def reduce(cls, oneline):
        return tuple(oneline)

    @classmethod
    def get(cls, w, verbose=False):
        assert cls.is_valid(w)
        oneline = cls.reduce(w.oneline)
        cache = cls.cache()
        if oneline not in cache:
            v, i = cls.get_ascent(Permutation(*oneline))
            if i is not None:
                s = cls.get(v)
                s = cls.divided_difference(s, i)
            elif v is not None:
                s = cls.top(v)
            else:
                s = cls.invalid_case()
            cache[oneline] = s
            if verbose:
                print(' . . .', cls.__name__, 'cache:', len(cache))
        return cache[oneline]

    @classmethod
    def from_code(cls, code):
        raise NotImplementedError

    @classmethod
    def least_term(cls, f):
        return min(f, key=lambda s: tuple(-i for k in sorted(s) for i in s[k] * [k]))

    @classmethod
    def decompose(cls, f):
        ans = Vector()
        while not f.is_zero():
            m = cls.least_term(f)
            code = max(m) * [0] if m else []
            for i, v in m.items():
                code[i - 1] = v
            code = tuple(code)
            w = cls.from_code(code)
            ans += Vector({w: f[m]})
            f = f - f[m] * cls.get(w)
        return ans


class Schubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return SCHUBERT_CACHE

    @classmethod
    def reduce(cls, oneline):
        while oneline and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        return tuple(oneline)

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if n == 0:
            return w.s_i(1), 1
        if w.is_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        return w * w.s_i(i), i

    @classmethod
    def from_code(cls, code):
        return Permutation.from_code(code)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.rothe_diagram():
            s *= x(i)
        return s


class Grothendieck(Schubert):

    beta = -1
    # beta = x(0)

    @classmethod
    def cache(cls):
        return GROTHENDIECK_CACHE

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)


class DoubleSchubert(Schubert):

    @classmethod
    def cache(cls):
        return DOUBLE_SCHUBERT_CACHE

    @classmethod
    def top(cls, w):
        s = one()
        for i, j in w.rothe_diagram():
                s *= (x(i) - y(j))
        return s


class DoubleGrothendieck(Grothendieck):

    @classmethod
    def cache(cls):
        return DOUBLE_GROTHENDIECK_CACHE

    @classmethod
    def top(cls, w):
        s = one()
        for i, j in w.rothe_diagram():
                s *= (x(i) + y(j) + cls.beta * x(i) * y(j))
        return s

    @classmethod
    def expand_double_reflection_chain(cls, start, chain, rank):
        length = lambda x: x.length()
        n = rank

        def act(w, a):
            return tuple(a[w.inverse()(i + 1) - 1] for i in range(n))

        def add(a, b):
            return tuple(a[i] + b[i] for i in range(n))

        def negate(a):
            return tuple(-a[i] for i in range(n))

        def alpha(i, j):
            ans = n * [0]
            ans[i - 1] = 1
            ans[j - 1] = -1
            return tuple(ans)

        rho = tuple(n - i for i in range(n))
        one = tuple(1 for i in range(n))

        def reduce(a):
            while min(a) < 0:
                a = add(one, a)
            while min(a) > 0:
                a = add(negate(one), a)
            assert all(v % n == 0 for v in a)
            a = [v // n for v in a]
            ans = y(0)**0
            for i in range(n):
                ans *= y(i + 1)**a[i]
            return ans

        ans = Vector()
        queue = [(1, act(start, negate(rho)), start, chain)]
        while queue:
            (sgn, bterm, z, c) = queue[0]
            queue = queue[1:]
            if len(c) == 0:
                bterm = add(bterm, act(z, rho))
                ans += Vector({z: sgn * reduce(bterm)})
            else:
                (i, j) = c[0]
                t = Permutation.t_ij(i, j)
                a = act(z, alpha(i, j))
                queue.append((sgn, add(a, bterm), z, c[1:]))
                zt = z * t
                if length(zt) == length(z) - 1:
                    queue.append((sgn if i < j else -sgn, bterm, zt, c[1:]))
        return ans


class FPFSchubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return FPF_SCHUBERT_CACHE

    @classmethod
    def reduce(cls, oneline):
        while len(oneline) >= 2 and oneline[-2] == len(oneline) and oneline[-1] == len(oneline) - 1:
            oneline = oneline[:-2]
        return tuple(oneline)

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.fpf_rothe_diagram():
            s *= x(i) + x(j)
        return s

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_fpf_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        s = w.s_i(i)
        return s * w * s, i

    @classmethod
    def is_valid(cls, w):
        return all(len(c) == 2 for c in w.cycles)

    @classmethod
    def from_code(cls, code):
        return Permutation.from_fpf_involution_code(code)


class FPFGrothendieck(FPFSchubert):

    beta = 1

    @classmethod
    def cache(cls):
        return FPF_GROTHENDIECK_CACHE

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.fpf_rothe_diagram():
            s *= (x(i) + x(j) + cls.beta * x(i) * x(j))
        return s

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)


class InvSchubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return INV_SCHUBERT_CACHE

    @classmethod
    def top(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= (x(i) + x(j)) if i != j else x(i)
        return s


    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        s = w.s_i(i)
        if s * w == w * s:
            return w * s, i
        else:
            return s * w * s, i

    @classmethod
    def is_valid(cls, w):
        return w.is_involution()

    @classmethod
    def from_code(cls, code):
        return Permutation.from_involution_code(code)


    @classmethod
    def vexillary_tableau_formula(cls, w):
        path = w.get_essential_path()
        
        t = [0, 0]
        for i in range(1, len(path)):
            a1, b1 = path[i - 1]
            a2, b2 = path[i]
            assert a1 == a2 or b1 == b2
            if a1 == a2:
                assert b2 == b1 - 1
                t.append(-X(b1))
            elif b1 == b2:
                assert a2 == a1 + 1
                t.append(X(a2))
        max_entry = path[0][0]
        mu = w.involution_shape()

        ans = 0
        for tab in Tableau.get_semistandard_shifted(mu, n=max_entry, diagonal_primes=True):
            summand = MPolynomial.one()
            for (i, j) in tab.mapping:
                a = abs(tab.get(i, j))
                b = j - i + 1
                s = -1 if tab.get(i, j).is_primed() else 1
                summand *= X(a) + s * t[b]
            ans += summand
        return ans

class InvGrothendieck(InvSchubert):

    beta = 1

    @classmethod
    def grothendieck_involution_words(cls, w):
        assert w.is_vexillary()
        f = cls.get(w)
        dec = Grothendieck.decompose(f)
        for u in dec:
            for word in u.get_reduced_words():
                yield word, dec[u]

    @classmethod
    def cache(cls):
        return INV_GROTHENDIECK_CACHE

    @classmethod
    def invalid_case(cls):
        return MPolynomial.zero()

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_dominant():
            return w, None
        if not w.is_vexillary():
            return None, None
        for i in set(range(1, n)) - w.right_descent_set:
            s = w.s_i(i)
            ws = w * s if  s * w == w * s else s * w *s
            if ws.is_vexillary():
                return ws, i
        return None, None

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= x(i) + x(j) + cls.beta * x(i) * x(j)
        return s

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)

    @classmethod
    def diagonal_product(cls, n):
        s = one()
        for i in range(1, n + 1):
            s *= (x(i) + x(i) + cls.beta * x(i) * x(i))
        return s

    @classmethod
    def diagonal_subproduct(cls, n):
        s = one()
        for i in range(1, n + 1):
            s *= (2 + cls.beta * x(i))
        return s


class AltInvGrothendieck(InvGrothendieck):

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_dominant():
            return w, None
        for i in set(range(1, n)) - w.right_descent_set:
            s = w.s_i(i)
            ws = w * s if  s * w == w * s else s * w *s
            return ws, i
        return None, None

    @classmethod
    def cache(cls):
        return ALT_INV_GROTHENDIECK_CACHE

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= (x(i) + x(j) + cls.beta * x(i) * x(j)) if i != j else x(i)
        return s


class GrothendieckC(Grothendieck):

    @classmethod
    def expand_reflection_chain(cls, start, chain, length):
        ans = Vector()
        queue = [(1, start, chain)]
        while queue:
            (sgn, z, c) = queue[0]
            queue = queue[1:]
            if len(c) == 0:
                ans += Vector({z: sgn})
            else:
                if len(c[0]) == 2:
                    a = 1
                    t, b = c[0]
                else:
                    t, a, b = c[0]
                queue.append((a * sgn, z, c[1:]))
                zt = z * t
                if length(zt) == length(z) + 1:
                    queue.append((b * sgn, zt, c[1:]))
        return ans

    @classmethod
    def get_hecke_compatible_sequences(cls, n, w, ell):
        return w.get_hecke_compatible_sequences_c(n, ell)

    @classmethod
    def beta_exponent(cls, a, w):
        return len(a) - len(w)

    @classmethod
    def exponent(cls, a, b):
        gamma = 0
        for i in range(len(a) - 1):
            if a[i] == a[i + 1] and b[i] == b[i + 1]:
                gamma += 1
        return len(set(b)) - gamma

    @classmethod
    def symmetric_is_zero(cls, n, w):
        return (n == 0 and len(w) > 0) or (0 < n < w.min_peaks() + 1)

    @classmethod
    def symmetric_with_unknown_ell(cls, n, w, x_not_y):
        ans = cls.symmetric(n, w, x_not_y, len(w))
        assert ans != 0
        ell = len(w) + 1
        while True:
            bns = cls.symmetric(n, w, x_not_y, ell)
            if ans == bns:
                break
            ans = bns
            ell += 1
        return ans

    @classmethod
    def symmetric(cls, n, w, x_not_y=True, ell=None):
        w = w.reduce()
        assert cls.is_valid(w)

        oneline = cls.reduce(w.oneline)
        key = (n, oneline, x_not_y, ell)
        cache = cls.symmetric_cache()
        
        if key not in cache:
            if cls.symmetric_is_zero(n, w):
                ans = MPolynomial.zero()
            elif ell is None:
                ans = cls.symmetric_with_unknown_ell(n, w, x_not_y)
            else:
                ansdict = {}
                dictionary = cls.get_hecke_compatible_sequences(n, w, ell)
                for a in dictionary:
                    for b in dictionary[a]:
                        term = one() * cls.beta**cls.beta_exponent(a, w) * 2**(len(b) + cls.exponent(a, b))
                        for i in b:
                            term *= x(i) if x_not_y else y(i)
                        ansdict[len(b)] = ansdict.get(len(b), 0) + term
                ans = sum([ansdict[i] // 2**i for i in ansdict])
            cache[key] = ans
        return cache[key]

    @classmethod
    def cache(cls):
        return C_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_cache(cls):
        return C_SYMMETRIC_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_simple(cls, w, verbose=False):
        cache = C_SYMMETRIC_GROTHENDIECK_SIMPLE

        w = w.reduce()
        n = w.rank
        key = tuple(w.oneline)
        
        if key not in cache:
            r = [i for i in range(1, n) if w(i) > w(i + 1)]
            if len(r) == 0:
                ans = MPolynomial.one()
                for i in range(1, n + 1):
                    if w(i) > 0:
                        break
                    ans *= y(i)**(-w(i))
            else:
                r = max(r)
                s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
                v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

                chain = []
                chain += [(SignedPermutation.reflection_s(i, r, n + 1), cls.beta) for i in range(n + 1, 0, -1) if i != r]
                chain += [(SignedPermutation.reflection_s(r, r, n + 1), cls.beta)]
                chain += [(SignedPermutation.reflection_t(i, r, n + 1), cls.beta) for i in range(1, r)]
        
                vec = Vector({v: 1}) - cls.expand_reflection_chain(v, chain, lambda x: x.length())
                vec *= -(cls.beta * MPolynomial.one())**-1

                ans = MPolynomial.zero()
                for (u, c) in vec.dictionary.items():
                    ans += cls.symmetric_simple(u) * c

                if verbose:
                    print(' . . .', cls.__name__, 'cache:', len(cache))
            cache[key] = ans
        
        return cache[key]

    @classmethod
    def get(cls, w, simple=True, verbose=False,):
        n = w.rank
        assert cls.is_valid(w)
        key = tuple(w.oneline)
        cache = cls.cache()
        if key not in cache:
            ans = MPolynomial.zero()
            for (u, v) in w.get_demazure_factorizations():
                sym = cls.symmetric_simple(u) if simple else cls.symmetric(n, u, x_not_y=False)
                ans += cls.beta**(u.length() + v.length() - w.length()) * sym * Grothendieck.get(v)
            cache[key] = ans
            if verbose:
                print(' . . .', cls.__name__, 'cache:', len(cache))
        return cache[key]


class GrothendieckB(GrothendieckC):

    @classmethod
    def exponent(cls, a, b):
        gamma = 0
        oB = 0
        for i in range(len(a)):
            if i < len(a) - 1 and a[i] == a[i + 1] and b[i] == b[i + 1]:
                gamma += 1
            if a[i] == 0:
                oB += 1
        return len(set(b)) - gamma - oB

    @classmethod
    def get_hecke_compatible_sequences(cls, n, w, ell):
        return w.get_hecke_compatible_sequences_b(n, ell)

    @classmethod
    def cache(cls):
        return B_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_cache(cls):
        return B_SYMMETRIC_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_simple(cls, w, verbose=False):
        cache = B_SYMMETRIC_GROTHENDIECK_SIMPLE

        w = w.reduce()
        n = w.rank
        key = tuple(w.oneline)
        
        if key not in cache:
            r = [i for i in range(1, n) if w(i) > w(i + 1)]
            if len(r) == 0:
                ans = MPolynomial.one()
                for i in range(1, n + 1):
                    if w(i) > 0:
                        break
                    ans *= y(i)**(-w(i))
            else:
                r = max(r)
                s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
                v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

                chain = []
                chain += [(SignedPermutation.reflection_s(r, r, n + 1), cls.beta)]
                chain += [(SignedPermutation.reflection_s(i, r, n + 1), cls.beta) for i in range(n + 1, 0, -1) if i != r]
                chain += [(SignedPermutation.reflection_s(r, r, n + 1), cls.beta)]
                chain += [(SignedPermutation.reflection_t(i, r, n + 1), cls.beta) for i in range(1, r)]
            
                vec = Vector({v: 1}) - cls.expand_reflection_chain(v, chain, lambda x: x.length())
                vec *= -(cls.beta * MPolynomial.one())**-1

                ans = MPolynomial.zero()
                for (u, c) in vec.dictionary.items():
                    ans += cls.symmetric_simple(u) * c

                if verbose:
                    print(' . . .', cls.__name__, 'cache:', len(cache))
            cache[key] = ans
        
        return cache[key]


class GrothendieckD(GrothendieckC):

    @classmethod
    def get_hecke_compatible_sequences(cls, n, w, ell):
        return w.get_hecke_compatible_sequences_d(n, ell)

    @classmethod
    def beta_exponent(cls, a, w):
        return len(a) - w.dlength()

    @classmethod
    def exponent(cls, a, b):
        gamma = 0
        oD = 0
        for i in range(len(a)):
            if i < len(a) - 1 and a[i] == a[i + 1] and b[i] == b[i + 1]:
                gamma += 1
            if abs(a[i]) == 1:
                oD += 1
        return len(set(b)) - gamma - oD

    @classmethod
    def symmetric_is_zero(cls, n, w):
        return (not w.is_even_signed()) or (n == 0 and len(w) > 0) or (0 < n < w.dtype_min_peaks() + 1)

    @classmethod
    def cache(cls):
        return D_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_cache(cls):
        return D_SYMMETRIC_GROTHENDIECK_CACHE

    @classmethod
    def symmetric_simple(cls, w, verbose=False):
        cache = D_SYMMETRIC_GROTHENDIECK_SIMPLE
        
        w = w.reduce()
        n = w.rank
        key = tuple(w.oneline)
        
        if key not in cache:
            r = [i for i in range(1, n) if w(i) > w(i + 1)]
            if len(r) == 0:
                ans = MPolynomial.one()
                for i in range(1, n + 1):
                    if w(i) >= -1:
                        break
                    ans *= y(i)**(-w(i) - 1)
            else:
                r = max(r)
                s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
                v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

                chain = []
                chain += [(SignedPermutation.reflection_s(i, r, n + 1), cls.beta) for i in range(n + 1, 0, -1) if i != r]
                chain += [(SignedPermutation.reflection_t(i, r, n + 1), cls.beta) for i in range(1, r)]
        
                vec = Vector({v: 1}) - cls.expand_reflection_chain(v, chain, lambda x: x.dlength())
                vec *= -(cls.beta * MPolynomial.one())**-1

                ans = MPolynomial.zero()
                for (u, c) in vec.dictionary.items():
                    ans += cls.symmetric_simple(u) * c

                if verbose:
                    print(' . . .', cls.__name__, 'cache:', len(cache))
            cache[key] = ans
        
        return cache[key]


class DoubleGrothendieckMixin:

    @classmethod
    def get(cls, w, simple=True, verbose=False,):
        n = w.rank
        assert cls.is_valid(w)
        key = tuple(w.oneline)
        cache = cls.cache()
        if key not in cache:
            ans = MPolynomial.zero()
            for (x, v2) in w.get_demazure_factorizations():
                for (u_inv, v1_inv) in x.inverse().get_demazure_factorizations():
                    u = u_inv.inverse()
                    sym = cls.symmetric_simple(u) if simple else cls.symmetric(n, u, x_not_y=False)
                    apart = Grothendieck.get(v2).odd_split_vars()
                    bpart = Grothendieck.get(v1_inv).even_split_vars()
                    ans += cls.beta**(u.length() + v1_inv.length() + v2.length() - w.length()) * sym * apart * bpart
            cache[key] = ans
            if verbose:
                print(' . . .', cls.__name__, 'cache:', len(cache))
        return cache[key]


class DoubleGrothendieckB(DoubleGrothendieckMixin,GrothendieckB):

    @classmethod
    def cache(cls):
        return B_DOUBLE_GROTHENDIECK_CACHE

    @classmethod
    def expand_double_reflection_chain(cls, start, chain, rank):
        length = lambda x: x.length()
        n = rank

        def act(w, a):
            return tuple(a[abs(w.inverse()(i + 1)) - 1] * (-1 if w.inverse()(i + 1) < 0 else 1) for i in range(n))

        def add(a, b):
            return tuple(a[i] + b[i] for i in range(n))

        def negate(a):
            return tuple(-a[i] for i in range(n))

        def alpha(i, j):
            ans = n * [0]
            ans[abs(i) - 1] += 2 if i > 0 else -2
            if i != j:
                ans[abs(j) - 1] += 2 if j > 0 else -2
            return tuple(ans)

        rho = tuple(reversed([2 * n - 1 - 2 * i for i in range(n)]))
        w0 = SignedPermutation.longest_element(n)

        def reduce(a):
            assert all(v % (4 * n) == 0 for v in a)
            a = [v // (4 * n) for v in a]
            ans = y(0)**0
            for i in range(n):
                ans *= y(i + 1)**a[i]
            return ans

        ans = Vector()
        queue = [(1, act(start, negate(rho)), start, chain)]
        while queue:
            (sgn, bterm, z, c) = queue[0]
            queue = queue[1:]
            if len(c) == 0:
                bterm = add(bterm, act(z, rho))
                ans += Vector({z: sgn * reduce(bterm)})
            else:
                (i, j) = c[0]
                assert abs(i) <= abs(j) and i + j != 0
                t = SignedPermutation.reflection_s(abs(i), abs(j), rank) if i * j > 0 else SignedPermutation.reflection_t(abs(i), abs(j), rank)
                a = act(z, alpha(i, j))
                queue.append((sgn, add(a, bterm), z, c[1:]))
                zt = z * t
                if length(zt) == length(z) - 1:
                    queue.append((sgn if j > 0 else -sgn, bterm, zt, c[1:]))
        return ans


class DoubleGrothendieckC(DoubleGrothendieckMixin,GrothendieckC):

    @classmethod
    def cache(cls):
        return C_DOUBLE_GROTHENDIECK_CACHE

    @classmethod
    def expand_double_reflection_chain(cls, start, chain, rank):
        length = lambda x: x.length()
        n = rank

        def act(w, a):
            return tuple(a[abs(w.inverse()(i + 1)) - 1] * (-1 if w.inverse()(i + 1) < 0 else 1) for i in range(n))

        def add(a, b):
            return tuple(a[i] + b[i] for i in range(n))

        def negate(a):
            return tuple(-a[i] for i in range(n))

        def alpha(i, j):
            ans = n * [0]
            ans[abs(i) - 1] += 1 if i > 0 else -1
            ans[abs(j) - 1] += 1 if j > 0 else -1
            return tuple(ans)

        rho = tuple(reversed([n - i for i in range(n)]))

        def reduce(a):
            assert all(v % (2 * n) == 0 for v in a)
            a = [v // (2 * n) for v in a]
            ans = y(0)**0
            for i in range(n):
                ans *= y(i + 1)**a[i]
            return ans

        ans = Vector()
        queue = [(1, act(start, negate(rho)), start, chain)]
        while queue:
            (sgn, bterm, z, c) = queue[0]
            queue = queue[1:]
            if len(c) == 0:
                bterm = add(bterm, act(z, rho))
                bterm = negate(bterm)
                ans += Vector({z: sgn * reduce(bterm)})
            else:
                (i, j) = c[0]
                assert abs(i) <= abs(j) and i + j != 0
                t = SignedPermutation.reflection_s(abs(i), abs(j), rank) if i * j > 0 else SignedPermutation.reflection_t(abs(i), abs(j), rank)
                a = act(z, alpha(i, j))
                queue.append((sgn, add(a, bterm), z, c[1:]))
                zt = z * t
                if length(zt) == length(z) + 1:
                    queue.append((sgn if j > 0 else -sgn, bterm, zt, c[1:]))
        return ans


class DoubleGrothendieckD(DoubleGrothendieckMixin,GrothendieckD):

    @classmethod
    def cache(cls):
        return D_DOUBLE_GROTHENDIECK_CACHE

    @classmethod
    def expand_double_reflection_chain(cls, start, chain, rank):
        length = lambda x: x.dlength()
        n = rank

        def act(w, a):
            return tuple(a[abs(w.inverse()(i + 1)) - 1] * (-1 if w.inverse()(i + 1) < 0 else 1) for i in range(n))

        def add(a, b):
            return tuple(a[i] + b[i] for i in range(n))

        def negate(a):
            return tuple(-a[i] for i in range(n))

        def alpha(i, j):
            ans = n * [0]
            ans[abs(i) - 1] += 1 if i > 0 else -1
            ans[abs(j) - 1] += 1 if j > 0 else -1
            return tuple(ans)

        rho = tuple(reversed([n - 1 - i for i in range(n)]))
        w0 = SignedPermutation.dtype_longest_element(n)

        def reduce(a):
            assert all(v % (2 * n - 2) == 0 for v in a)
            a = [v // (2 * n - 2) for v in a]
            ans = y(0)**0
            for i in range(n):
                ans *= y(i + 1)**a[i]
            return ans

        ans = Vector()
        queue = [(1, act(start, negate(rho)), start, chain)]
        while queue:
            (sgn, bterm, z, c) = queue[0]
            queue = queue[1:]
            if len(c) == 0:
                bterm = add(bterm, act(z, rho))
                bterm = act(w0, bterm)
                ans += Vector({z: sgn * reduce(bterm)})
            else:
                (i, j) = c[0]
                assert abs(i) < abs(j) and i + j != 0
                t = SignedPermutation.reflection_s(abs(i), abs(j), rank) if i * j > 0 else SignedPermutation.reflection_t(abs(i), abs(j), rank)
                a = act(z, alpha(i, j))
                queue.append((sgn, add(a, bterm), z, c[1:]))
                zt = z * t
                if length(zt) == length(z) + 1:
                    queue.append((sgn if j > 0 else -sgn, bterm, zt, c[1:]))
        return ans
