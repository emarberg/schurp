from partitions import StrictPartition
from permutations import Permutation
from vectors import Vector


involution_stanley_cache = {(): {(): 1}}


class InvStanleyExpander:

    @property
    def cache(self):
        return involution_stanley_cache

    def __init__(self, w):
        assert type(w) == Permutation
        self.w = w
        self._key = None

    def __repr__(self):
        return 'Expand[ %s ]' % str(self.w)

    @property
    def key(self):
        if self._key is None:
            k = self.w.oneline
            while k and k[-1] == len(k):
                k = k[:-1]
            self._key = tuple(k)
        return self._key

    def set_cache(self, vector):
        self.cache[self.key] = {
            tuple(q.mu.parts): v
            for q, v in vector.items()
        }

    def get_cache(self):
        return Vector({
            SchurP(StrictPartition(*k)): v
            for k, v in self.cache[self.key].items()
        })

    def expand(self):
        if self.key in self.cache:
            return self.get_cache()
        else:
            children = self.get_children()
            if children:
                ans = Vector()
                for u in children:
                    ans += InvStanleyExpander(u).expand()
            else:
                q = SchurP(self.get_grassmannian_shape())
                ans = Vector({q: 1})
            self.set_cache(ans)
            return ans

    def get_children(self):
        u, r, s = self.descend()
        if u is None:
            return []
        # sanity check
        # assert [w] == cls.phi_plus(u, r) and u(r) <= r
        children = InvStanleyExpander(u).phi_minus(u(r))
        if len(children) == 0:
            n = len(self.w.oneline)
            c = Permutation.cycle(range(1, n + 2))
            new_w = c**-1 * self.w * c
            return InvStanleyExpander(new_w).get_children()
        return children

    def descend(self):
        if self.is_grassmannian():
            return None, None, None
        n = len(self.w.oneline)
        for r in reversed(range(1, n + 1)):
            for s in reversed(range(r + 1, n + 1)):
                if self.w(s) <= min([self.w(r), r]):
                    y = self.tau_inverse(r, s)
                    return y, r, s

    def get_descents(self):
        n = len(self.w.oneline)
        for i in range(1, n + 1):
            if self.w(i) > self.w(i + 1) and (self.w(i + 1) <= i):
                yield i

    def is_grassmannian(self):
        return len(list(self.get_descents())) <= 1

    def get_grassmannian_shape(self):
        des = list(self.get_descents())
        assert len(des) <= 1
        n = des[0] if des else 0
        r = max(self.w.support()) - n if des else 0
        return StrictPartition(*[n + 1 - self.w(n + 1 + i) for i in range(r)])

    def tau(self, i, j):
        w = self.w
        a_tup = tuple(sorted(set([i, j, w(i), w(j)])))
        if len(a_tup) == 2 and w(i) == i:
            r = Permutation.cycle([i, j])
            return w * r
        elif len(a_tup) == 3:
            a, b, c = a_tup
            if (i, j) in [(b, c), (a, c)] and w(a) == b and w(c) == c:
                r = Permutation.cycle([b, c])
                return r * w * r
            elif (i, j) in [(a, b), (a, c)] and w(a) == a and w(b) == c:
                r = Permutation.cycle([a, b])
                return r * w * r
        elif len(a_tup) == 4:
            a, b, c, d = a_tup
            if (i, j) == (b, c) and w(a) == b and w(c) == d:
                r = Permutation.cycle([i, j])
                return r * w * r
            elif (i, j) in [(a, c), (b, d), (a, d)] and w(a) == b and w(c) == d:
                s = Permutation.cycle([a, b])
                r = Permutation.cycle([a, c])
                return r * w * s * r
            elif (i, j) in [(a, b), (c, d), (a, d)] and w(a) == c and w(b) == d:
                r = Permutation.cycle([a, b])
                return r * w * r
        return w

    def tau_inverse(self, i, j):
        w = self.w
        a_tup = tuple(sorted(set([i, j, w(i), w(j)])))
        if len(a_tup) == 2 and w(i) == j:
            r = Permutation.cycle([i, j])
            return w * r
        elif len(a_tup) == 3:
            a, b, c = a_tup
            if (i, j) in [(b, c), (a, c)] and w(a) == c and w(b) == b:
                r = Permutation.cycle([b, c])
                return r * w * r
            elif (i, j) in [(a, b), (a, c)] and w(a) == c and w(b) == b:
                r = Permutation.cycle([a, b])
                return r * w * r
        elif len(a_tup) == 4:
            a, b, c, d = a_tup
            if (i, j) == (b, c) and w(a) == c and w(b) == d:
                r = Permutation.cycle([i, j])
                return r * w * r
            elif (i, j) in [(a, c), (b, d), (a, d)] and w(a) == d and w(b) == b and w(c) == c:
                s = Permutation.cycle([a, b])
                r = Permutation.cycle([a, c])
                return r * w * r * s
            elif (i, j) in [(a, b), (c, d), (a, d)] and w(a) == d and w(b) == c:
                r = Permutation.cycle([a, b])
                return r * w * r
        return w

    def phi_plus(self, r):
        ans = []
        n = max(r, len(self.w.oneline))
        for i in range(r + 1, n + 2):
            tau = self.tau(r, i)
            if (tau.involution_length() - self.w.involution_length()) == 1:
                ans += [tau]
        return sorted(set(ans))

    def phi_minus(self, r):
        ans = []
        for i in range(1, r):
            tau = self.tau(i, r)
            if (tau.involution_length() - self.w.involution_length()) == 1:
                ans += [tau]
        return sorted(set(ans))


class SchurP:

    def __init__(self, mu):
        assert type(mu) == StrictPartition
        self.mu = mu

    # def skew(self, nu):
    #     assert type(nu) == StrictPartition
    #     assert self.mu.contains(nu)

    #     m = len(self.mu)
    #     n = len(nu)
    #     if (m + n) % 2 != 0:
    #         n += 1

        
    #     matrix = []

    def __hash__(self):
        return hash(self.mu)

    def __eq__(self, other):
        assert type(other) == type(self)
        return self.mu == other.mu

    def __lt__(self, other):
        assert type(other) == type(self)
        return self.mu < other.mu

    def __repr__(self):
        return 'P(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __add__(self, other):
        if type(other) == Vector:
            return Vector.base(self) + other
        else:
            assert type(other) == type(self)
            return Vector.base(self) + Vector.base(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == Vector:
            return Vector.base(self) - other
        else:
            assert type(other) == SchurP
            return Vector.base(self) - Vector.base(other)

    def __rsub__(self, other):
        return other.__sub__(Vector.base(self))

    def __neg__(self):
        return -Vector.base(self)

    def __mul__(self, other):
        assert type(other) == type(self)
        # return Vector({SchurP(p): v for p, v in self.mu.pieri(other.mu(1)).items()})
        u = self.mu.to_grassmannian()
        v = other.mu.to_grassmannian()
        n = len(u.oneline)
        w = Permutation(u.oneline + [i + n for i in v.oneline])
        return InvStanleyExpander(w).expand()

    def __rmul__(self, i):
        self.__mul__(i)


class SchurQ(SchurP):

    def __repr__(self):
        return 'Q(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __mul__(self, other):
        assert type(other) == type(self)
        mult = len(self.mu) + len(other.mu)
        u = self.mu.to_grassmannian()
        v = other.mu.to_grassmannian()
        n = len(u.oneline)
        w = Permutation(u.oneline + [i + n for i in v.oneline])
        vector = InvStanleyExpander(w).expand()
        return Vector({
            SchurQ(p.mu): 2**(mult - len(p.mu)) * v
            for p, v in vector.items()
        })
