from signed import SignedPermutation
from permutations import Permutation
import itertools


CLAN_WORDS_CACHE = {}
CLAN_ATOMS_CACHE = {}


class Clan:

    TYPE_A = 'clans for G = GL(p + q), K = GL(p) x GL(q)'
    TYPE_B = 'clans for G = SO(2n + 1), K = S(O(2p) x O(2q + 1))'
    TYPE_C1 = 'clans for G = Sp(2n), K = GL(n)'
    TYPE_C2 = 'clans for G = Sp(2n), K = Sp(2p) x Sp(2q)'

    def __init__(self, oneline, family=TYPE_A):
        pairs = [
            (i, j)
            for i in range(len(oneline))
            for j in range(i + 1, len(oneline))
            if type(oneline[i]) == int and type(oneline[j]) == int and oneline[i] == oneline[j]
        ]
        oneline = list(oneline)
        for i, j in pairs:
            oneline[i] = j + 1
            oneline[j] = i + 1
        self.oneline = tuple(oneline)
        self.family = family

    @classmethod
    def create_a(cls, oneline):
        return Clan(oneline, cls.TYPE_A)

    @classmethod
    def create_b(cls, oneline):
        return Clan(oneline, cls.TYPE_B)

    @classmethod
    def create_c1(cls, oneline):
        return Clan(oneline, cls.TYPE_C1)

    @classmethod
    def create_c2(cls, oneline):
        return Clan(oneline, cls.TYPE_C2)

    def __repr__(self):
        l = [str(min(i, self(i))) if type(i) == int else '+' if i else '-' for i in self.oneline]
        return ' '.join(l)

    def __eq__(self, other):
        assert type(other) == Clan and self.family == other.family
        return self.oneline == other.oneline

    def __hash__(self):
        return hash((self.family, self.oneline))

    def is_matchless(self):
        return not any(type(i) == int for i in self.oneline)

    def clan_type(self):
        p = len([x for x in self.oneline if x is True])
        q = len([x for x in self.oneline if x is False])
        return p - q

    def is_aligned(self, matching, verbose=False):
        for (i, j) in matching:
            if self.family in [self.TYPE_B, self.TYPE_C2] and i + j == 0:
                continue
            if self.family == self.TYPE_A:
                pair = (self(i), self(j))
            elif self.family == self.TYPE_B:
                i, j = self.rank() + 1 + i, self.rank() + 1 + j
                pair = (self(i), (self(j)))
            elif self.family in [self.TYPE_C1, self.TYPE_C2]:
                i = self.rank() + i + (0 if i > 0 else 1)
                j = self.rank() + j + (0 if j > 0 else 1)
                pair = (self(i), self(j))
            else:
                raise Exception
            assert type(pair[0]) == bool and type(pair[1]) == bool
            if pair in [(True, True), (False, False)]:
                return False

        trivial = [a for (a, b) in sorted(matching) if -a == b]

        if self.family == self.TYPE_C2:
            return len(trivial) == abs(self.clan_type()) // 2
        if self.family == self.TYPE_B:
            k = abs(self.clan_type()) // 2
            if len(trivial) < k:
                return False
            signs = [self(self.rank() + 1 + i) for i in trivial]
            signs += [self(self.rank() + 1)]
            if any(signs[i] == signs[i + 1] for i in range(k, len(signs) - 1)):
                return False
        return True

    @classmethod
    def _all_a(cls, p, q):
        for w in Permutation.involutions(p + q):
            fixed = [i - 1 for i in range(1, p + q + 1) if w(i) == i]
            if len(fixed) < abs(p - q) or len(fixed) % 2 != (p - q) % 2:
                continue
            base = [i if i < w(i) else w(i) if w(i) < i else False for i in range(1, p + q + 1)]
            for subset in itertools.combinations(fixed, (p - q + len(fixed)) // 2):
                oneline = base[:]
                for i in subset:
                    oneline[i] = True
                cl = Clan(oneline, cls.TYPE_A)
                assert cl.clan_type() == p - q
                yield cl

    @classmethod
    def _all_b(cls, p, q):
        n = p + q
        for w in SignedPermutation.involutions(n):
            fixed = [i for i in range(1, n + 1) if w(i) == i]
            f = len(fixed)

            base = [i if i < w(i) else w(i) if w(i) < i else False for i in range(-n, 0)]
            base += [(2 * p - 2 * q + 2 * f - 1) % 4 == 1]
            base += [i if i < w(i) else w(i) if w(i) < i else False for i in range(1, n + 1)]

            e = 1 if base[n] else -1
            f = len(fixed)
            # 2 * k + e - 2 * (f - k) == 4 * k - 2 * f + e == 2 * p - 2 * q - 1
            # 4 * k + e == 2 * p - 2 * q + 2 * f - 1
            k = (2 * p - 2 * q + 2 * f - 1 - e)
            assert k % 4 == 0
            k = k // 4

            if 0 <= k <= f:
                for subset in itertools.combinations(fixed, k):
                    oneline = base[:]
                    for i in subset:
                        oneline[n + i] = True
                        oneline[n - i] = True
                    cl = Clan(oneline, cls.TYPE_B)
                    assert cl.clan_type() == 2 * p - 2 * q - 1
                    yield cl

    @classmethod
    def all_a(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls.all_a(p, n - p):
                    yield clan
        else:
            for c in cls._all_a(p, q):
                yield c

    @classmethod
    def all_b(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls.all_b(p, n - p):
                    yield clan
        for c in cls._all_b(p, q):
            yield c

    @classmethod
    def all_c1(cls, n):
        for w in SignedPermutation.involutions(n):
            fixed = [i for i in range(1, n + 1) if w(i) == i]
            base = [i if i < w(i) else w(i) if w(i) < i else True for i in range(-n, 0)]
            base += [i if i < w(i) else w(i) if w(i) < i else False for i in range(1, n + 1)]
            for p in range(len(fixed) + 1):
                for subset in itertools.combinations(fixed, p):
                    oneline = base[:]
                    for i in subset:
                        # the tuple is 0-indexed
                        oneline[n + i - 1] = True
                        oneline[n - i] = False
                    yield Clan(oneline, cls.TYPE_C1)

    @classmethod
    def _all_c2(cls, p, q):
        n = p + q
        for w in SignedPermutation.involutions(n):
            if any(w(i) == -i for i in range(1, n + 1)):
                continue
            fixed = [i for i in range(1, n + 1) if w(i) == i]
            base = [i if i < w(i) else w(i) if w(i) < i else False for i in range(-n, 0)]
            base += [i if i < w(i) else w(i) if w(i) < i else False for i in range(1, n + 1)]

            f = len(fixed)
            # 2 * k - 2 * (f - k) == 4 * k - 2 * f == 2 * p - 2 * q
            # 4 * k == 2 * p - 2 * q + 2 * f
            k = (2 * p - 2 * q + 2 * f)
            assert k % 4 == 0
            k = k // 4

            if 0 <= k <= f:
                for subset in itertools.combinations(fixed, k):
                    oneline = base[:]
                    for i in subset:
                        oneline[n + i - 1] = True
                        oneline[n - i] = True
                    cl = Clan(oneline, cls.TYPE_C2)
                    assert cl.clan_type() == 2 * p - 2 * q
                    yield cl

    @classmethod
    def all_c2(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls.all_b(p, n - p):
                    yield clan
        for c in cls._all_c2(p, q):
            yield c

    def rank(self):
        if self.family == self.TYPE_A:
            return len(self.oneline)
        if self.family == self.TYPE_B:
            return (len(self.oneline) - 1) // 2
        if self.family == self.TYPE_C1:
            return len(self.oneline) // 2
        if self.family == self.TYPE_C2:
            return len(self.oneline) // 2

    def generators(self):
        start = 1 if self.family == self.TYPE_A else 0
        for i in range(start, self.rank()):
            yield i

    def get_clan_words(self):
        if self not in CLAN_WORDS_CACHE:
            ans = []
            for i in self.generators():
                other = i * self
                if other != self:
                    ans += [w + (i,) for w in other.get_clan_words()]
            CLAN_WORDS_CACHE[self] = ans if len(ans) > 0 else [()]
        return CLAN_WORDS_CACHE[self]

    def get_atoms(self):
        def s(i):
            if self.family == self.TYPE_A:
                return Permutation.s_i(i)
            elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
                return SignedPermutation.s_i(i, self.rank())

        if self.family == self.TYPE_A:
            identity = Permutation()
        elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            identity = SignedPermutation.identity(self.rank())

        if self not in CLAN_ATOMS_CACHE:
            ans = set()
            for i in self.generators():
                other = i * self
                if other != self:
                    for w in other.get_atoms():
                        ans.add(w * s(i))
            CLAN_ATOMS_CACHE[self] = ans if len(ans) > 0 else {identity}
        return CLAN_ATOMS_CACHE[self]

    def cycles(self):
        cycles = []
        for i, a in enumerate(self.oneline):
            if type(a) == int and i + 1 < a:
                cycles.append((i + 1, a))
        return cycles

    def richardson_springer_map(self):
        def phi_a(cycles):
            w = Permutation()
            for (i, j) in cycles:
                w *= Permutation.t_ij(i, j)
            return w

        def phi_bcd(cycles, n):
            w = SignedPermutation.identity(n)
            for (i, j) in cycles:
                if abs(i) <= j:
                    if i < 0:
                        s = SignedPermutation.reflection_s(-i, j, n)
                    else:
                        s = SignedPermutation.reflection_t(i, j, n)
                    w *= s
            return w

        if self.family == self.TYPE_A:
            return phi_a(self.cycles())

        if self.family == self.TYPE_B:
            n = self.rank()
            cycles = [(i - n - 1, j - n - 1) for (i, j) in self.cycles()]
            return phi_bcd(cycles, n)

        if self.family in [self.TYPE_C1, self.TYPE_C2]:
            n = self.rank()
            cycles = [(
                i - n - 1 if i <= n else i - n,
                j - n - 1 if j <= n else j - n) for (i, j) in self.cycles()]
            return phi_bcd(cycles, n)

    def __call__(self, i):
        return self.oneline[i - 1]

    def _conjugate(self, i, j=None):
        if j is None:
            j = i + 1
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            newline = list(self.oneline)
            a, b = self(i), self(j)
            if type(a) == int:
                newline[a - 1] = j
            if type(b) == int:
                newline[b - 1] = i
            newline[i - 1], newline[j - 1] = b, a
            return Clan(newline, self.family)

    def _translate(self, i):
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            assert type(self(i)) != int and type(self(i + 1)) != int
            newline = list(self.oneline)
            newline[i - 1], newline[i] = i + 1, i
            return Clan(newline, self.family)

    def _multiply_a(self, i):
        assert 1 <= i < self.rank()
        a, b = self(i), self(i + 1)
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(i)
        if type(a) == int and type(b) != int and a < i:
            return self._conjugate(i)
        if type(a) != int and type(b) == int and i + 1 < b:
            return self.conjugate(i)
        if type(a) != int and type(b) != int and a != b:
            return self._translate(i)
        return self

    def _multiply_b(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 2)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 2)
            if type(a) != int and type(self(n + 1)) != int and a != self(n + 1):
                oneline = list(self.oneline)
                oneline[n - 1], oneline[n + 1] = n + 2, n
                oneline[n] = not oneline[n]
                return Clan(oneline, self.TYPE_B)
            return self

        a, b = self(n + 1 + i), self(n + 2 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self._conjugate(n + 1 + i)
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + 1 + i)._conjugate(n - i)
        if type(a) == int and type(b) != int and a < n + 1 + i:
            return self._conjugate(n + 1 + i)._conjugate(n - i)
        if type(a) != int and type(b) == int and n + 2 + i < b:
            return self._conjugate(n + 1 + i)._conjugate(n - i)
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + 1 + i)._translate(n - i)
        return self

    def _multiply_c1(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 1)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 1)
            if type(a) != int and type(b) != int and a != b:
                oneline = list(self.oneline)
                oneline[n - 1], oneline[n] = n + 1, n
                return Clan(oneline, self.TYPE_C1)
            return self

        a, b = self(n + i), self(n + 1 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self._conjugate(n + i)
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) == int and type(b) != int and a < n + i:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) != int and type(b) == int and n + 1 + i < b:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + i)._translate(n - i)
        return self

    def _multiply_c2(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 1)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 1)
            return self

        a, b = self(n + i), self(n + 1 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) == int and type(b) != int and a < n + i:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) != int and type(b) == int and n + 1 + i < b:
            return self._conjugate(n + i)._conjugate(n - i)
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + i)._translate(n - i)
        return self

    def __rmul__(self, i):
        assert type(i) == int
        if self.family == self.TYPE_A:
            return self._multiply_a(i)
        if self.family == self.TYPE_B:
            return self._multiply_b(i)
        if self.family == self.TYPE_C1:
            return self._multiply_c1(i)
        if self.family == self.TYPE_C2:
            return self._multiply_c2(i)

