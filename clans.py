from signed import SignedPermutation
from permutations import Permutation
import itertools
import subprocess


BASE_DIRECTORY = '/Users/emarberg/examples/clans/'


CLAN_WORDS_CACHE = {}
CLAN_ATOMS_CACHE = {}
CLAN_HECKE_ATOMS_CACHE = {}
CLAN_HECKE_ATOMS_EXTENDED_CACHE = {}


class Clan:

    TYPE_A = 'clans for G = GL(p+q), K = GL(p) x GL(q)'
    TYPE_B = 'clans for G = SO(2n+1), K = S(O(2p) x O(2q+1))'
    TYPE_C1 = 'clans for G = Sp(2n), K = GL(n)'
    TYPE_C2 = 'clans for G = Sp(2n), K = Sp(2p) x Sp(2q)'
    TYPE_D1 = 'clans for G = SO(2n), K = S(O(2p) x O(2q))'
    TYPE_D2 = 'clans for G = SO(2n), K = S(O(2p-1) x O(2q-1)) with p+q = n+1'
    TYPE_D3 = 'clans for G = SO(2n), K = GL(n)'

    def __init__(self, oneline, family=TYPE_A):
        s = sorted(oneline)
        assert not any(type(s[i]) == int and type(s[i + 1]) == int and type(s[i + 2]) == int and s[i] == s[i + 2] for i in range(len(s) - 2))

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
        assert all(type(oneline[i]) == bool or oneline[oneline[i] - 1] - 1 == i for i in range(len(oneline)))

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

    @classmethod
    def create_d1(cls, oneline):
        return Clan(oneline, cls.TYPE_D1)

    @classmethod
    def create_d2(cls, oneline):
        return Clan(oneline, cls.TYPE_D2)

    @classmethod
    def create_d3(cls, oneline):
        return Clan(oneline, cls.TYPE_D3)

    @classmethod
    def _draw(cls, clans, folder, filename):
        assert len(clans) > 0

        def printer(c):
            return str(c) + '\n' + c.richardson_springer_involution().cycle_repr()

        edges = []
        for c in clans:
            for i in c.generators():
                d, doubled = c.weak_order_action(i)
                if doubled is not None:
                    edges.append((c, d, i, doubled))
        
        s = []
        s += ['digraph G {']
        s += ['    rankdir=BT;']
        s += ['    overlap=false;']
        s += ['    splines=true;']
        
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        for x in clans:
            s += ['    "%s"' % printer(x) + (' [fillcolor=white];' if not x.is_alternating() else ';')]
        for c, d, i, doubled in edges:
            s += ['    "%s" -> "%s" [label="%s",color="%s"];' % (printer(c), printer(d), str(i), 'blue' if doubled else 'black')]
        s += ['}']
        s = '\n'.join(s)
        
        dot_filename = BASE_DIRECTORY + 'dot/' + folder + '/%s.dot' % filename
        png_filename = BASE_DIRECTORY + 'png/' + folder + '/%s.png' % filename
        with open(dot_filename, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])

    @classmethod
    def draw_a(cls, p, q=None):
        clans = list(cls.all_a(p, q))
        folder = 'A'
        filename = folder + '_'+ ('n' + str(p) if q is None else 'p' + str(p) + 'q' + str(q)) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_b(cls, p, q=None):
        clans = list(cls.all_b(p, q))
        folder = 'B'
        filename = folder + '_'+ ('n' + str(p) if q is None else 'p' + str(p) + 'q' + str(q)) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_c1(cls, n):
        clans = list(cls.all_c1(n))
        folder = 'C1'
        filename = folder + '_'+ 'n' + str(n) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_c2(cls, p, q=None):
        clans = list(cls.all_c2(p, q))
        folder = 'C2'
        filename = folder + '_'+ ('n' + str(p) if q is None else 'p' + str(p) + 'q' + str(q)) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_d1(cls, p, q=None):
        clans = list(cls.all_d1(p, q))
        folder = 'D1'
        filename = folder + '_'+ ('n' + str(p) if q is None else 'p' + str(p) + 'q' + str(q)) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_d2(cls, p, q=None):
        clans = list(cls.all_d2(p, q))
        folder = 'D2'
        filename = folder + '_'+ ('n' + str(p) if q is None else 'p' + str(p) + 'q' + str(q)) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    @classmethod
    def draw_d3(cls, n):
        clans = list(cls.all_d3(n))
        folder = 'D3'
        filename = folder + '_'+ 'n' + str(n) + '_' + str(len(clans))
        cls._draw(clans, folder, filename)

    def __repr__(self):
        l = [str(min(i, self(i))) if type(i) == int else '+' if i else '-' for i in self.oneline]
        return ' '.join(l)

    def __eq__(self, other):
        assert type(other) == Clan and self.family == other.family
        return self.oneline == other.oneline

    def __hash__(self):
        return hash((self.family, self.oneline))

    def is_alternating(self):
        s = [a for a in self.oneline if type(a) is bool]
        if self.family == self.TYPE_A:
            return all(s[i] == s[i + 1] for i in range(len(s) - 1)) or not any(s[i] == s[i + 1] for i in range(len(s) - 1))
        elif self.family in [self.TYPE_B]:
            return all(s[i] == s[i + 1] for i in range(len(s)//2 - 1)) or not any(s[i] == s[i + 1] for i in range(len(s) - 1))
        elif self.family == self.TYPE_C1:
            return not any(s[i] == s[i + 1] for i in range(len(s) - 1))
        elif self.family in [self.TYPE_C2, self.TYPE_D1, self.TYPE_D2]:
            return all(s[i] == s[i + 1] for i in range(len(s)//2 - 1)) or not any(s[i] == s[i + 1] for i in range(len(s)//2 - 1))
        elif self.family == self.TYPE_D3:
            return len([i for i in self.oneline if type(i) == bool]) <= 2 * (2 if self.rank() % 2 == 0 else 1)
        else:
            raise Exception

    def is_strongly_alternating(self):
        s = [a for a in self.oneline if type(a) is bool]
        if self.family in [self.TYPE_A, self.TYPE_C1]:
            return not any(s[i] == s[i + 1] for i in range(len(s) - 1))
        elif self.family in [self.TYPE_B]:
            return not any(s[i] == s[i + 1] for i in range(len(s) - 1))
        elif self.family in [self.TYPE_C2, self.TYPE_D1, self.TYPE_D2]:
            return not any(s[i] == s[i + 1] for i in range(len(s)//2 - 1))
        elif self.family == self.TYPE_D3:
            return len([i for i in self.oneline if type(i) == bool]) <= 1
        else:
            raise Exception

    def is_signless(self):
        return not any(type(i) == bool for i in self.oneline)

    def is_matchless(self):
        return not any(type(i) == int for i in self.oneline)

    def plus(self):
        return len([x for x in self.oneline if x is True])

    def minus(self):
        return len([x for x in self.oneline if x is False])

    def clan_type(self):
        return self.plus() - self.minus()

    def is_aligned(self, matching, verbose=False):
        for (i, j) in matching:
            if self.family in [self.TYPE_B, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2] and i + j == 0:
                continue
            if self.family == self.TYPE_A:
                if i < 0:
                    assert j < 0 or i + j == 0
                    continue
                pair = (self(i), self(j))
            elif self.family == self.TYPE_B:
                i, j = self.rank() + 1 + i, self.rank() + 1 + j
                pair = (self(i), self(j))
            elif self.family in [self.TYPE_C1, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
                i = self.rank() + i + (0 if i > 0 else 1)
                j = self.rank() + j + (0 if j > 0 else 1)
                pair = (self(i), self(j))
            else:
                raise Exception
            assert type(pair[0]) == bool and type(pair[1]) == bool
            if pair in [(True, True), (False, False)]:
                return False

        if self.family == self.TYPE_A:
            support = [a for pair in matching for a in pair if min(pair) > 0]
            trivial = [a for a in range(1, self.rank() + 1) if type(self(a)) == bool and a not in support]
            signs = [self(i) for i in trivial]
            if any(signs[i] != signs[i + 1] for i in range(len(signs) - 1)):
                return False
            return True

        trivial = [a for (a, b) in sorted(matching) if -a == b]
        k = abs(self.clan_type()) // 2

        if self.family == self.TYPE_C2:
            return len(trivial) == k
        if self.family == self.TYPE_B:
            if len(trivial) < k:
                return False
            signs = [self(self.rank() + 1 + i) for i in trivial]
            signs += [self(self.rank() + 1)]
            if any(signs[i] == signs[i + 1] for i in range(k, len(signs) - 1)):
                return False
        if self.family in [self.TYPE_D1, self.TYPE_D2]:
            if len(trivial) < k:
                return False
            signs = [self(self.rank() + 1 + i) for i in reversed(trivial)]
            if any(signs[i] != signs[i + 1] for i in range(0, len(signs) - 1)):
                return False
        if self.family == self.TYPE_D3:
            if len(trivial) % 2 != self.rank() % 2:
                return False
            signs = [self(self.rank() + 1 + i) for i in reversed(trivial)]
            
            # inner-most 2 trivial blocks have same sign,
            # then next inner-most 2 trivial blocks have same sign,
            # and so on

            if any(signs[i] != signs[i + 1] for i in range(0, len(signs) - 1, 2)):
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
    def symmetric_clans(cls, p, q, disallow_negations=False):
        n = (p + q) // 2
        mod = 1 if (p + q) % 2 == 0 else 0
        for w in SignedPermutation.involutions(n):
            if disallow_negations and any(w(i) == -i for i in range(1, n + 1)):
                continue

            fixed = [i for i in range(1, n + 1) if w(i) == i]
            f = len(fixed)

            # suppose p + q is even. let k be number of positive +.
            # then p - q = 2*k - 2*(f-k) = 4*k - 2*f
            # so k = (p - q + 2*f) / 4
            # 
            # instead suppose p + q is odd.
            # let e = 1 if central + and e = -1 if central -
            # if k is number of positive + then must have
            #
            #   p - q = 2*k + e - 2*(f-k) = 4*k - 2*f + e
            #
            # thus e = 1 iff p - q + 2 * f is 1 mod 4
            # and k = (p - q + 2*f - e) / 4

            base = [i if i < w(i) else w(i) if w(i) < i else False for i in range(-n, 0)]
            
            if (p + q) % 2 != 0:
                base += [(p - q + 2 * f) % 4 == 1]
                e = 1 if base[n] else -1
                k = p - q + 2 * f - e
                assert k % 4 == 0
                k = k // 4
            else:
                k = p - q + 2 * f
                if k % 4 != 0:
                    continue
                k = k // 4
            
            base += [i if i < w(i) else w(i) if w(i) < i else False for i in range(1, n + 1)]

            if 0 <= k <= f:
                for subset in itertools.combinations(fixed, k):
                    oneline = base[:]
                    for i in subset:
                        oneline[n - mod + i] = True
                        oneline[n - i] = True
                    yield oneline

    @classmethod
    def all_a(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls._all_a(p, n - p):
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
        else:
            for oneline in cls.symmetric_clans(2 * p, 2 * q + 1):
                cl = Clan(oneline, cls.TYPE_B)
                assert cl.clan_type() == 2 * p - 2 * q - 1
                yield cl

    @classmethod
    def skew_symmetric_clans(cls, p, q, disallow_negations=False):
        assert p == q
        n = p
        for w in SignedPermutation.involutions(n):
            if disallow_negations and any(w(i) == -i for i in range(1, n + 1)):
                continue

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
                    yield oneline

    @classmethod
    def all_c1(cls, n):
        for oneline in cls.skew_symmetric_clans(n, n):
            yield Clan(oneline, cls.TYPE_C1)

    @classmethod
    def all_c2(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls.all_c2(p, n - p):
                    yield clan
        else:
            for oneline in cls.symmetric_clans(2 * p, 2 * q, disallow_negations=True):
                yield Clan(oneline, cls.TYPE_C2)

    @classmethod
    def all_d1(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n):
                for clan in cls.all_d1(p, n - p):
                    yield clan
        else:
            for oneline in cls.symmetric_clans(2 * p, 2 * q):
                cl = Clan(oneline, cls.TYPE_D1)
                assert cl.clan_type() == 2 * p - 2 * q
                yield cl

    @classmethod
    def all_d2(cls, p, q=None):
        if q is None:
            n = p
            for p in range(1, n + 1):
                q = n + 1 - p
                for clan in cls.all_d2(p, q):
                    yield clan
        else:
            for oneline in cls.symmetric_clans(2 * p - 1, 2 * q - 1):
                cl = Clan(oneline, cls.TYPE_D2)
                assert cl.clan_type() == (2 * p - 1) - (2 * q - 1)
                yield cl

    @classmethod
    def all_d3(cls, n):
        for oneline in cls.skew_symmetric_clans(n, n, disallow_negations=True):
            a = len([i for i in range(n) if oneline[i] is False])
            b = len([
                (i, j)
                for i in range(n)
                for j in range(i + 1, n)
                if type(oneline[i]) == int and oneline[i] == oneline[j]
            ])
            if (a + b) % 2 == 0:
                yield Clan(oneline, cls.TYPE_D3)

    def rank(self):
        if self.family == self.TYPE_A:
            return len(self.oneline)
        else:
            return len(self.oneline) // 2

    def generators(self):
        if self.rank() >= 2 and self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            yield -1
        elif self.rank() >= 1 and self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            yield 0
        for i in range(1, self.rank()):
            yield i

    def simple_generator(self, i):
        if self.family == self.TYPE_A:
            return Permutation.s_i(i)
        elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            return SignedPermutation.s_i(i, self.rank())
        elif self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return SignedPermutation.ds_i(i, self.rank())
        else:
            raise Exception

    def weyl_group(self):
        if self.family == self.TYPE_A:
            return Permutation.all(self.rank())
        elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            return SignedPermutation.all(self.rank())
        elif self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return SignedPermutation.all(self.rank(), dtype=True)
        else:
            raise Exception

    def weyl_group_identity(self):
        if self.family == self.TYPE_A:
            return Permutation()
        elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return SignedPermutation.identity(self.rank())
        else:
            raise Exception

    def weyl_group_longest_element(self):
        n = self.rank()
        if self.family == self.TYPE_A:
            return Permutation.longest_element(n)
        elif self.family in [self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            return SignedPermutation.longest_element(n)
        elif self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return SignedPermutation.dtype_longest_element(n)
        else:
            raise Exception

    def weyl_group_bruhat_leq(self, a, b):
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            return a.strong_bruhat_less_equal(b)
        elif self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return a.dbruhat_less_equal(b)
        else:
            raise Exception

    def weyl_group_length(self, w):
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2]:
            return w.length()
        elif self.family in [self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            return w.dlength()
        else:
            raise Exception

    def weyl_group_shape(self, w):
        n = self.rank()
        k = abs(self.clan_type()) // 2

        if self.family == self.TYPE_A:
            k = abs(self.clan_type())
            sh = w.twisted_shape(n, k)
        elif self.family == self.TYPE_B:
            g = w.bbase_atom(n, k)
            sh = (g * w).shape()
            sh = tuple(sorted(sh))
        elif self.family == self.TYPE_C1:
            sh = w.shape()
        elif self.family == self.TYPE_C2:
            sh = w.fpf_shape(offset=k)
        elif self.family in [self.TYPE_D1, self.TYPE_D2]:
            g = w.dbase_atom(n, k % 2 != 0, k)
            sh = (g * w).dshape(k)
        elif self.family == self.TYPE_D3:
            sh = w.fpf_dshape()
        else:
            raise Exception

        return tuple(sorted(sh))

    def weyl_group_weight(self, w):
        top = True
        for i in self.generators():
            c, doubled = self.weak_order_action(i)
            if self != c:
                top = False
                v = w * self.simple_generator(i)
                if self.weyl_group_length(v) < self.weyl_group_length(w):
                    return c.weyl_group_weight(v) + (1 if doubled else 0)
        if top and self.weyl_group_length(w) == 0:
            return 0
        else:
            return None

    def get_clan_words(self):
        if self not in CLAN_WORDS_CACHE:
            ans = []
            for i in self.generators():
                other = i * self
                if other != self:
                    ans += [w + (i,) for w in other.get_clan_words()]
            CLAN_WORDS_CACHE[self] = ans if len(ans) > 0 else [()]
        return CLAN_WORDS_CACHE[self]

    def get_hecke_atoms(self):
        if self not in CLAN_HECKE_ATOMS_CACHE:
            atoms = self.get_atoms()
            upper_poset = {}
            for w in self.weyl_group():
                for a in atoms:
                    if self.weyl_group_bruhat_leq(a, w):
                        upper_poset[w] = set()
            for x in upper_poset:
                for y in upper_poset:
                    if self. weyl_group_bruhat_leq(x, y):
                        upper_poset[y].add(x)

            ans = {}
            def mobius(x):
                if x not in ans:
                    mu = 0
                    for z in upper_poset[x] - {x}:
                        mu += mobius(z)
                    ans[x] = 1 - mu
                return ans[x]
            for x in upper_poset:
                mobius(x)

            a = next(iter(atoms))
            CLAN_HECKE_ATOMS_CACHE[self] = set()
            for w in ans:
                if ans[w] == 0:
                    continue
                assert ans[w] == (-1)**(self.weyl_group_length(w) - self.weyl_group_length(a))
                CLAN_HECKE_ATOMS_CACHE[self].add(w)
        return CLAN_HECKE_ATOMS_CACHE[self]

    def get_atoms(self):
        if self not in CLAN_ATOMS_CACHE:
            ans = set()
            for i in self.generators():
                other = i * self
                if other != self:
                    for w in other.get_atoms():
                        ans.add(w * self.simple_generator(i))
            CLAN_ATOMS_CACHE[self] = ans if len(ans) > 0 else {self.weyl_group_identity()}
        return CLAN_ATOMS_CACHE[self]

    @classmethod
    def get_pseudo_hecke_atoms(cls, minimal_element, simple, length, conjugate, translate=None):
        m = minimal_element
        one = m * m.inverse()
        action = {one: m}
        
        level = {one}
        while level:
            newlevel = set()
            for w in level:
                z = action[w]
                for s in simple:
                    ws = w * s
                    if length(ws) > length(w):
                        newlevel.add(ws)
                        szs = None if z is None else conjugate(z, s)
                        if szs is None or (translate is None and length(szs) == length(z)):
                            action[ws] = None
                        elif translate is not None and length(szs) == length(z):
                            zs = translate(z, s)
                            action[ws] = None if zs is None else z if length(zs) < length(z) else zs
                        else:
                            action[ws] = z if length(szs) < length(z) else szs
            level = newlevel

        ans = {}
        for (w, z) in action.items():
            if z not in ans:
                ans[z] = set()
            ans[z].add(w)
        return ans

    def get_hecke_atoms_extended(self):
        z = self.richardson_springer_involution()
        v = self.richardson_springer_base()
        
        key = (self.family, self.rank(), v)
        if key not in CLAN_HECKE_ATOMS_EXTENDED_CACHE:
            n = self.rank()
            k = abs(self.clan_type()) // 2

            simple = [self.simple_generator(i) for i in self.generators()]
            length = lambda w: w.length()
            translate = lambda x,s: x * s

            t = SignedPermutation.identity(n)
            
            if self.family == self.TYPE_A:
                t = Permutation.longest_element(n)
            elif self.family in [self.TYPE_B, self.TYPE_C1]:
                pass
            elif self.family == self.TYPE_C2:
                translate = lambda x,s: x * s if (x * s).is_fpf_involution() else None
            elif self.family in [self.TYPE_D1, self.TYPE_D2]:
                length = lambda x: x.dlength()
                if k % 2 != 0:
                    t = SignedPermutation.s_i(0, n)
            elif self.family == self.TYPE_D3:
                length = lambda x: x.dlength()
                if n % 2 != 0:
                    t = SignedPermutation.s_i(0, n)
                translate = lambda x,s: x * s #x * s # if (t * x*s).is_fpf_involution() else None
            else:
                raise Exception

            conjugate = lambda x,s: t * s * t * x * s
            CLAN_HECKE_ATOMS_EXTENDED_CACHE[key] = self.get_pseudo_hecke_atoms(v, simple, length, conjugate, translate)

        dictionary = CLAN_HECKE_ATOMS_EXTENDED_CACHE[key]
        return dictionary[z]
            
    def get_atoms_extended(self):
        z = self.richardson_springer_involution()
        n = self.rank()
        offset = abs(self.clan_type()) // 2
            
        if self.family == self.TYPE_A:
            offset = abs(self.clan_type())
            return z.get_twisted_atoms(n, offset)

        elif self.family == self.TYPE_B:
            return z.get_atoms(offset)

        elif self.family == self.TYPE_C1:
            return z.get_atoms()

        elif self.family == self.TYPE_C2:
            return z.get_fpf_atoms(offset)

        elif self.family == self.TYPE_D1:
            twisted = offset % 2 != 0
            return z.get_atoms_d(twisted, offset)

        elif self.family == self.TYPE_D2:
            twisted = offset % 2 != 0
            return z.get_atoms_d(twisted, offset)

        elif self.family == self.TYPE_D3:
            return z.get_fpf_atoms_d()

        else:
            raise Exception

    def cycles(self):
        cycles = []
        for i, a in enumerate(self.oneline):
            if type(a) == int and i + 1 < a:
                cycles.append((i + 1, a))
        return cycles

    def richardson_springer_involution(self):
        w0 = self.weyl_group_longest_element()
        phi = self.richardson_springer_map()
        return w0 * phi.inverse()

    def richardson_springer_base(self):
        n = self.rank()
        k = abs(self.clan_type()) // 2

        if self.family == self.TYPE_A:
            k = (n - abs(self.clan_type())) // 2
            w = Permutation()
            for i in range(1, k + 1):
                w *= Permutation.t_ij(i, n + 1 - i)
            return self.weyl_group_longest_element() * w

        elif self.family == self.TYPE_B:
            return SignedPermutation.longest_element(n, k)

        elif self.family == self.TYPE_C1:
            return self.weyl_group_identity()

        elif self.family == self.TYPE_C2:
            w = SignedPermutation.longest_element(n, k)
            for i in range(k + 1, n, 2):
                w *= SignedPermutation.s_i(i, n)
            return w

        elif self.family in [self.TYPE_D1, self.TYPE_D2]:
            return SignedPermutation.dbase(n, k)

        elif self.family == self.TYPE_D3:
            return SignedPermutation.one_fpf_d(n)

        else:
            raise Exception

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
        elif self.family == self.TYPE_B:
            n = self.rank()
            cycles = [(i - n - 1, j - n - 1) for (i, j) in self.cycles()]
            return phi_bcd(cycles, n)
        if self.family in [self.TYPE_C1, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            n = self.rank()
            cycles = [(
                i - n - 1 if i <= n else i - n,
                j - n - 1 if j <= n else j - n) for (i, j) in self.cycles()]
            w = phi_bcd(cycles, n)
            if self.family == self.TYPE_D2:
                return w * SignedPermutation.s_i(0, n)
            else:
                return w
        else:
            raise Exception

    def __call__(self, i):
        return self.oneline[i - 1]

    def _conjugate(self, i, j=None):
        if j is None:
            j = i + 1
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            newline = list(self.oneline)
            a, b = self(i), self(j)
            if type(a) == int:
                newline[a - 1] = j
            if type(b) == int:
                newline[b - 1] = i
            newline[i - 1], newline[j - 1] = b, a
            return Clan(newline, self.family)
        else:
            raise Exception

    def _translate(self, i):
        if self.family in [self.TYPE_A, self.TYPE_B, self.TYPE_C1, self.TYPE_C2, self.TYPE_D1, self.TYPE_D2, self.TYPE_D3]:
            assert type(self(i)) != int and type(self(i + 1)) != int
            newline = list(self.oneline)
            newline[i - 1], newline[i] = i + 1, i
            return Clan(newline, self.family)
        else:
            raise Exception

    def _flip(self):
        assert len(self.oneline) % 2 == 0
        k = len(self.oneline) // 2
        if self(k) == k + 1:
            return self
        else:
            return self._conjugate(k, k + 1)

    def _multiply_a(self, i):
        assert 1 <= i < self.rank()
        a, b = self(i), self(i + 1)
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(i), False
        if type(a) == int and type(b) != int and a < i:
            return self._conjugate(i), False
        if type(a) != int and type(b) == int and i + 1 < b:
            return self._conjugate(i), False
        if type(a) != int and type(b) != int and a != b:
            return self._translate(i), False
        return self, None

    def _multiply_b(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 2)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 2), False
            if type(a) != int and type(self(n + 1)) != int and a != self(n + 1):
                oneline = list(self.oneline)
                oneline[n - 1], oneline[n + 1] = n + 2, n
                oneline[n] = not oneline[n]
                return Clan(oneline, self.TYPE_B), True
            return self, None

        a, b = self(n + 1 + i), self(n + 2 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self._conjugate(n + 1 + i), True
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + 1 + i)._conjugate(n - i), False
        if type(a) == int and type(b) != int and a < n + 1 + i:
            return self._conjugate(n + 1 + i)._conjugate(n - i), False
        if type(a) != int and type(b) == int and n + 2 + i < b:
            return self._conjugate(n + 1 + i)._conjugate(n - i), False
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + 1 + i)._translate(n - i), False
        return self, None

    def _multiply_c1(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 1)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 1), False
            if type(a) != int and type(b) != int and a != b:
                oneline = list(self.oneline)
                oneline[n - 1], oneline[n] = n + 1, n
                return Clan(oneline, self.TYPE_C1), False
            return self, None

        a, b = self(n + i), self(n + 1 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self._conjugate(n + i), True
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) == int and type(b) != int and a < n + i:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) == int and n + 1 + i < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + i)._translate(n - i), False
        return self, None

    def _multiply_c2(self, i):
        assert 0 <= i < self.rank()
        n = self.rank()
        if i == 0:
            a, b = self(n), self(n + 1)
            if type(a) == int and type(b) == int and a < b:
                return self._conjugate(n, n + 1), False
            return self, None

        a, b = self(n + i), self(n + 1 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self, None
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) == int and type(b) != int and a < n + i:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) == int and n + 1 + i < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + i)._translate(n - i), False
        return self, None

    def _multiply_d1(self, i):
        assert i in set(self.generators())
        n = self.rank()

        if i == -1:
            active = False
            a, b, c, d = self(n - 1), self(n), self(n + 1), self(n + 2)

            if type(a) == bool and type(b) == bool and a != b:
                return self._translate(n - 1)._translate(n + 1)._conjugate(n), False
            elif type(a) == int and type(b) == int and a == n:
                return self._conjugate(n - 1, n + 1), True
            elif type(a) != int and type(b) == int:
                active = (b == n + 1) or (b < n - 1)
            elif type(a) == int and type(b) != int:
                active = a < n - 1
            elif type(a) == int and type(b) == int:
                if b == n + 1 and a < n - 1:
                    active = True
                elif a == n + 2 and b < n - 1:
                    active = True
                elif a < n - 1 and b < n - 1:
                    active = True
                elif a < c < n - 1:
                    active = True 
                elif n + 2 < a < c:
                    active = True 

            if active:
                return self._conjugate(n - 1, n + 1)._conjugate(n, n + 2), False
            else:
                return self, None

        a, b = self(n + i), self(n + 1 + i)
        if type(a) == int and type(b) == int and a == n - i and b == n - i + 1:
            return self._conjugate(n + i), True
        if type(a) == int and type(b) == int and a < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) == int and type(b) != int and a < n + i:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) == int and n + 1 + i < b:
            return self._conjugate(n + i)._conjugate(n - i), False
        if type(a) != int and type(b) != int and a != b:
            return self._translate(n + i)._translate(n - i), False
        return self, None

    def _multiply_d2(self, i):
        return self._multiply_d1(i)

    def _multiply_d3(self, i):
        assert i in set(self.generators())
        if i == -1:
            ans, doubled = self._flip()._multiply_c2(1)
            return ans._flip(), doubled
        else:
            return self._multiply_c2(i)

    def weak_order_action(self, i):
        assert type(i) == int
        if self.family == self.TYPE_A:
            return self._multiply_a(i)
        elif self.family == self.TYPE_B:
            return self._multiply_b(i)
        elif self.family == self.TYPE_C1:
            return self._multiply_c1(i)
        elif self.family == self.TYPE_C2:
            return self._multiply_c2(i)
        elif self.family == self.TYPE_D1:
            return self._multiply_d1(i)
        elif self.family == self.TYPE_D2:
            return self._multiply_d2(i)
        elif self.family == self.TYPE_D3:
            return self._multiply_d3(i)
        else:
            raise Exception


    def __rmul__(self, i):
        assert type(i) == int
        if self.family == self.TYPE_A:
            return self._multiply_a(i)[0]
        elif self.family == self.TYPE_B:
            return self._multiply_b(i)[0]
        elif self.family == self.TYPE_C1:
            return self._multiply_c1(i)[0]
        elif self.family == self.TYPE_C2:
            return self._multiply_c2(i)[0]
        elif self.family == self.TYPE_D1:
            return self._multiply_d1(i)[0]
        elif self.family == self.TYPE_D2:
            return self._multiply_d2(i)[0]
        elif self.family == self.TYPE_D3:
            return self._multiply_d3(i)[0]
        else:
            raise Exception

