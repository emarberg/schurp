from partitions import StrictPartition, Partition
from permutations import Permutation
from vectors import Vector


involution_stanley_cache = {(): {(): 1}}
fpf_stanley_cache = {(): {(): 1}}
stanley_cache = {(): {(): 1}}


class StanleyExpander:

    @property
    def cache(self):
        return stanley_cache

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
            Schur(Partition(*k)): v
            for k, v in self.cache[self.key].items()
        })

    def expand(self):
        if self.key in self.cache:
            return self.get_cache()
        else:
            children, _ = self.get_children()
            if children:
                ans = Vector()
                for u in children:
                    ans += StanleyExpander(u).expand()
            else:
                q = Schur(self.get_grassmannian_shape())
                ans = Vector({q: 1})
            self.set_cache(ans)
            return ans

    def get_children(self):
        u, r, s = self.descend()
        if u is None:
            return [], None
        children = StanleyExpander(u).phi_minus(r)
        if len(children) == 0:
            n = len(self.w.oneline)
            c = Permutation.cycle(range(1, n + 2))
            new_w = c * self.w * c.inverse()
            return StanleyExpander(new_w).get_children()
        return children, u

    def descend(self):
        if self.is_grassmannian():
            return None, None, None
        n = len(self.w.oneline)
        for r in reversed(range(1, n + 1)):
            for s in reversed(range(r + 1, n + 1)):
                if self.w(s) < self.w(r):
                    t = Permutation.cycle([r, s])
                    y = self.w * t
                    return y, r, s

    def get_descents(self):
        n = len(self.w.oneline)
        for i in range(1, n + 1):
            if self.w(i) > self.w(i + 1):
                yield i

    def is_grassmannian(self):
        return len(list(self.get_descents())) <= 1

    def get_grassmannian_shape(self):
        des = list(self.get_descents())
        assert len(des) <= 1
        n = des[0] if des else 0
        return Partition(*reversed(sorted([self.w(i) - i for i in range(1, n + 1)])))

    def phi_plus(self, r):
        ans = []
        n = max(r, len(self.w.oneline))
        for i in range(r + 1, n + 2):
            tau = self.tau(r, i)
            if (tau.length() - self.w.length()) == 1:
                ans += [tau]
        return sorted(set(ans))

    def phi_minus(self, r):
        ans = []
        for i in range(1, r):
            tau = self.tau(i, r)
            if (tau.length() - self.w.length()) == 1:
                ans += [tau]
        return sorted(set(ans))

    def tau(self, i, j):
        t = Permutation.cycle([i, j])
        return self.w * t


class FPFStanleyExpander:

    @property
    def cache(self):
        return fpf_stanley_cache

    def __init__(self, w):
        assert type(w) == Permutation
        assert w.is_fpf_involution()
        self.w = w
        self._key = None

    def __repr__(self):
        return 'FPF Expand[ %s ]' % str(self.w)

    @property
    def key(self):
        if self._key is None:
            k = self.w.oneline
            while k and k[-2] == len(k) and k[-1] == len(k) - 1:
                k = k[:-2]
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
            children, _ = self.get_children()
            if children:
                ans = Vector()
                for u in children:
                    ans += FPFStanleyExpander(u).expand()
            else:
                q = SchurP(self.get_grassmannian_shape())
                ans = Vector({q: 1})
            self.set_cache(ans)
            return ans

    def get_children(self):
        u, r, s = self.descend()
        if u is None:
            return [], None
        # sanity check
        # assert [w] == cls.phi_plus(u, r) and u(r) <= r
        children = FPFStanleyExpander(u).phi_minus(u(r))
        if len(children) == 0:
            n = len(self.w.oneline)
            c = Permutation.cycle(range(1, n + 3)) ** 2
            new_w = c * self.w * c.inverse() * Permutation.s_i(1)
            return FPFStanleyExpander(new_w).get_children()
        return children, u

    def descend(self):
        if self.is_grassmannian():
            return None, None, None
        n = len(self.w.oneline)
        for r in reversed(range(1, n + 1)):
            for s in reversed(range(r + 1, n + 1)):
                if self.w(s) < min([self.w(r), r]):
                    t = Permutation.cycle([r, s])
                    y = t * self.w * t
                    return y, r, s

    def get_dearc_descents(self):
        v = self.dearc()
        n = len(v.oneline)
        for i in range(1, n + 1):
            if v(i) > v(i + 1) and v(i + 1) <= i:
                yield i

    def dearc(self):
        cycles = [tuple(sorted(c)) for c in self.w.cycles]
        cycles = [(p, q) for (p, q) in cycles if any((a, b) for (a, b) in cycles if p < b < q)]
        v = Permutation()
        for p, q in cycles:
            v *= Permutation.cycle((p, q))
        return v

    def is_grassmannian(self):
        return len(list(self.get_dearc_descents())) <= 1

    def get_grassmannian_shape(self):
        des = list(self.get_dearc_descents())
        assert len(des) <= 1
        n = des[0] if des else 0
        r = max(self.dearc().support()) - n if des else 0
        return StrictPartition(*[n - self.dearc()(n + 1 + i) for i in range(r)])

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

    def tau(self, i, j):
        t = Permutation.cycle([i, j])
        return t * self.w * t


class InvStanleyExpander:

    @property
    def cache(self):
        return involution_stanley_cache

    def __init__(self, w):
        assert type(w) == Permutation
        self.w = w
        self._key = None

    def __repr__(self):
        return 'Inv Expand[ %s ]' % str(self.w)

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
            children, _ = self.get_children()
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
            return [], None
        # sanity check
        # assert [w] == cls.phi_plus(u, r) and u(r) <= r
        children = InvStanleyExpander(u).phi_minus(u(r))
        if len(children) == 0:
            n = len(self.w.oneline)
            c = Permutation.cycle(range(1, n + 2))
            new_w = c * self.w * c.inverse()
            return InvStanleyExpander(new_w).get_children()
        return children, u

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
            return r * w
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
                return r * s * w * r
            elif (i, j) in [(a, b), (c, d), (a, d)] and w(a) == c and w(b) == d:
                r = Permutation.cycle([a, b])
                return r * w * r
        return w

    def tau_inverse(self, i, j):
        w = self.w
        a_tup = tuple(sorted(set([i, j, w(i), w(j)])))
        if len(a_tup) == 2 and w(i) == j:
            r = Permutation.cycle([i, j])
            return r * w
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
                return s * r * w * r
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


class Determinant:

    def __init__(self, matrix):
        self.matrix = matrix
        n = len(matrix)
        assert all(len(matrix[i]) == n for i in range(n))

    def __repr__(self):
        m = max([len(str(i)) for row in self.matrix for i in row])

        def pad(s):
            return str(s) + (m - len(str(s))) * ' '

        return '\n'.join([' '.join([pad(s) for s in row]) for row in self.matrix])

    def evaluate(self):
        ans = Vector()
        n = len(self.matrix)
        for w in Permutation.all(n):
            term = (-1)**len(w)
            for i in range(1, n + 1):
                term = term * self.matrix[i - 1][w(i) - 1]
            ans = ans + term
        return ans


class Pfaffian:

    def __init__(self, matrix):
        self.matrix = matrix
        n = len(matrix)
        assert all(len(matrix[i]) == n for i in range(n))
        assert all(matrix[i][j] == -matrix[j][i] for i in range(n) for j in range(n))

    def __repr__(self):
        m = max([len(str(i)) for row in self.matrix for i in row])

        def pad(s):
            return str(s) + (m - len(str(s))) * ' '

        return '\n'.join([' '.join([pad(s) for s in row]) for row in self.matrix])

    @classmethod
    def fpf_involutions(cls, n):
        if n % 2 != 0 or n <= 0:
            return
        m = n // 2
        indices = m * [0]
        while indices[-1] != 1:
            z = Permutation()
            numbers = list(range(1, n + 1))
            for i in indices:
                a = numbers[-1]
                b = numbers[i]
                numbers = numbers[:i] + numbers[i + 1:-1]
                z = z * Permutation.cycle([a, b])
            yield z

            indices[0] += 1
            for i in range(m - 1):
                if indices[i] == n - 1 - 2 * i:
                    indices[i] = 0
                    indices[i + 1] += 1

    def evaluate(self):
        ans = Vector()
        n = len(self.matrix)
        for w in self.fpf_involutions(n):
            term = (-1)**w.involution_length()
            for c in w.cycles:
                a, b = tuple(sorted(c))
                term = term * self.matrix[b - 1][a - 1]
            ans = ans + term
        return ans


class SchurP:

    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) in [Partition, StrictPartition]:
            args = args[0].parts
        self.mu = StrictPartition(*args)

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
        if type(other) == int:
            return Vector({self: other})
        assert type(other) == type(self)
        # return Vector({SchurP(p): v for p, v in self.mu.pieri(other.mu(1)).items()})
        u = self.mu.to_grassmannian()
        v = other.mu.to_grassmannian()
        n = len(u.oneline)
        w = Permutation(u.oneline + [i + n for i in v.oneline])
        return InvStanleyExpander(w).expand()

    def __rmul__(self, i):
        return self.__mul__(i)

    @classmethod
    def _check_vector(cls, vector):
        if type(vector) == cls:
            vector = Vector({vector: 1})
        assert type(vector) == Vector
        assert all(type(k) == cls for k in vector.keys())
        return vector

    @classmethod
    def to_s_basis(cls, vector):
        return SchurQ.to_s_basis(cls.to_q_basis(vector))

    @classmethod
    def to_q_basis(cls, vector):
        vector = cls._check_vector(vector)
        assert all(v % 2**len(k.mu) == 0 for k, v in vector.items())
        ans = Vector()
        for s, v in vector.items():
            ans += Vector({SchurQ(s.mu): v // 2**len(s.mu)})
        return ans

    @classmethod
    def to_schur_basis(cls, vector):
        vector = cls._check_vector(vector)
        ans = Vector()
        for s, val in vector.items():
            assert type(s.mu) == StrictPartition
            for w in s.mu.to_grassmannian().get_atoms():
                ans += StanleyExpander(w).expand() * val
        return ans

    def skew(self, nu, verbose=False):
        mu = self.mu
        ans = SchurQ(mu).skew(nu, verbose)
        e = len(mu) - len(nu)
        return SchurQ.to_p_basis(ans) / 2**e

s_lambda_cache = {}


class SchurQ(SchurP):

    @classmethod
    def s_lambda(cls, *args):
        if len(args) == 1 and type(args[0]) in [Partition, StrictPartition]:
            args = args[0].parts
        mu = Partition(*args)
        if len(mu) == 0:
            return Vector({SchurQ(): 1})
        if tuple(mu.parts) not in s_lambda_cache:
            n = len(mu)
            matrix = [[Vector() for i in range(n)] for j in range(n)]
            for i in range(1, n + 1):
                for j in range(1, n + 1):
                    q = mu(i) - i + j
                    if q >= 0:
                        matrix[i - 1][j - 1] = Vector({SchurQ(q): 1})
            det = Determinant(matrix)
            s_lambda_cache[tuple(mu.parts)] = det.evaluate()
        return s_lambda_cache[tuple(mu.parts)]

    @classmethod
    def to_p_basis(cls, vector):
        vector = cls._check_vector(vector)
        ans = Vector()
        for s, v in vector.items():
            ans += Vector({SchurP(s.mu): 2**len(s.mu) * v})
        return ans

    @classmethod
    def to_s_basis(cls, vector):
        return cls.decompose_s_lambda(vector)

    @classmethod
    def to_schur_basis(cls, vector):
        return SchurP.to_schur_basis(cls.to_p_basis(vector))

    @classmethod
    def decompose_s_lambda(cls, vector):
        vector = cls._check_vector(vector)
        if vector.is_zero():
            return Vector()
        else:
            q = min(vector)
            c = vector[q]
            s = cls.s_lambda(q.mu) * c
            return Vector({SchurS(q.mu): c}) + cls.decompose_s_lambda(vector - s)

    def skew(self, nu, verbose=False):
        assert type(nu) == StrictPartition
        assert self.mu.contains(nu)

        m = len(self.mu)
        n = len(nu)

        if m == n == 0:
            return Vector({SchurQ(StrictPartition()): 1})

        if (m + n) % 2 != 0:
            n += 1

        matrix = [[Vector() for i in range(m + n)] for j in range(m + n)]
        for i in range(m):
            for j in range(i + 1, m):
                x = StrictPartition(self.mu(i + 1), self.mu(j + 1))
                matrix[i][j] = Vector({SchurQ(x): 1})
                matrix[j][i] = Vector({SchurQ(x): -1})
        for i in range(m):
            for j in range(n):
                q = self.mu(i + 1) - nu(n - j)
                if q >= 0:
                    x = StrictPartition(q)
                    matrix[i][m + j] = Vector({SchurQ(x): 1})
                    matrix[m + j][i] = Vector({SchurQ(x): -1})

        if verbose:
            lengths = [len(str(matrix[i][j])) for i in range(m + n) for j in range(m + n)]

            def pad(t):
                return (max(lengths) - len(str(t))) * ' ' + str(t)

            print('')
            for row in matrix:
                print(' '.join(pad(t) for t in row))

        return Pfaffian(matrix).evaluate()

    def __repr__(self):
        return 'Q(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __mul__(self, other):
        if type(other) == int:
            return Vector({self: other})
        assert type(other) == type(self)
        mult = len(self.mu) + len(other.mu)
        u = self.mu.to_grassmannian()
        v = other.mu.to_grassmannian()
        assert type(self.mu) == StrictPartition
        assert type(other.mu) == StrictPartition
        n = len(u.oneline)
        w = Permutation(u.oneline + [i + n for i in v.oneline])
        vector = InvStanleyExpander(w).expand()
        return Vector({
            SchurQ(p.mu): 2**(mult - len(p.mu)) * v
            for p, v in vector.items()
        })


class SchurS(SchurP):

    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) in [Partition, StrictPartition]:
            args = args[0].parts
        self.mu = Partition(*args)

    def __repr__(self):
        return 'S(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __mul__(self, other):
        if type(other) == int:
            return Vector({self: other})
        assert type(other) == type(self)
        return SchurQ.to_s_basis(SchurS.to_q_basis(self) * SchurS.to_q_basis(other))

    def __add__(self, other):
        if type(other) == Vector:
            return Vector.base(self) + other
        else:
            assert type(other) == type(self)
            return SchurQ.to_s_basis(self.to_q_basis(Vector.base(self) + Vector.base(other)))

    def __radd__(self, other):
        return self.__add__(other)

    @classmethod
    def to_p_basis(cls, vector):
        return SchurQ.to_p_basis(cls.to_q_basis(vector))

    @classmethod
    def to_q_basis(cls, vector):
        vector = cls._check_vector(vector)
        ans = Vector()
        for s, v in vector.items():
            for x, u in SchurQ.s_lambda(s.mu).items():
                ans += Vector({x: u * v})
        return ans

    @classmethod
    def to_schur_basis(cls, vector):
        return SchurP.to_schur_basis(cls.to_p_basis(vector))


class Schur(SchurP):

    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) in [Partition, StrictPartition]:
            args = args[0].parts
        self.mu = Partition(*args)

    def __repr__(self):
        return 'schur(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __mul__(self, other):
        if type(other) == int:
            return Vector({self: other})
        assert type(other) == type(self)
        u = self.mu.to_grassmannian()
        v = other.mu.to_grassmannian()
        assert type(self.mu) == Partition
        assert type(other.mu) == Partition
        n = len(u.oneline)
        w = Permutation(u.oneline + [i + n for i in v.oneline])
        return StanleyExpander(w).expand()

    def __add__(self, other):
        if type(other) == Vector:
            return Vector.base(self) + other
        else:
            assert type(other) == type(self)
            return Vector.base(self) + other.to_schur_basis(other)

    def __radd__(self, other):
        return self.__add__(other)
