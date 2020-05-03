from permutations import Permutation
from signed import SignedPermutation
from polynomials import q
from vectors import Vector
import math
import time


class QPModuleElement:

    def __init__(self, qpmodule, n):
        self.qpmodule = qpmodule
        self.n = n

    def __eq__(self, other):
        assert type(self) == type(other)
        assert self.qpmodule == other.qpmodule
        return self.n == other.n

    def __hash__(self):
        return hash((self.n, self.qpmodule))

    def __repr__(self):
        return self.qpmodule.printer(self.qpmodule.reduced_word(self.n))

    def __mul__(self, i):
        result = self.qpmodule.operate(self.n, i)
        if result == self.n:
            if i in self.qpmodule.weak_ascents(self.n):
                return Vector({self: q(1)})
            else:
                return Vector({self: -q(-1)})
        new = QPModuleElement(self.qpmodule, result)
        if result < self.n:
            return Vector({self: q(1) - q(-1), new: 1})
        else:
            return Vector({new: 1})


class QPModule:

    HECKE_A = 'HECKE_A'

    @classmethod
    def GELFAND_A(cls, k): # noqa
        return 'GELFAND_A(%s)' % k

    @classmethod
    def GELFAND_BC(cls, k): # noqa
        return 'GELFAND_BC(%s)' % k

    @classmethod
    def GELFAND_D(cls, k): # noqa
        return 'GELFAND_D(%s)' % k

    def __eq__(self, other):
        return self.label == other.label

    def __hash__(self):
        return hash(self.label)

    def __len__(self):
        return self.size

    def __repr__(self):
        k = 24
        s = []
        for i in range(min(self.size, k)):
            w = self.reduced_word(i)
            s += ['[ height ' + str(len(w)) + ' : ' + self.printer(w) + ' ]']
        if self.size > k:
            s += ['', '(... and %s more)' % (self.size - k)]
        return '\n' + '\n'.join(s) + '\n'

    def reduced_word(self, n):
        for i in self.strict_descents(n):
            return self.reduced_word(self.operate(n, i)) + (i,)
        return ()

    def __iter__(self):
        return iter(range(self.size))

    def __getitem__(self, n):
        start = n * self.rank * self.stepsize
        return [self._int(self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize]) for i in range(self.rank)]

    @classmethod
    def _is_weak_ascent(cls, frame):
        return all(b == 0xFF for b in frame)

    @classmethod
    def _is_weak_descent(cls, frame):
        return all(b == 0xFF for b in frame[:-1]) and frame[-1] == 0xFE

    @classmethod
    def _int(cls, step):
        return int.from_bytes(step, byteorder='big', signed=False)

    def element(self, n):
        assert 0 <= n < self.size
        return Vector({QPModuleElement(self, n): 1})

    def operate(self, n, *args):
        for i in args:
            assert i in range(self.rank)
            start = n * self.rank * self.stepsize
            step = self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize]
            if self._is_weak_descent(step) or self._is_weak_ascent(step):
                pass
            else:
                n = self._int(step)
        return n

    def _steps(self, n):
        start = n * self.rank * self.stepsize
        for _ in range(self.rank):
            yield self.frame[start:start + self.stepsize]
            start += self.stepsize

    def strict_ascents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_descent(step):
                if self._int(step) > n:
                    yield i

    def strict_descents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_descent(step):
                if self._int(step) < n:
                    yield i

    def weak_ascents(self, n):
        for i, step in enumerate(self._steps(n)):
            if self._is_weak_ascent(step):
                yield i

    def weak_descents(self, n):
        for i, step in enumerate(self._steps(n)):
            if self._is_weak_descent(step):
                yield i

    def __init__(self, family, rank, size, frame=None, printer=None, label=None):
        self.rank = rank
        self.size = size
        self.stepsize = max(1, math.ceil(math.log2((size + 2) / 8.0)))
        self.frame = frame
        self.printer = printer if printer is not None else str
        self.label = (rank, size, family)

    @classmethod
    def create(cls, rank, size, stepsize, minima, operate, verbose=True):
        assert size * rank * stepsize < 0.8e+10

        t0 = time.time()
        if verbose:
            print()
            print('creating bytearray of size %s bytes' % (size * rank * stepsize))
            print()

        frame = bytearray(size * rank * stepsize)

        t1 = time.time()
        if verbose:
            print('* initialized, time elapsed: %s milliseconds' % int(1000 * (t1 - t0)))

        level = {w: [] for w in minima}
        progress, position, start = 0, 0, 0
        while level:
            nextlevel = {}
            for w, origins in level.items():
                assert position < size
                descents = {}
                for (i, n, o) in origins:
                    frame[o:o + stepsize] = position.to_bytes(stepsize, byteorder='big')
                    descents[i] = n
                for i in range(rank):
                    if i in descents:
                        frame[start:start + stepsize] = descents[i].to_bytes(stepsize, byteorder='big')
                    else:
                        y = operate(w, i)
                        if y == (w, True):
                            frame[start:start + stepsize] = bytes(stepsize * [0xFF])
                        elif y == (w, False):
                            frame[start:start + stepsize] = bytes((stepsize - 1) * [0xFF] + [0xFE])
                        else:
                            nextlevel[y] = nextlevel.get(y, []) + [(i, position, start)]
                    start += stepsize
                position += 1
            level = nextlevel

            if verbose:
                print('* level %s done @ position %s of %s' % (progress, position, size))

            progress += 1

        t1 = time.time()
        if verbose:
            print()
            print('finished, time elapsed: %s milliseconds' % int(1000 * (t1 - t0)))
            print()
        return frame

    @classmethod
    def create_hecke_a(cls, n):
        size = math.factorial(n)
        rank = max(0, n - 1)
        module = cls(cls.HECKE_A, rank, size)
        stepsize = module.stepsize

        def multiply(w, i):
            return w * Permutation.s_i(i + 1)

        def printer(word):
            return str(Permutation.from_word([i + 1 for i in word]))

        module.frame = cls.create(rank, size, stepsize, [Permutation()], multiply)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_a(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n

        size = math.factorial(n) // math.factorial(k) // 2**k // math.factorial(n - 2 * k)
        rank = max(0, n - 1)
        module = cls(cls.GELFAND_A(k), rank, size)
        stepsize = module.stepsize

        def conjugate(w, i):
            i += 1
            if w(i) == i and w(i + 1) == i + 1:
                return (w, plus)
            if w(i) == i + 1 and w(i + 1) == i:
                return (w, not plus)
            s = Permutation.s_i(i)
            return s * w * s

        w = Permutation(*[1 + i + (-1)**i for i in range(2 * k)])

        def printer(word):
            v = w
            for i in word:
                s = Permutation.s_i(i + 1)
                v = s * v * s
            return v.cycle_repr()

        module.frame = cls.create(rank, size, stepsize, [w], conjugate)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_bc(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        rank = max(0, n)
        module = cls(cls.GELFAND_BC(k), rank, size)
        stepsize = module.stepsize

        def conjugate(w, i):
            s = SignedPermutation.s_i(i, w.rank)
            if i == 0 and w(1) in [-1, 1]:
                return w * s
            if i > 0 and abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
                return (w, plus)
            if i > 0 and w(i) == i and w(i + 1) == i + 1:
                return (w, not plus)
            if i > 0 and w(i) == -i and w(i + 1) == -i - 1:
                return (w, not plus)
            return s * w * s

        w = SignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [i for i in range(2 * k + 1, n + 1)]
        ))

        def printer(word):
            v = w
            for i in word:
                s = SignedPermutation.s_i(i, w.rank)
                v = s * v * s if not (i == 0 and v(1) == 1) else v * s
            return v.cycle_repr()

        module.frame = cls.create(rank, size, stepsize, [w], conjugate)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_d(cls, n, k, plus=True):
        assert n >= 2
        assert n % 2 != 0
        assert 0 <= 2 * k <= n

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        assert size % 2 == 0
        size = size // 2

        rank = max(0, n)
        module = cls(cls.GELFAND_D(k), rank, size)
        stepsize = module.stepsize

        def conjugate(w, i):
            if i == 0:
                s = SignedPermutation.ds_i(-1, w.rank)
                t = SignedPermutation.ds_i(1, w.rank)
                if abs(w(1)) != 1 and abs(w(2)) != 2:
                    if w * s == s * w:
                        return (w, plus)
                    return s * w * s
                if (w(1) == 1 and w(2) == 2) or (w(1) == -1 and w(2) == -2):
                    return s * w * t
                if (w(1) == 1 and w(2) == -2) or (w(1) == -1 and w(2) == 2):
                    return (w, not plus)
                if (abs(w(1)) == 1 and abs(w(2)) != 2) or (abs(w(1)) != 1 and abs(w(2)) == 2):
                    return (s * w * s).dstar()

            if i > 0 and abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
                return (w, plus)
            if i > 0 and w(i) == i and w(i + 1) == i + 1:
                return (w, not plus)
            if i > 0 and w(i) == -i and w(i + 1) == -i - 1:
                return (w, not plus)

            s = SignedPermutation.ds_i(i, w.rank)
            return s * w * s

        w = SignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [i for i in range(2 * k + 1, n + 1)]
        ))

        def printer(word):
            v = w
            for i in word:
                if i == 0:
                    s = SignedPermutation.ds_i(-1, v.rank)
                    t = SignedPermutation.ds_i(1, v.rank)
                    if abs(v(1)) != 1 and abs(v(2)) != 2:
                        v = s * v * s
                    elif (v(1) == 1 and v(2) == 2) or (v(1) == -1 and v(2) == -2):
                        v = s * v * t
                    else:
                        v = (s * v * s).dstar()
                else:
                    s = SignedPermutation.ds_i(i, v.rank)
                    v = s * v * s
                assert v.is_involution()
            return v.cycle_repr()

        module.frame = cls.create(rank, size, stepsize, [w], conjugate)
        module.printer = printer
        return module
