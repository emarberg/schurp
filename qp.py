from permutations import Permutation
from signed import SignedPermutation
from even import EvenSignedPermutation
import polynomials
from vectors import Vector
from heapq import heappush, heappop
import json
import math
import time


class QPWGraph:

    def __init__(self, qpmodule, nbytes=8, verbose=True):
        self.qpmodule = qpmodule
        self.nbytes = 8
        self.frame = None

        t0 = time.time()
        if verbose:
            a = self.qpmodule.get_filename().split('/')[-1]
            b = len(self.qpmodule)
            print()
            print('QPWGraph for %s (%s elements)' % (a, b))
            print()
            print('* computing heights', end='') # noqa

        self.heights = []
        for n in self.qpmodule:
            nht = self.qpmodule.height(n)
            if nht >= len(self.heights):
                self.heights.append(1)
            else:
                self.heights[nht] += 1

        t1 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t1 - t0)))
            print('* computing size', end='') # noqa

        self.size = 0
        for i in range(len(self.heights)):
            for j in range(i + 1, len(self.heights)):
                self.size += self._space(i, j) * self.heights[i] * self.heights[j]
        assert self.size < 0.8e+10

        t2 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t2 - t1)))
            print('* computing suboffsets', end='') # noqa

        self.hbytes = max(1, math.ceil(math.log2(max(self.heights) / 8.0)))
        self.suboffsets = bytearray(self.hbytes * len(self.qpmodule))

        offset = 0
        start = 0
        for n in self.qpmodule:
            if offset > 0 and self.qpmodule.height(n) > self.qpmodule.height(n - 1):
                offset = 0
            v = self.heights[self.qpmodule.height(n)] - 1 - offset
            self.suboffsets[start:start + self.hbytes] = v.to_bytes(self.hbytes, byteorder='big')
            start += self.hbytes
            offset += 1

        t3 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t3 - t2)))
            print('* computing height offsets', end='') # noqa

        self.offsets = {}
        total_offsets = {}
        for j in range(len(self.heights)):
            self.offsets[j] = (j + 1) * [0]
            for i in range(j - 2, -1, -1):
                self.offsets[j][i] = self.offsets[j][i + 1] + self.heights[i + 1] * self._space(i + 1, j)
            total_offsets[j] = self.offsets[j][0] + self.heights[0] * self._space(0, j)

        t4 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t4 - t3)))
            print('* computing addresses', end='') # noqa

        self.abytes = max(1, math.ceil(math.log2((self.size) / 8.0)))
        self.addresses = bytearray(self.abytes * len(self.qpmodule))

        a = 0
        start = 0
        for j in self.qpmodule:
            self.addresses[start:start + self.abytes] = a.to_bytes(self.abytes, byteorder='big')
            start += self.abytes
            a += total_offsets[self.qpmodule.height(j)]

        t5 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t5 - t4)))
            print('* finished')
            print()

    def height(self, n):
        return self.qpmodule.height(n)

    def _space(self, i, j):
        return (1 + (j - i - 1) // 2) * self.nbytes

    @classmethod
    def _int(cls, step, signed=False):
        return int.from_bytes(step, byteorder='big', signed=signed)

    def __repr__(self):
        a = self.qpmodule.get_filename().split('/')[-1]
        b = len(self.qpmodule)
        s = ['QPWGraph for %s (%s elements)' % (a, b)]
        s += ['*        nbytes = %s' % str(self.nbytes)]
        s += ['*        abytes = %s' % str(self.abytes)]
        s += ['*  address size = %s' % str(self.abytes * b)]
        s += ['*          size = %s' % str(self.size)]
        return '\n'.join(s)

    def address_cbasis(self, i, j):
        hi = self.height(i)
        hj = self.height(j)

        assert hi < hj
        space = self._space(hi, hj)

        start_i = i * self.hbytes
        start_j = j * self.abytes

        addr = self._int(self.addresses[start_j:start_j + self.abytes]) + \
            self.offsets[hj][hi] + \
            self._int(self.suboffsets[start_i:start_i + self.hbytes]) * space
        return addr, addr + space, space

    def get_cbasis_leading(self, i, j):
        if self.height(i) >= self.height(j):
            return 0
        _, stop, _ = self.address_cbasis(i, j)
        return self._int(self.frame[stop - self.nbytes:stop], signed=True)

    def get_cbasis(self, i, j, return_bytes=True):
        hi = self.height(i)
        hj = self.height(j)
        space = self._space(hi, hj)

        if i == j:
            return (1).to_bytes(self.nbytes, byteorder='big', signed=True) if return_bytes else 1
        elif hi >= hj:
            return (0).to_bytes(self.nbytes, byteorder='big', signed=True) if return_bytes else 0
        start, stop, _ = self.address_cbasis(i, j)
        return self.frame[start:stop] if return_bytes else self._int(self.frame[start:stop], signed=True)

    def get_cbasis_polynomial(self, i, j):
        if i == j:
            return polynomials.q(0)
        c = self.get_cbasis(i, j)
        ans = 0 * polynomials.q(0)
        start = 0
        for e in range(1 + (j - i - 1) // 2):
            ans += polynomials.q(e) * self._int(c[start:start + self.nbytes], signed=True)
            start += self.nbytes
        return ans

    def set_cbasis(self, i, j, v, set_bytes=True):
        start, stop, size = self.address_cbasis(i, j)
        if set_bytes:
            self.frame[start:stop] = v
        else:
            self.frame[start:stop] = v.to_bytes(size, byteorder='big', signed=True)

    def _safe_add(self, start, f, shift=0, mu=1):
        if mu == 0:
            return
        astart = start + shift * self.nbytes
        fstart = 0
        while fstart < len(f):
            a = self._int(self.frame[astart:astart + self.nbytes], signed=True)
            b = self._int(f[fstart:fstart + self.nbytes], signed=True)
            v = (a + mu * b).to_bytes(self.nbytes, byteorder='big', signed=True)
            self.frame[astart:astart + self.nbytes] = v
            astart += self.nbytes
            fstart += self.nbytes

    def compute(self, verbose=True):
        t0 = time.time()
        self.frame = bytearray(self.size)

        if verbose:
            print('Compututing canonical basis:')
            progress = 0

        for j in self.qpmodule:
            hj = self.height(j)
            wdes = set(self.qpmodule.weak_descents(j))
            sdes = set(self.qpmodule.strict_descents(j))
            des = wdes | sdes

            for i in range(j - 1, -1, -1):
                hi = self.height(i)
                if hi == hj:
                    continue

                if des & set(self.qpmodule.weak_ascents(i)):
                    continue

                start, _, _ = self.address_cbasis(i, j)

                if verbose:
                    newprogress = int(100 * start / self.size)
                    if newprogress > progress:
                        print('*', newprogress, 'percent done (%s milliseconds elapsed)' % int(1000 * (time.time() - t0)))
                        progress = newprogress

                t = set(self.qpmodule.strict_ascents(i)) & des
                if t:
                    x = self.qpmodule.operate(i, next(iter(t)))
                    self._safe_add(start, self.get_cbasis(x, j))
                    continue

                s = next(iter(sdes))
                si = self.qpmodule.operate(i, s)
                sj = self.qpmodule.operate(j, s)

                self._safe_add(start, self.get_cbasis(si, sj))
                self._safe_add(start, self.get_cbasis(i, sj), shift=1)
                for x in range(i, sj):
                    if s in self.qpmodule.weak_ascents(x) or s in self.qpmodule.strict_ascents(x):
                        continue
                    hx = self.height(x)
                    if (hj - hx) % 2 != 0:
                        continue
                    mu = self.get_cbasis_leading(x, sj)
                    self._safe_add(start, self.get_cbasis(i, x), (hj - hx) // 2, -mu)

        if verbose:
            print('Done computing (%s milliseconds elapsed)' % int(1000 * (time.time() - t0)))


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
        return str(self.qpmodule.permutation(self.n))

    def __mul__(self, i):
        result = self.qpmodule.operate(self.n, i)
        if result == self.n:
            if i in self.qpmodule.weak_ascents(self.n):
                return Vector({self: polynomials.q(1)})
            else:
                return Vector({self: -polynomials.q(-1)})
        new = QPModuleElement(self.qpmodule, result)
        if result < self.n:
            return Vector({self: polynomials.q(1) - polynomials.q(-1), new: 1})
        else:
            return Vector({new: 1})


class QPModule:

    DIRECTORY = '/Users/emarberg/Desktop/examples/models/'

    HECKE_A = 'HECKE_A'
    GELFAND_A = 'GELFAND_A'
    GELFAND_BC = 'GELFAND_BC'
    GELFAND_D = 'GELFAND_D'

    def expand(self, n):
        if self.printer:
            return self.printer(self.reduced_word(n))
        else:
            word = self.reduced_word(n)
            return len(word), word

    def permutation(self, n):
        return self.expand(n)[1]

    def length(self, n):
        return self.height(n)

    def height(self, n):
        return self[n][0]

    def __eq__(self, other):
        if self.layer != other.layer:
            return False
        if self.family != other.family or self.rank != other.rank or self.size != other.size:
            return False
        if self.stepsize != other.stepsize or self.height_bytes != other.height_bytes:
            return False
        return self.frame == other.frame

    def __hash__(self):
        return hash((self.family, self.rank, self.layer))

    def __len__(self):
        return self.size

    def __repr__(self):
        k = 24
        s = [self.family + str(self.rank) + ' : LAYER ' + str(self.layer)]
        s += ['']
        if self.frame:
            for i in range(min(self.size, k)):
                ht, w = self.expand(i)
                des = str(list(self.weak_descents(i)))
                asc = str(list(self.weak_ascents(i)))
                s += ['[ height ' + str(ht) + ' : ' + str(w) + ' : weak descents ' + des + ' : weak ascents ' + asc + ' ]']
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
        start = n * (self.height_bytes + self.rank * self.stepsize)
        ans = [self._int(self.frame[start:start + self.height_bytes])]
        start += self.height_bytes
        return ans + [
            self._int(self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize])
            for i in range(self.rank)
        ]

    @classmethod
    def _is_weak_ascent(cls, frame):
        return all(b == 0xFF for b in frame)

    @classmethod
    def _is_weak_descent(cls, frame):
        return all(b == 0xFF for b in frame[:-1]) and frame[-1] == 0xFE

    @classmethod
    def _int(cls, step, signed=False):
        return int.from_bytes(step, byteorder='big', signed=signed)

    def element(self, n):
        assert 0 <= n < self.size
        return Vector({QPModuleElement(self, n): 1})

    def operate(self, n, *args):
        for i in args:
            assert i in range(self.rank)
            start = n * (self.height_bytes + self.rank * self.stepsize) + self.height_bytes
            step = self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize]
            if self._is_weak_descent(step) or self._is_weak_ascent(step):
                pass
            else:
                n = self._int(step)
        return n

    def _steps(self, n):
        start = n * (self.height_bytes + self.rank * self.stepsize) + self.height_bytes
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

    def __init__(self, family, rank, layer, size, height_bytes, frame=None, printer=None):
        self.family = family
        self.rank = rank
        self.layer = layer
        self.size = size
        self.stepsize = max(1, math.ceil(math.log2((size + 2) / 8.0)))
        self.height_bytes = height_bytes
        self.frame = frame
        self.printer = printer

    def get_filename(self):
        return self.filename(self.family, self.rank, self.layer)

    @classmethod
    def filename(cls, family, rank, layer):
        file = cls.DIRECTORY + family + str(rank)
        if layer is not None:
            file += '_' + str(layer)
        return file

    def metadata(self):
        return {
            'family': self.family,
            'rank': self.rank,
            'layer': self.layer,
            'size': self.size,
            'stepsize': self.stepsize,
            'height_bytes': self.height_bytes,
        }

    def write(self):
        filename = self.get_filename()
        bytefile = filename + '.b'
        with open(bytefile, 'wb') as file:
            file.write(self.frame)
            file.close()
        metafile = filename + '.meta'
        with open(metafile, 'w') as file:
            file.write(json.dumps(self.metadata()))
            file.close()

    @classmethod
    def read_gelfand_a(cls, n, k):
        filename = cls.filename(cls.GELFAND_A, n, k)
        return cls.read(filename)

    @classmethod
    def read_gelfand_bc(cls, n, k):
        filename = cls.filename(cls.GELFAND_BC, n, k)
        return cls.read(filename)

    @classmethod
    def read_gelfand_d(cls, n, k):
        filename = cls.filename(cls.GELFAND_D, n, k)
        return cls.read(filename)

    @classmethod
    def read(cls, filename):
        metafile = filename + '.meta'
        with open(metafile, 'r') as file:
            dictionary = json.loads(file.read())
            file.close()

        family = dictionary['family']
        rank = int(dictionary['rank'])
        layer = dictionary['layer']
        layer = None if layer == 'None' else int(layer)
        size = int(dictionary['size'])
        stepsize = int(dictionary['stepsize'])
        height_bytes = int(dictionary['height_bytes'])
        module = cls(family, rank, layer, size, height_bytes)
        assert module.stepsize == stepsize

        bytefile = filename + '.b'
        with open(bytefile, 'rb') as file:
            module.frame = bytearray(file.read())
            file.close()
        return module

    @classmethod
    def create(cls, rank, size, stepsize, height_bytes, minima, operate, verbose=True):
        assert size * (rank * stepsize + height_bytes) < 0.8e+10

        t0 = time.time()
        if verbose:
            print()
            print('RANK', rank, 'SIZE', size)
            print()
            print('* creating bytearray of size %s bytes' % (size * rank * stepsize))

        frame = bytearray(size * (rank * stepsize + height_bytes))
        level = [(wht, time.time(), w) for (wht, w) in minima]
        origins = {w: [] for (_, w) in minima}
        position, start = 0, 0

        while level:
            wht, _, w = heappop(level)
            assert position < size
            frame[start:start + height_bytes] = wht.to_bytes(height_bytes, byteorder='big')
            start += height_bytes
            descents = {}
            for (i, n, o) in origins[w]:
                frame[o:o + stepsize] = position.to_bytes(stepsize, byteorder='big')
                descents[i] = n
            for i in range(rank):
                if i in descents:
                    frame[start:start + stepsize] = descents[i].to_bytes(stepsize, byteorder='big')
                else:
                    ht, y = operate(wht, w, i)
                    if y == (w, True):
                        frame[start:start + stepsize] = bytes(stepsize * [0xFF])
                    elif y == (w, False):
                        frame[start:start + stepsize] = bytes((stepsize - 1) * [0xFF] + [0xFE])
                    else:
                        if y not in origins:
                            heappush(level, (ht, time.time(), y))
                        origins[y] = origins.get(y, []) + [(i, position, start)]
                start += stepsize
            position += 1
            del origins[w]

        t1 = time.time()
        if verbose:
            print('* finished, time elapsed: %s milliseconds' % int(1000 * (t1 - t0)))
            print()
        return frame

    @classmethod
    def create_hecke_a(cls, n):
        assert n >= 0
        size = math.factorial(n + 1)
        height_bytes = 1
        module = cls(cls.HECKE_A, n, None, size, height_bytes)
        stepsize = module.stepsize

        def multiply(ht, w, i):
            ht += (1 if w(i + 1) < w(i + 2) else -1)
            return ht, w * Permutation.s_i(i + 1)

        def printer(word):
            ht, v = 0, Permutation()
            for i in word:
                ht, v = multiply(ht, v, i)
            return ht, v

        minima = [(0, Permutation())]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, multiply)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_a(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n + 1

        size = math.factorial(n + 1) // math.factorial(k) // 2**k // math.factorial(n + 1 - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_A, n, k, size, height_bytes)
        stepsize = module.stepsize

        def conjugate(ht, w, i):
            i += 1
            if w(i) == i and w(i + 1) == i + 1:
                return ht, (w, plus)
            if w(i) == i + 1 and w(i + 1) == i:
                return ht, (w, not plus)
            s = Permutation.s_i(i)
            return ht + 1, s * w * s

        w = Permutation(*[1 + i + (-1)**i for i in range(2 * k)])

        def printer(word):
            ht, v = 0, w
            for i in word:
                ht, v = conjugate(ht, v, i)
            return ht, v

        minima = [(0, w)]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, conjugate)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_bc(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_BC, n, k, size, height_bytes)
        stepsize = module.stepsize

        def conjugate(ht, w, i):
            s = SignedPermutation.s_i(i, w.rank)
            if i == 0 and w(1) in [-1, 1]:
                return ht + 1, w * s
            if i > 0 and abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
                return ht, (w, not plus)
            if i > 0 and w(i) == i and w(i + 1) == i + 1:
                return ht, (w, plus)
            if i > 0 and w(i) == -i and w(i + 1) == -i - 1:
                return ht, (w, plus)
            return ht + 1, s * w * s

        w = SignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [i for i in range(2 * k + 1, n + 1)]
        ))

        def printer(word):
            ht, v = 0, w
            for i in word:
                ht, v = conjugate(ht, v, i)
            return ht, v

        minima = [(0, w)]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, conjugate)
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

        height_bytes = 1
        module = cls(cls.GELFAND_D, n, k, size, height_bytes)
        stepsize = module.stepsize

        def conjugate(ht, w, i):
            if i == 0:
                s = SignedPermutation.ds_i(-1, w.rank)
                t = SignedPermutation.ds_i(1, w.rank)
                if abs(w(1)) != 1 and abs(w(2)) != 2:
                    if w * s == s * w:
                        return ht, (w, not plus)
                    return ht + 1, s * w * s
                if (w(1) == 1 and w(2) == 2) or (w(1) == -1 and w(2) == -2):
                    return ht + 1, s * w * t
                if (w(1) == 1 and w(2) == -2) or (w(1) == -1 and w(2) == 2):
                    return ht, (w, plus)
                if (abs(w(1)) == 1 and abs(w(2)) != 2) or (abs(w(1)) != 1 and abs(w(2)) == 2):
                    return ht + 1, (s * w * s).dstar()
                raise Exception

            if abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
                return ht, (w, not plus)
            if w(i) == i and w(i + 1) == i + 1:
                return ht, (w, plus)
            if w(i) == -i and w(i + 1) == -i - 1:
                return ht, (w, plus)

            s = SignedPermutation.ds_i(i, w.rank)
            return ht + 1, s * w * s

        w = SignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [i for i in range(2 * k + 1, n + 1)]
        ))

        def printer(word):
            ht, v = 0, w
            for i in word:
                ht, v = conjugate(ht, v, i)
            return ht, v

        minima = [(0, w)]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, conjugate)
        module.printer = printer
        return module

    @classmethod
    def slow_create_gelfand_a(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n + 1

        m = 2 * n + 2 - 2 * k
        a = [Permutation.s_i(i) for i in range(n + 2, m)]
        s = {i: Permutation.s_i(i + 1) for i in range(n + 1)}
        w = Permutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [n + 2 + i for i in range(n + 1 - 2 * k)] +
            [2 * k + 1 + i for i in range(n + 1 - 2 * k)]
        ))

        size = math.factorial(n + 1) // math.factorial(k) // 2**k // math.factorial(n + 1 - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_A, n, k, size, height_bytes)
        return cls.create_gelfand_classical(plus, module, a, s, w)

    @classmethod
    def slow_create_gelfand_bc(cls, n, k, plus=True):
        assert 0 <= 2 * k <= n

        m = 2 * n - 2 * k
        a = [SignedPermutation.s_i(i, m) for i in range(n + 1, m)]
        s = {i: SignedPermutation.s_i(i, m) for i in range(n)}
        w = SignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [n + 1 + i for i in range(n - 2 * k)] +
            [2 * k + 1 + i for i in range(n - 2 * k)]
        ))

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_BC, n, k, size, height_bytes)
        return cls.create_gelfand_classical(plus, module, a, s, w)

    @classmethod
    def slow_create_gelfand_d(cls, n, k, plus=True):
        assert n >= 2
        assert n % 2 != 0
        assert 0 <= 2 * k <= n

        m = 2 * n - 2 * k
        a = [EvenSignedPermutation.s_i(i, m) for i in range(n + 1, m)]
        s = {i: EvenSignedPermutation.s_i(i, m) for i in range(n)}
        w = EvenSignedPermutation(*(
            [1 + i + (-1)**i for i in range(2 * k)] +
            [n + 1 + i for i in range(n - 2 * k)] +
            [2 * k + 1 + i for i in range(n - 2 * k)]
        ))

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        assert size % 2 == 0
        size = size // 2

        height_bytes = 1
        module = cls(cls.GELFAND_D, n, k, size, height_bytes)
        return cls.create_gelfand_classical(plus, module, a, s, w)

    @classmethod
    def create_gelfand_classical(cls, plus, module, a, s, w):

        def conjugate(ht, u, i):
            if u * s[i] * u in a:
                return ht, (u, plus)
            v = s[i] * u * s[i]
            if v == u:
                return ht, (u, not plus)
            return ht + 1, v

        def printer(word):
            ht, v = 0, w
            for i in word:
                ht, v = conjugate(ht, v, i)
            return ht, v

        minima = [(0, w)]
        module.frame = cls.create(
            module.rank,
            module.size,
            module.stepsize,
            module.height_bytes,
            minima,
            conjugate
        )
        module.printer = printer
        return module
