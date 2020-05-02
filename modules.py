from permutations import Permutation
import math
import time


class QPModule:

    def __repr__(self):
        k = 10
        if self.size > k:
            s = '\n'.join([self.printer(self.reduced_word(i)) for i in range(k)]) + '\n\n(... and %s more)' % (self.size - k)
        else:
            s = '\n'.join([self.printer(self.reduced_word(i)) for i in self])
        return '\n' + s + '\n'

    def reduced_word(self, n):
        for i in self.strict_descents(n):
            return (i,) + self.reduced_word(self.act(i, n))
        return ()

    def __iter__(self):
        return iter(range(self.size))

    def __getitem__(self, n):
        start = n * self.rank * self.stepsize
        return [self._int(self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize]) for i in range(self.rank)]

    @classmethod
    def _is_weak_ascent(cls, frame):
        return all(b == 0b1111111 for b in frame)

    @classmethod
    def _is_weak_decent(cls, frame):
        return all(b == 0b0000000 for b in frame[:-1]) and frame[-1] == 0b11111110

    @classmethod
    def _int(cls, step):
        return int.from_bytes(step, byteorder='big', signed=False)

    def act(self, i, n):
        start = n * self.rank * self.stepsize
        step = self.frame[start  + i * self.stepsize:start + (i + 1) * self.stepsize]
        if self._is_weak_decent(step) or self._is_weak_ascent(step):
            return n
        return self._int(step)

    def _steps(self, n):
        start = n * self.rank * self.stepsize
        for _ in range(self.rank):
            yield self.frame[start:start + self.stepsize]
            start += self.stepsize

    def strict_ascents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_decent(step):
                if self._int(step) > n:
                    yield i

    def strict_descents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_decent(step):
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

    def __init__(self, rank, size, frame=None, printer=None):
        self.rank = rank
        self.size = size
        self.stepsize = max(1, math.ceil(math.log2((size + 2) / 8.0)))
        self.frame = frame
        self.printer = printer if printer is not None else str

    @classmethod
    def create(cls, rank, size, stepsize, minima, act, verbose=True):
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
            print()

        level = {w: [] for w in minima}
        progress, position, start = 0, 0, 0
        while level:
            nextlevel = {}
            for w, origins in level.items():
                descents = {}
                for (i, n, o) in origins:
                    frame[o:o + stepsize] = position.to_bytes(stepsize, byteorder='big')
                    descents[i] = n
                for i in range(rank):
                    if i in descents:
                        frame[start:start + stepsize] = descents[i].to_bytes(stepsize, byteorder='big')
                    else:
                        y = act(i, w)
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
                print('* level %s done @ position %s' % (progress, position))

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
        module = cls(rank, size)
        stepsize = module.stepsize

        def multiply(i, w):
            return Permutation.s_i(i + 1) * w

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
        module = cls(rank, size)
        stepsize = module.stepsize
        w = Permutation(*[1 + i + (-1)**i for i in range(2 * k)])

        def conjugate(i, w):
            i += 1
            if w(i) == i and w(i + 1) == i + 1:
                return (w, plus)
            if w(i) == i + 1 and w(i + 1) == i:
                return (w, not plus)
            s = Permutation.s_i(i)
            return s * w * s

        def printer(word):
            v = w
            for i in reversed(word):
                s = Permutation.s_i(i + 1)
                v = s * v * s
            return v.cycle_repr()

        module.frame = cls.create(rank, size, stepsize, [w], conjugate)
        module.printer = printer
        return module

    @classmethod
    def create_gelfand_bc(cls, n):
        pass

    @classmethod
    def create_gelfand_d(cls, n):
        pass
