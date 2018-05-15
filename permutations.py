import itertools


class Permutation:

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield Permutation(args)

    @classmethod
    def involutions(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            w = Permutation(args)
            if all(w(w(i)) == i for i in range(1, n + 1)):
                yield w

    def shift(self, n):
        assert n >= 0
        oneline = [i + 1 for i in range(n)] + [i + n for i in self.oneline]
        return Permutation(oneline)

    def standardize(self, e):
        index = {b: a + 1 for a, b in enumerate(sorted(map(self, e)))}
        oneline = [index[self(i)] for i in sorted(e)]
        return Permutation(oneline)

    def get_descentset_fpf_visible(self):
        return sorted(
            [ij[0] for ij in self.get_fpf_visible_inversions() if ij[1] == ij[0] + 1]
        )

    def get_descentset_visible(self):
        return sorted(
            [ij[0] for ij in self.get_visible_inversions() if ij[1] == ij[0] + 1]
        )

    def psi(self):
        """See transitions paper, Schur P-positivity section."""
        max_vis_inv = self.max_visible_inversion()
        if not max_vis_inv:
            return None
        q, r = max_vis_inv
        a = set([self(q), self(r), q, r])
        t = Permutation.cycle([q, r])
        if len(a) == 2:
            return self * t
        else:
            return t * self * t

    def get_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j):
                    ans.add((i, j))
        return ans

    def get_visible_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j) and i >= self(j):
                    ans.add((i, j))
        return ans

    def get_fpf_visible_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j) and i > self(j):
                    ans.add((i, j))
        return ans

    def max_inversion(self):
        return max(self.get_inversions())

    def max_visible_inversion(self):
        return max(self.get_visible_inversions())

    def support(self):
        return [i + 1 for i in range(len(self.oneline)) if self.oneline[i] != i + 1]

    def is_dominant(self):
        n = len(self.oneline)
        for i in range(1, n - 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    if self(i) < self(k) < self(j):
                        return False
        return True

    def rothe_diagram(self):
        ans = []
        n = len(self.oneline)
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                if self(i) > self(j):
                    ans += [(i, self(j))]
        return sorted(ans)

    def code_helper(self, diag):
        n = len(self.oneline)
        ans = n * [0]
        for i, j in diag:
            ans[i - 1] += 1
        return ans

    def code(self):
        return self.code_helper(self.rothe_diagram())

    def involution_code(self):
        return self.code_helper(self.involution_rothe_diagram())

    def involution_code_fpf(self):
        return self.code_helper(self.involution_rothe_diagram(True))

    def fpf_rothe_diagram(self, fpf=False):
        return self.involution_rothe_diagram(True)

    def involution_rothe_diagram(self, fpf=False):
        diag = self.rothe_diagram()
        ans = []
        for i, j in diag:
            if i > j or (not fpf and i == j):
                ans += [(i, j)]
        return ans

    def __init__(self, oneline=[]):
        while len(oneline) > 0 and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        self.oneline = oneline
        self.cycles = self.get_cycles(oneline)

    def is_left_descent(self, i):
        return self(i) > self(i + 1)

    def get_descent_L(self):
        for i in range(1, len(self.oneline)):
            if(self(i) > self(i + 1)):
                return i
        return 0

    def get_descent_R(self):
        return (self.inverse()).get_descent_L()

    def get_descentset_T(self):
        ans = []
        for i in range(1, len(self.oneline)):
            if(self(i) > self(i + 1) and self(i) <= i and self(i + 1) < i + 1):
                ans.append(i)
        return ans

    def left_descent_set(self):
        return self.get_descentset_L()

    def right_descent_set(self):
        return self.get_descentset_R()

    def get_descentset_L(self):
        ans = []
        for i in range(1, len(self.oneline)):
            if(self(i) > self(i + 1)):
                ans.append(i)
        return ans

    def get_descentset_R(self):
        return (self.inverse()).get_descentset_L()

    def number_two_cycles(self):
        ans = 0
        for c in self.cycles:
            if len(c) == 2:
                ans += 1
        return ans

    def contains_atom(self, w):
        expr = w.reduced_expr()
        return self.expr_to_involution(expr) == self

    def is_atom(self):
        expr = self.reduced_expr()
        return self.expr_to_involution(expr).involution_length() == len(expr)

    def get_min_atom(self):
        cycles = sorted([list(reversed(sorted(c))) for c in self.cycles], key=lambda x: x[-1])
        return Permutation([i for cycle in cycles for i in cycle])**-1

    def is_involution(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 2

    def is_identity(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 1

    # Input is [i_1, i_2, ... , i_k], returns permutation (i_1 i_2 ... i_k)
    @staticmethod
    def cycle(cyc):
        return Permutation.cycles([cyc])

    @staticmethod
    def cycles(cyc, oneline=None):

        if len(cyc) == 0:
            return Permutation()

        if oneline is None:
            n = max(map(max, cyc))
            oneline = list(range(1, 1 + n))
            for c in cyc:
                oneline[c[len(c) - 1] - 1] = c[0]
                for i in range(len(c) - 1):
                    oneline[c[i] - 1] = c[i + 1]
        ans = Permutation()
        ans.oneline = oneline
        ans.cycles = cyc
        return ans

    @classmethod
    def s_i(cls, i):
        return cls.cycle([i, i + 1])

    def inverse(self):
        oneline = list(range(1, 1 + len(self.oneline)))
        for i in range(1, 1 + len(self.oneline)):
            oneline[self(i) - 1] = i
        return Permutation(oneline)

    def __pow__(self, e):
        if e < 0:
            return self.inverse()**(-e)
        m = 1
        for ell in list(map(len, self.cycles)):
                m = m * ell
        e = e % m
        x = Permutation()
        for i in range(e):
            x = x * self
        return x

    # other is a permutation
    def __mul__(self, other):
        newline = []
        n = max(len(self.oneline), len(other.oneline))
        for i in range(1, 1 + n):
            newline.append(other(self(i)))
        while len(newline) > 0 and newline[-1] == len(newline):
            newline = newline[:-1]
        return Permutation(newline)

    def __len__(self):
        if(len(self.cycles) == 0):
            return 0
        return len(self.get_inversions())

    def length(self):
        return len(self)

    def __call__(self, i):
        if i < 1:
            return i
        if i > len(self.oneline):
            return i
        return self.oneline[i - 1]

    def get_cycles(self, oneline):
        cycles = []
        numbers = list(range(1, len(oneline) + 1))
        while len(numbers) > 0:
            cycle = [numbers[0]]
            next = oneline[numbers[0] - 1]
            numbers.remove(numbers[0])

            while next != cycle[0]:
                cycle.append(next)
                numbers.remove(next)
                next = oneline[cycle[len(cycle) - 1] - 1]
            cycles.append(cycle)
        return cycles

    # Lexicographic order on one-line representation
    def __lt__(self, other):
        for i in range(1, 1 + max(max(self.oneline), max(other.oneline))):
            if self(i) >= other(i):
                return False
        return True

    def __eq__(self, other):
        if not isinstance(other, Permutation):
            return False
        return len(self * (other**-1)) == 0

    def __ne__(self, other):
        return not (self == other)

    def strong_bruhat_less_than(self, other):
        if self.length() >= other.length():
            return False
        des = self.get_descentset_L()
        for k in des:
            x = sorted(map(self, range(1, k + 1)))
            y = sorted(map(other, range(1, k + 1)))
            for i in range(k):
                if x[i] > y[i]:
                    return False
        return True

    def __repr__(self):
        # return str(self.oneline)
        return self.cycle_repr()

    def cycle_repr(self):
        if len(self) == 0:
            return '1'
        EXCLUDE_SINGLE_CYCLES = True
        SPACE = ' '
        DELIM = ''  # ','

        s = ''
        for c in self.cycles:
            if not EXCLUDE_SINGLE_CYCLES or len(c) > 1:
                s += '(' + (DELIM + SPACE).join([str(x) for x in c]) + ')'
        return s

    def __hash__(self):
        return hash(str(self))

    def involution_length(self):
        return (self.length() + len(list(filter(lambda i: len(i) > 1, self.cycles)))) // 2
