from .cached import cached_value
from collections import defaultdict
import itertools

PARTITIONS = {}
RIM_CACHE = {}

DECOMPOSE_RIM_CACHE = {}
DECOMPOSE_VSTRIP_CACHE = {}


class Partition:

    FRENCH = True

    @classmethod
    def stabilizer_order(cls, mu):
        mu = cls.sort(mu, trim=True)
        m = 1
        end = 0
        for i in range(1, len(mu)):
            if mu[i] == mu[end]:
                m *= (i + 1 - end)
            else:
                end = i
        return m

    @classmethod
    def combine(cls, mu, nu):
        n = max(len(mu), len(nu))
        ans = [cls.get(mu, i) + cls.get(nu, i) for i in range(1, n + 1)]
        return Partition.trim(ans)

    @classmethod
    def add(cls, mu, row):
        if row <= len(mu):
            mu = list(mu)
            mu[row - 1] += 1
            return tuple(mu)
        elif row == len(mu) + 1:
            mu = tuple(mu) + (1,)
            return mu
        raise Exception

    @classmethod
    def merge_shifted(cls, mu, shape):
        if type(mu) in [list, tuple]:
            mu = cls.shifted_shape(mu)
        return cls.from_shape(mu | set(shape))

    @classmethod
    def from_shape(cls, sh):
        ans = []
        for (i, j) in sh:
            while i > len(ans):
                ans += [0]
            ans[i - 1] += 1
        return tuple(ans)

    @classmethod
    def shape(cls, mu):
        ans = set()
        for i, a in enumerate(mu):
            for j in range(a):
                ans.add((i + 1, j + 1))
        return ans

    @classmethod
    def shifted_shape(cls, mu, nu=None):
        if nu is None:
            ans = set()
            for i, a in enumerate(mu):
                for j in range(a):
                    ans.add((i + 1, i + j + 1))
            return ans
        return cls.shifted_shape(mu) - cls.shifted_shape(nu)

    @classmethod
    def remove_inner_corners(cls, mu):
        rows = []
        for i in range(1, len(mu) + 1):
            try:
                cls.remove_inner_corner(mu, i)
                rows += [i]
            except:
                continue
        for k in range(len(rows) + 1):
            for subset in itertools.combinations(rows, k):
                ans = mu
                for row in subset:
                    ans = cls.remove_inner_corner(ans, row)
                yield ans

    @classmethod
    def remove_inner_corner(cls, mu, row):
        if row == 1 and len(mu) > 1 and mu[1] < mu[0]:
            return (mu[0] - 1,) + mu[1:]
        elif 1 < row < len(mu) and mu[row] < mu[row - 1]:
            ans = list(mu)
            ans[row - 1] -= 1
            return tuple(ans)
        elif row == len(mu) and mu[-1] > 1:
            return mu[:-1] + (mu[-1] - 1,)
        elif row == len(mu) and mu[-1] == 1:
            return mu[:-1]
        raise Exception

    @classmethod
    def find_shifted_outer_corner(cls, mu, diagonal):
        shape = cls.shifted_shape(mu)
        i, j = 1, 1 + diagonal
        while (i, j) in shape:
            i, j = i + 1, j + 1
        if (i == 1 or (i - 1, j) in shape) and (i == j or (i, j - 1) in shape):
            return i
        return None

    @classmethod
    def find_shifted_inner_corner(cls, mu, diagonal):
        shape = cls.shifted_shape(mu)
        i, j = 0, diagonal
        while (i + 1, j + 1) in shape:
            i, j = i + 1, j + 1
        if (i, j) not in shape:
            return None
        if (i + 1, j) in shape or (i, j + 1) in shape:
            return None
        return i

    @classmethod
    def remove_shifted_inner_corners(cls, mu):
        rows = []
        for i in range(1, len(mu) + 1):
            if cls.remove_shifted_inner_corner(mu, i) is not None:
                rows += [i]
        for k in range(len(rows) + 1):
            for subset in itertools.combinations(rows, k):
                ans = mu
                for row in subset:
                    ans = cls.remove_shifted_inner_corner(ans, row)
                yield ans

    @classmethod
    def remove_shifted_inner_corner(cls, mu, row):
        if row == 1 and len(mu) > 1 and mu[1] + 1 < mu[0]:
            return (mu[0] - 1,) + mu[1:]
        elif 1 < row < len(mu) and mu[row] + 1 < mu[row - 1]:
            ans = list(mu)
            ans[row - 1] -= 1
            return tuple(ans)
        elif row == len(mu) and mu[-1] > 1:
            return mu[:-1] + (mu[-1] - 1,)
        elif row == len(mu) and mu[-1] == 1:
            return mu[:-1]

    @classmethod
    def complement(cls, n, mu):
        assert all(mu[i] <= n - i for i in range(len(mu)))
        dictionary = defaultdict(int)
        for i in range(n):
            start = i + (mu[i] if i < len(mu) else 0) + 1
            for j in range(start, n + 1):
                dictionary[j] += 1
        return Partition.sort(dictionary.values(), trim=True)

    @classmethod
    def print_diagram(cls, mu, nu=None):
        print()
        print(cls.printable(mu, nu))

    @classmethod
    def print_shifted(cls, mu, nu=None):
        print()
        print(cls.printable(mu, nu, shifted=True))

    @classmethod
    def printable(cls, nu, mu=None, shifted=False):
        if nu is None:
            return
        s = []
        for i, a in enumerate(nu):
            x = 0 if mu is None else cls.get(mu, i + 1)
            b = [(i * '  ' if shifted else '') + x * '. ' + (a - x) * '* ']
            s = (b + s) if cls.FRENCH else (s + b)
        if s:
            m = max(map(len, s))
            s = [line + (m - len(line)) * ' ' for line in s]
        return '\n'.join(s)

    @classmethod
    def printables(cls, *args):
        return cls._printables(list(args))

    @classmethod
    def shifted_printables(cls, *args):
        return cls._printables(list(args), shifted=True)

    @classmethod
    def _printables(cls, partitions, shifted=False):
        gap = '   ->    '
        lines = []
        for mu in partitions:
            diagram = cls.printable(mu, shifted).split('\n')
            if lines:
                m = max(len(lines), len(diagram))
                if len(lines) < m:
                    b = (m - len(lines)) * [len(lines[0]) * ' ']
                    lines = (b + lines) if cls.FRENCH else (lines + b)
                if len(diagram) < m:
                    b = (m - len(diagram)) * [(len(diagram[0]) if diagram else 0) * ' ']
                    diagram = (b + diagram) if cls.FRENCH else (diagram + b)
                lines = [lines[i] + gap + diagram[i] for i in range(m)]
            else:
                lines = diagram
        return '\n'.join(lines)

    @classmethod
    def trim(cls, mu):
        while mu and mu[-1] == 0:
            mu = mu[:-1]
        return tuple(mu)

    @classmethod
    def sort(cls, mu, trim=False):
        ans = tuple(reversed(sorted(mu)))
        if trim:
            ans = cls.trim(ans)
        return ans

    @classmethod
    def get(cls, mu, i):
        return mu[i - 1] if 0 <= i - 1 < len(mu) else 0

    @classmethod
    def isin(cls, mu, i, j, shifted=False):
        return 1 <= j - (i - 1 if shifted else 0) <= cls.get(mu, i)

    @classmethod
    def isin_shifted(cls, mu, i, j):
        return cls.isin(mu, i, j, True)

    @classmethod
    def is_partition(cls, mu):
        return all(mu[i - 1] >= mu[i] for i in range(1, len(mu))) and (mu == () or mu[-1] >= 0)

    @classmethod
    def is_strict(cls, mu):
        return cls.is_strict_partition(mu)

    @classmethod
    def is_strict_partition(cls, mu):
        return all(mu[i - 1] > mu[i] for i in range(1, len(mu)))

    @classmethod
    def transpose(cls, mu):
        if mu:
            return tuple(len([i for i in range(len(mu)) if mu[i] > j]) for j in range(mu[0]))
        else:
            return mu

    @classmethod
    def from_skew(cls, positions, shifted=False):
        nu, mu = [], []
        i = 1
        while any(a >= i for (a, b) in positions):
            row = [(a, b) for (a, b) in positions if a == i]
            if row:
                mu.append(min([b for (a, b) in positions if a == i]) - 1 - (i - 1 if shifted else 0))
                nu.append(max([b for (a, b) in positions if a == i]) - (i - 1 if shifted else 0))
            else:
                mu.append(-1)
                nu.append(-1)
            i += 1
        nu.append(0)
        mu.append(0)
        for i in range(len(nu) - 2, -1, -1):
            if nu[i] == -1:
                nu[i] = nu[i + 1] + (1 if shifted else 0)
                mu[i] = nu[i]
        nu, mu = cls.trim(nu), cls.trim(mu)
        while len(nu) == len(mu):
            nu, mu = cls.trim([a - 1 for a in nu]), cls.trim([a - 1 for a in mu])
        return nu, mu

    @classmethod
    def successors(cls, mu, strict=False):
        mu = Partition.trim(mu)
        for i in range(len(mu)):
            if i == 0 or mu[i - 1] + (0 if strict else 1) > mu[i] + 1:
                for nu in cls.successors(mu[i + 1:], strict):
                    yield Partition.trim(mu[:i] + (mu[i] + 1,) + nu)
        yield mu

    @classmethod
    def skew(cls, mu, nu, shifted=False):
        """Assumes nu is contained in mu."""
        ans = set()
        for i, part in enumerate(mu):
            subpart = nu[i] if i < len(nu) else 0
            for j in range(subpart, part):
                ans.add((i + 1, j + 1 + (i if shifted else 0)))
        return ans

    @classmethod
    def skew_key(cls, mu, nu, shifted=False):
        skew = cls.skew(mu, nu, shifted)
        while skew:
            if all(i > 1 and j > 1 for (i, j) in skew):
                skew = {(i - 1, j - 1) for (i, j) in skew}
                continue

            if all(i > 1 for (i, j) in skew) and not (shifted and any(i == j for (i, j) in skew)):
                skew = {(i - 1, j) for (i, j) in skew}
                continue

            if all(j > 1 for (i, j) in skew) and not (shifted and any(i == j - 1 for (i, j) in skew)):
                skew = {(i, j - 1) for (i, j) in skew}
                continue
            break
        return tuple(sorted(skew))

    @classmethod
    def contains(cls, bigger, smaller):
        """Returns true if mu subseteq nu as partitions."""
        if len(smaller) > len(bigger):
            return False
        return all(0 <= smaller[i] <= bigger[i] for i in range(len(smaller)))

    @classmethod
    def compositions(cls, n):
        if n == 0:
            return [()]
        ans = []
        for i in range(1, n + 1):
            for e in cls.compositions(n - i):
                ans.append((i,) + e)
        return ans

    @classmethod
    def peak_compositions(cls, n):
        return [alpha for alpha in cls.compositions(n) if all(a >= 2 for a in alpha[:-1])]

    @classmethod
    def all(cls, n, max_part=None, max_row=None, strict=False, even_parts=False):
        for i in range(n + 1):
            for mu in cls.generate(i, max_part=max_part, max_row=max_row, strict=strict, even_parts=even_parts):
                yield mu

    @classmethod
    def generate(cls, n, max_part=None, max_row=None, strict=False, even_parts=False):
        if n == 0:
            yield ()
        else:
            if (n, max_part, max_row, strict, even_parts) not in PARTITIONS:
                ans = []
                max_part = n if (max_part is None or max_part > n) else max_part
                max_row = n if (max_row is None or max_row > n) else max_row
                if max_row > 0:
                    for i in range(1, max_part + 1):
                        if even_parts and i % 2 != 0:
                            continue
                        for mu in cls.generate(n - i, i, max_row - 1, strict=strict, even_parts=even_parts):
                            nu = (i,) + mu
                            if not strict or Partition.is_strict_partition(nu):
                                ans.append(nu)
                                yield nu
                PARTITIONS[(n, max_part, strict, even_parts)] = ans
            else:
                for mu in PARTITIONS[(n, max_part, strict, even_parts)]:
                    yield mu

    @classmethod
    def subpartitions(cls, mu, nu=(), strict=False):
        lowerbound = Partition.trim(nu)

        def _subpartitions(mu, strict, i):
            if mu:
                for nu in _subpartitions(mu[1:], strict, i + 1):
                    lb = (nu[0] + (1 if strict else 0)) if (nu and nu[0] > 0) else 0
                    lb = max(lb, Partition.get(lowerbound, i))
                    ub = mu[0]
                    for a in range(lb, ub + 1):
                        yield (a,) + nu
            else:
                yield ()

        for nu in _subpartitions(mu, strict, 1):
            yield cls.trim(nu)

    @classmethod
    def shifted_ribbon_complements(cls, mu):
        nu = Partition.trim(mu[1:])
        for kappa in cls.subpartitions(mu, nu, strict=True):
            yield kappa

    @classmethod
    def shifted_ribbon(cls, alpha):
        nu = []
        lam = []
        for i in range(len(alpha)):
            nu.append(sum(alpha[i:]))
            lam.append(sum(alpha[i + 1:]))
        return Partition.trim(nu), Partition.trim(lam)

    @classmethod
    def shifted_ribbons(cls, n):
        def generate(nu, lam, ub=None):
            if len(nu) == 0:
                yield nu, lam
            else:
                lb = lam[0]
                ub = nu[0] if ub is None else min(nu[0], ub)
                for a in range(lb, ub + 1):
                    for x, y in generate(nu[1:], lam[1:], a - 1):
                        yield [nu[0]] + x, [a] + y

        seen = set()
        for alpha in cls.compositions(n):
            nu, lam = cls.shifted_ribbon(alpha)
            if (tuple(nu), tuple(lam)) not in seen:
                for x, y in generate(nu, lam):
                    x = cls.trim(x)
                    y = cls.trim(y)
                    if (x, y) not in seen:
                        yield x, y
                        seen.add((x, y))

    @classmethod
    def decrement_one(cls, nu):
        i = 0
        while i + 1 < len(nu) and nu[i] - 1 <= nu[i + 1]:
            i += 1
        return cls.trim(nu[:i] + (nu[i] - 1,) + nu[i + 1:])

    @classmethod
    def decrement_strict(cls, nu):
        if nu == ():
            return nu
        i = 0
        mu = list(cls.trim(nu))
        mu[i] -= 1
        while i + 1 < len(mu) and 0 < mu[i] == mu[i + 1]:
            i += 1
            mu[i] -= 1
        return cls.trim(mu)

    @cached_value(DECOMPOSE_VSTRIP_CACHE)
    def decompose_shifted_shape_by_vertical_strips(cls, mu): # noqa
        mu = cls.trim(mu)
        if mu == ():
            return [()]

        i = 0
        while i < len(mu) and mu[i] + i == mu[0]:
            i += 1

        return [
            cls.trim(mu[:a] + tuple(mu[b] - 1 for b in range(a, i)) + nu)
            for a in range(i + 1)
            for nu in cls.decompose_shifted_shape_by_vertical_strips(mu[i:])
        ]

    @cached_value(DECOMPOSE_RIM_CACHE)
    def decompose_shifted_shape_by_rims(cls, mu): # noqa
        mu = cls.trim(mu)
        if mu == ():
            return [()]
        ans = []
        for nu in cls.decompose_shifted_shape_by_rims(mu[1:]):
            a = (nu[0] + 1) if nu else 0
            b = mu[1] if len(mu) > 1 else 0
            for a in range(max(a, b), mu[0] + 1):
                ans.append(cls.trim((a,) + nu))
        return ans

    @cached_value(RIM_CACHE)
    def _rims_helper(cls, border):  # noqa
        if border:
            ret = []
            (a, b), next_border = border[0], border[1:]
            # skip (a, b)
            while next_border and next_border[0][0] == a:
                next_border = next_border[1:]
            for ans in cls._rims_helper(next_border):
                ret.append(ans)
            # include (a, b)
            ans, next_border = [(a, b)], border[1:]
            while next_border and next_border[0][1] == b:
                ans.append(next_border[0])
                next_border = next_border[1:]
            for bns in cls._rims_helper(next_border):
                ret.append(ans + bns)
            return ret
        else:
            return [[]]

    @classmethod
    def rims(cls, mu, first_row_bound=1):
        assert cls.is_strict_partition(mu)
        mu = cls.trim(mu) + (0,)
        mu = (mu[0],) + mu
        border = []
        for i in range(len(mu) - 1, 0, -1):
            border += [(i, i + mu[i] + j) for j in range(mu[i - 1] - mu[i])]
        border += [(1, 1 + mu[0])]
        border = tuple(border)
        for r in cls._rims_helper(border):
            if border[-1] in r:
                for i in range(first_row_bound):
                    yield tuple(sorted(r + [(1, 1 + mu[0] + j) for j in range(1, i + 1)]))
            else:
                yield tuple(sorted(r))

    @classmethod
    def rim_terminal_part(cls, rim):
        t = tuple(sorted(rim, key=lambda x: (-x[0], x[1])))[-3:]
        if len(t) >= 2 and t[-2][0] != t[-1][0] and t[-2][1] != t[-1][1]:
            t = t[-1:]
        if len(t) == 3 and t[0][0] != t[1][0] and t[0][1] != t[1][1]:
            t = t[1:]
        return t

    @classmethod
    def add_box_to_row(cls, p, row):
        parts = [cls.get(p, i) + (1 if i == row else 0) for i in range(1, max(len(p), row) + 1)]
        return Partition.trim(parts)

    @classmethod
    def add_box_to_column(cls, p, col, shift):
        shape = cls.shifted_shape(p) if shift else cls.shape(p)
        a = [i for (i, j) in shape if j == col]
        row = max(a) + 1 if a else 1
        shape.add((row, col))
        if shift:
            assert all((i, j - 1) in shape or i == j for (i, j) in shape)
        else:
            assert all((i, j - 1) in shape or j == 1 for (i, j) in shape)
        parts = []
        for (i, j) in shape:
            while i - 1 >= len(parts):
                parts += [0]
            parts[i - 1] += 1
        return Partition.trim(parts)

    @classmethod
    def union(cls, p, q):
        parts = [max(cls.get(p, i), cls.get(q, i)) for i in range(1, max(len(p), len(q)) + 1)]
        return Partition.trim(parts)

    @classmethod
    def shifted_growth_diagram(cls, dictionary, m=None, n=None):
        def shdiff(nu, lam):
            return next(iter(cls.skew(nu, lam, shifted=True)))

        dictionary = {(a, i + 1) for i, a in enumerate(dictionary)} if type(dictionary) in [list, tuple] else dictionary
        dictionary = {k: 1 for k in dictionary} if type(dictionary) == set else dictionary
        n = max([0] + [i for i, _ in dictionary]) if n is None else n
        m = max([0] + [j for _, j in dictionary]) if m is None else m

        g = [[tuple() for _ in range(m + 1)] for _ in range(n + 1)]
        edges = [[False for _ in range(m + 1)] for _ in range(n + 1)]
        corners = [[None for _ in range(m + 1)] for _ in range(n + 1)]
        reasons = [['' for _ in range(m + 1)] for _ in range(n + 1)]

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                v = dictionary.get((i, j), 0)
                assert v in [0, 1]
                lam, nu, mu = g[i - 1][j - 1], g[i - 1][j], g[i][j - 1]
                if v == 1 and cls.get(lam, 1) == cls.get(mu, 1):
                    # case (1)
                    gamma = cls.add_box_to_row(mu, row=1)
                    reasons[i][j] = '01'
                elif v == 1:
                    # case (2)
                    assert cls.get(lam, 1) + 1 == cls.get(mu, 1)
                    gamma = mu
                    corners[i][j] = 1
                    reasons[i][j] = '02'
                elif mu == lam:
                    # case (3a)
                    gamma = nu
                    edges[i][j] = edges[i - 1][j]
                    corners[i][j] = corners[i - 1][j]
                    reasons[i][j] = '03a'
                elif nu == lam and not edges[i - 1][j] and corners[i - 1][j] is None:
                    # case (3b)
                    gamma = mu
                    reasons[i][j] = '03b'
                elif not cls.contains(mu, nu):
                    # case (4)
                    gamma = cls.union(nu, mu)
                    edges[i][j] = edges[i - 1][j]
                    reasons[i][j] = '04'
                elif lam != nu:
                    a, b = shdiff(nu, lam)
                    row = {(x, y) for (x, y) in cls.skew(mu, lam, shifted=True) if x == a + 1}
                    col = {(x, y) for (x, y) in cls.skew(mu, lam, shifted=True) if y == b + 1}

                    if not edges[i - 1][j] and a != b:
                        if len(row) == 0:
                            # case (5)
                            gamma = cls.add_box_to_row(mu, a + 1)
                            reasons[i][j] = '05'
                        else:
                            # case (6)
                            gamma = mu
                            corners[i][j] = a + 1
                            reasons[i][j] = '06'
                    else:
                        if len(col) == 0:
                            # case (7)
                            gamma = cls.add_box_to_column(mu, b + 1, shift=True)
                            edges[i][j] = True
                            reasons[i][j] = '07'
                        else:
                            # case (8)
                            gamma = mu
                            edges[i][j] = True
                            corners[i][j] = b + 1
                            reasons[i][j] = '08'
                elif lam == nu and not edges[i - 1][j]:
                    a = corners[i - 1][j]
                    b = cls.get(nu, a) + a - 1
                    skew = cls.skew(mu, lam, shifted=True)
                    if (a, b + 1) not in skew and (a + 1, b) not in skew:
                        # case (9) --
                        gamma = mu
                        corners[i][j] = a
                        reasons[i][j] = '09'
                    elif (a + 1, b) in skew:
                        # case (10)
                        gamma = mu
                        corners[i][j] = a + 1
                        reasons[i][j] = '10'
                    elif (a, b + 1) in skew:
                        if any(x == a + 1 for (x, y) in skew):
                            # case (??)
                            gamma = mu
                            corners[i][j] = a + 1
                            reasons[i][j] = '??a'
                        elif a == b:
                            # case (12)
                            gamma = mu
                            edges[i][j] = True
                            corners[i][j] = b + 1
                            reasons[i][j] = '12'
                        else:
                            # case (11)
                            gamma = cls.add_box_to_row(mu, a + 1)
                            reasons[i][j] = '11'
                    else:
                        raise Exception

                elif lam == nu and edges[i - 1][j]:
                    b = corners[i - 1][j]
                    a = max([x for (x, y) in cls.shifted_shape(nu) if y == b])
                    skew = cls.skew(mu, lam, shifted=True)

                    if (a, b + 1) not in skew and (a + 1, b) not in skew:
                        # case (13) --
                        gamma = mu
                        corners[i][j] = b
                        edges[i][j] = True
                        reasons[i][j] = '13'
                    elif (a, b + 1) in skew:
                        # case (14)
                        gamma = mu
                        corners[i][j] = b + 1
                        edges[i][j] = True
                        reasons[i][j] = '14'
                    elif (a + 1, b) in skew:
                        if any(y == b + 1 for (x, y) in skew):
                            # case (??)
                            gamma = mu
                            corners[i][j] = b + 1
                            edges[i][j] = True
                            reasons[i][j] = '??b'
                        else:
                            # case (15)
                            gamma = cls.add_box_to_column(mu, b + 1, shift=True)
                            edges[i][j] = True
                            reasons[i][j] = '15'
                    else:
                        raise Exception
                else:
                    raise Exception
                g[i][j] = gamma

        return g, edges, corners, reasons

    @classmethod
    def growth_diagram(cls, dictionary, m=None, n=None):
        dictionary = {k: 1 for k in dictionary} if type(dictionary) == set else dictionary
        n = max([0] + [i for i, _ in dictionary]) if n is None else n
        m = max([0] + [j for _, j in dictionary]) if m is None else m
        g = [[tuple() for _ in range(m + 1)] for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                carry, k = dictionary.get((i, j), 0), 1
                rho, mu, nu = g[i - 1][j - 1], g[i - 1][j], g[i][j - 1]
                lam = []
                while True:
                    v = max(cls.get(mu, k), cls.get(nu, k)) + carry
                    if v == 0:
                        break
                    lam.append(v)
                    carry, k = min(cls.get(mu, k), cls.get(nu, k)) - cls.get(rho, k), k + 1
                g[i][j] = Partition.trim(lam)
        return g

    @classmethod
    def print_growth_diagram(cls, g):
        g = cls.growth_diagram(g) if type(g) != list else g
        g = [[str(mu) for mu in row] for row in g]
        m = max([6] + [len(mu) for row in g for mu in row])
        g = [[mu + (m - len(mu)) * ' ' for mu in row] for row in g]
        print()
        for row in reversed(g):
            print(' '.join(row))
        print()
