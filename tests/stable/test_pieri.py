from stable.tableaux import Tableau
from stable.partitions import Partition
from stable.vectors import Vector
import pytest
import stable.utils as utils
import time
from collections import defaultdict

LSHAPES = {}
RSHAPES = {}

LCACHE = {}
RCACHE = {}

SPLIT_CACHE = {}


class RowVector:

    def __init__(self, row=()):
        self.row = list(row)
        self.size = len(row)

    def __iter__(self):
        return iter(self.row)

    def __getitem__(self, i):
        return self.row[i]

    def __len__(self):
        return self.size

    def __add__(self, other):
        assert self.size == other.size
        return RowVector([self[i] + other[i] for i in range(self.size)])

    def __sub__(self, other):
        assert self.size == other.size
        return RowVector([self[i] - other[i] for i in range(self.size)])

    def __mul__(self, other):
        assert type(other) == int
        return RowVector([self[i] * other for i in range(self.size)])

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return RowVector([-i for i in self.row])

    def __rshift__(self, n):
        assert type(n) == int and n >= 0
        return RowVector(n * [0] + self.row[:self.size - n])

    def __eq__(self, other):
        if type(other) == int and other == 0:
            return all(self[i] == 0 for i in range(self.size))
        assert self.size == other.size
        return all(self[i] == other[i] for i in range(self.size))

    def append(self, x):
        self.row.append(x)
        self.size += 1

    def subvector(self, n):
        assert n <= self.size
        return RowVector(self.row[:n])

    def __hash__(self):
        return hash(tuple(self.row))

    def __repr__(self):
        m = max(1 + max([0] + [len(str(i)) for i in self.row]), 6)

        def pad(i):
            s = str(i)
            return (m - len(s)) * ' ' + s + ','

        return 'RowVector([ ' + (''.join([pad(i) for i in self.row])) + ' ])'

    def copy(self):
        return RowVector(self.row.copy())

    def __or__(self, other):
        return RowVector(self.row + other.row)

    @classmethod
    def linear_combination(cls, coefficients, vectors):
        assert len(vectors) > 0
        ans = coefficients[0] * vectors[0]
        for i in range(1, len(vectors)):
            ans += coefficients[i] * vectors[i]
        return ans

    @classmethod
    def print_matrix(cls, vectors):
        print()
        for v in vectors:
            print(str(v))
        print()

    @classmethod
    def zero(cls, length):
        row = length * [0]
        return cls(row)

    @classmethod
    def elementary(cls, index, length, coeff=1):
        row = length * [0]
        row[index] = coeff
        return cls(row)

    @classmethod
    def solve(cls, vectors, target, integral=False):
        coeffs = [-1, 1]
        n = len(vectors)
        i = [
            (a, c1)
            for a in range(n)
            for c1 in coeffs
            if target == vectors[a]]
        if i:
            a, c1 = i[0]
            return [cls.elementary(a, n, c1)]
        i = [
            (a, b, c1, c2)
            for a in range(n)
            for b in range(a + 1, n)
            for c1 in coeffs for c2 in coeffs
            if c1 * vectors[a] + c2 * vectors[b] == target]
        if i:
            a, b, c1, c2 = i[0]
            return [cls.elementary(a, n, c1) + cls.elementary(b, n, c2)]
        i = [
            (a, b, c, c1, c2, c3)
            for a in range(n)
            for b in range(a + 1, n)
            for c in range(b + 1, n)
            for c1 in coeffs for c2 in coeffs for c3 in coeffs
            if c1 * vectors[a] + c2 * vectors[b] + c3 * vectors[c] == target]
        if i:
            a, b, c, c1, c2, c3 = i[0]
            return [cls.elementary(a, n, c1) + cls.elementary(b, n, c2) + cls.elementary(c, n, c3)]
        i = [
            (a, b, c, d, c1, c2, c3, c4)
            for a in range(n)
            for b in range(a + 1, n)
            for c in range(b + 1, n)
            for d in range(c + 1, n)
            for c1 in coeffs for c2 in coeffs for c3 in coeffs for c4 in coeffs
            if c1 * vectors[a] + c2 * vectors[b] + c3 * vectors[c] + c4 * vectors[d] == target]
        if i:
            a, b, c, d, c1, c2, c3, c4 = i[0]
            return [cls.elementary(a, n, c1) + cls.elementary(b, n, c2) + cls.elementary(c, n, c3) + cls.elementary(d, n, c4)]
        i = [_ for _ in range(n) if target == -vectors[_]]
        if i:
            return [cls.elementary(i[0], n, -1)]

        def find_nonzero(arr, j):
            for i in range(j, len(arr[j])):
                if arr[j][i] != 0:
                    return i

        def cancel(arr, i, j, k):
            a = arr[j][i]
            b = arr[j][k]
            if b != 0:
                try:
                    import numpy
                    lcm = int(numpy.lcm(a, b))
                except:
                    return []
                for t in range(len(arr)):
                    arr[t].row[k] = arr[t][k] * (lcm // b) - arr[t][i] * (lcm // a)

        def kernel(arr):
            return []

        m, n = len(target), len(vectors)
        arr = [v.copy() for v in vectors] + [target.copy()]

        for j in range(n):
            i = find_nonzero(arr, j)
            if i is None:
                continue
            for k in range(m):
                if i != k:
                    cancel(arr, i, j, k)
            for v in range(n + 1):
                arr[v].row[j], arr[v].row[i] = arr[v][i], arr[v][j]

        solution = {}
        for i in range(m - 1, -1, -1):
            for j in range(n + 1):
                if arr[j][i] != 0:
                    if j == n:
                        return []
                    v = arr[n][i]
                    for k in range(j + 1, n):
                        v -= arr[k][i] * solution.get(k, 0)
                    if integral and v % arr[j][i] != 0:
                        return[]
                    solution[j] = v // arr[j][i] if v % arr[j][i] == 0 else v / arr[j][i]
                    break
        s = RowVector([solution.get(i, 0) for i in range(n)])
        return [s] + [s + vec for vec in kernel(arr[:-1])]


def split(nu, mu):
    if (nu, mu) not in SPLIT_CACHE:
        i, j = None, None
        for lam in lshapes(nu, mu):
            term = Partition.rim_terminal_part(left_rim(nu, mu, lam))
            i = max(([] if i is None else [i]) + [a for (a, _) in term])
            j = min(([] if j is None else [j]) + [b for (_, b) in term])
        for lam in rshapes(nu, mu):
            term = Partition.rim_terminal_part(right_rim(nu, mu, lam))
            i = max(([] if i is None else [i]) + [a for (a, _) in term])
            j = min(([] if j is None else [j]) + [b for (_, b) in term])
        skew = Partition.skew(nu, mu, shifted=True)
        # new_nu = Partition.trim(tuple(min(nu[t], j - t - 1) for t in range(i)) + nu[i:])
        # new_mu = Partition.trim(tuple(min(mu[t], j - t - 1) for t in range(i)) + mu[i:])
        while (i + 1, j) in skew:
            i += 1
        rest = {(a, b) for (a, b) in skew if (a <= i or b >= j)}
        split_nu, split_mu = Partition.from_skew(rest, shifted=False)
        diagonal = any(a == b for (a, b) in skew)
        SPLIT_CACHE[(nu, mu)] = (split_nu, split_mu, diagonal)
    return SPLIT_CACHE[(nu, mu)]


def left_rim(nu, mu, lam):
    return Partition.skew(nu, lam, shifted=True)


def right_rim(nu, mu, lam):
    return Partition.skew(lam, mu, shifted=True)


def print_lshape(nu, mu, lam):
    lines = []
    for i in range(1, len(nu) + 1):
        line = (i - 1) * [' '] + Partition.get(mu, i) * ['.']
        line += (Partition.get(lam, i) - Partition.get(mu, i)) * ['*']
        line += (Partition.get(nu, i) - Partition.get(lam, i)) * ['R']
        lines += [' '.join(line)]
    print('\n'.join(reversed(lines)) + '\n')
    print('term:', Partition.rim_terminal_part(left_rim(nu, mu, lam)))
    print()


def print_rshape(nu, mu, lam):
    lines = []
    for i in range(1, len(nu) + 1):
        line = (i - 1) * [' '] + Partition.get(mu, i) * ['.']
        line += (Partition.get(lam, i) - Partition.get(mu, i)) * ['R']
        line += (Partition.get(nu, i) - Partition.get(lam, i)) * ['*']
        lines += [' '.join(line)]
    print('\n'.join(reversed(lines)) + '\n')
    print('term:', Partition.rim_terminal_part(right_rim(nu, mu, lam)))
    print()


def columns(lam, mu):
    return len({j for (i, j) in Partition.skew(lam, mu, True)})


def lvector(nu, mu, p_max):
    ans = LCACHE.get((nu, mu), RowVector())
    if len(ans) >= p_max:
        return ans.subvector(p_max)
    for p in range(len(ans), p_max):
        v = 0
        for lam in lshapes(nu, mu):
            c = sum(lam) - sum(mu)

            if p == 0:
                d = Tableau.KOG_counts(nu, lam, p)
            else:
                d = 2 * Tableau.KOG_counts(nu, lam, p) + Tableau.KOG_counts(nu, lam, p + 1)
            v += d * 2 ** (len(mu) - c) * (-1) ** (columns(lam, mu) + c)

        if p == 0 and columns(nu, mu) == 1:
            v *= -1

        ans.append(v)
    LCACHE[(nu, mu)] = ans
    return ans


def rvector(nu, mu, p_max):
    ans = RCACHE.get((nu, mu), RowVector())
    if len(ans) >= p_max:
        return ans.subvector(p_max)
    for p in range(len(ans), p_max):
        v = 0
        for lam in rshapes(nu, mu):
            c = sum(nu) - sum(lam)
            v += 2 ** (len(lam) - c) * (-1) ** (columns(nu, lam) + c) * Tableau.KLG_counts(lam, mu, p)

        if p == 0 and columns(nu, mu) == 1:
            v *= -1

        ans.append(v)
    RCACHE[(nu, mu)] = ans
    return ans


def lshapes(nu, mu=None):
    if nu not in LSHAPES:
        LSHAPES[nu] = {}
        for lam in Partition.decompose_shifted_shape_by_rims(nu):
            for x in Partition.decompose_shifted_shape_by_vertical_strips(lam):
                if len(x) == len(lam):
                    if x not in LSHAPES[nu]:
                        LSHAPES[nu][x] = []
                    LSHAPES[nu][x].append(lam)
    if mu is None:
        return LSHAPES[nu]
    else:
        return LSHAPES[nu].get(mu, [])


def rshapes(nu, mu=None):
    if nu not in RSHAPES:
        RSHAPES[nu] = {}
        for lam in Partition.decompose_shifted_shape_by_vertical_strips(nu):
            if len(lam) == len(nu):
                for x in Partition.decompose_shifted_shape_by_rims(lam):
                    if x not in RSHAPES[nu]:
                        RSHAPES[nu][x] = []
                    RSHAPES[nu][x].append(lam)
    if mu is None:
        return RSHAPES[nu]
    else:
        return RSHAPES[nu].get(mu, [])


def is_reducible(nu, mu):
    if len(nu) == 0:
        return True
    if any(nu[i] == mu[i] for i in range(len(mu))):
        return True
    if sum(nu) - sum(mu) <= 3:
        return True
#    if any(nu[i] + 1 < mu[i - 1] for i in range(1, len(nu))):
#        return True
#    if nu == mu:
#        return True
    # if Partition.get(nu, 1) <= Partition.get(mu, 1) + 1:
    #    return True

    skew = Partition.skew(nu, mu, shifted=True)
    (i, j) = min(skew, key=lambda ij: (-ij[0], ij[1]))
    while True:
        if (i - 1, j) in skew:
            (i, j) = (i - 1, j)
        elif (i, j + 1) in skew:
            (i, j) = (i, j + 1)
        else:
            break
    if (i, j) != max(skew, key=lambda ij: (-ij[0], ij[1])):
        return True


def deep_dive(nu=(6, 5, 3, 2, 1), mu=(4, 3, 2, 1)):
    split_nu, split_mu, diagonal = split(nu, mu)
    print()
    print('TYPE:')
    print()
    print(split_nu, split_mu, diagonal)
    print()
    print(Partition.printable(split_nu, split_mu, shifted=False))
    print()

    print('nu =', nu, 'mu =', mu)
    print()
    print(Partition.printable(nu, mu, shifted=True))
    print()

    m = sum(nu) - sum(mu) + 2
    target = lvector(nu, mu, m)

    nu1, mu1 = Partition.decrement_one(nu), mu
    nu2, mu2 = Partition.decrement_one(Partition.decrement_one(nu)), mu
    nu3, mu3 = Partition.decrement_one(Partition.decrement_one(Partition.decrement_one(nu))), Partition.decrement_one(mu)
    nu4, mu4 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
    nu5, mu5 = Partition.trim(nu[2:]), Partition.trim(mu[2:])
    nu6, mu6 = tuple(i - 1 for i in nu), tuple(i - 1 for i in mu)

    print(Partition.printable(nu1, mu1, shifted=True), '\n')
    print(Partition.printable(nu2, mu2, shifted=True), '\n')
    print(Partition.printable(nu3, mu3, shifted=True), '\n')
    print(Partition.printable(nu4, mu4, shifted=True), '\n')
    print(Partition.printable(nu5, mu5, shifted=True), '\n')

    m = sum(nu) - sum(mu) + 2

    print('nu =', nu, '; mu =', mu)
    print()
    print(Partition.printable(nu, mu, shifted=True))
    print()

    v = [
        lvector(nu1, mu1, m),
        lvector(nu1, mu1, m) >> 1,
        lvector(nu2, mu2, m),
        lvector(nu2, mu2, m) >> 1,
        lvector(nu3, mu3, m),
        lvector(nu3, mu3, m) >> 1,
        lvector(nu4, mu4, m),
        lvector(nu4, mu4, m) >> 1,
        lvector(nu5, mu5, m),
        lvector(nu5, mu5, m) >> 1,
        lvector(nu5, mu5, m) >> 2,
    ]
    print('target =')
    print(target)
    print()
    RowVector.print_matrix(v)
    print()
    print('partial solutions =', RowVector.solve(v, target, integral=False))
    print()
    # for mu2 in Partition.subpartitions(mu, strict=True):
    #     for nu2 in Partition.subpartitions(nu, strict=True):
    #         if (nu, mu) != (nu2, mu2) and Partition.contains(nu2, mu2):
    #             if not Partition.skew(nu2, mu2).issubset(Partition.skew(nu, mu)):
    #                 continue
    #             v = lvector(nu2, mu2, m)
    #             w = lvector(nu2, mu2, m) >> 1
    #             print(v, nu2, mu2)
    #             print(w)
    #             print()
    #             print(Partition.printable(nu2, mu2, shifted=True))
    #             print()
    #             vectors += [v, w]


def get_types(n):
    types = defaultdict(list)
    for nu in Partition.all(n, strict=True):
        for mu in set(lshapes(nu)) | set(rshapes(nu)):
            if is_reducible(nu, mu):
                continue
            triple = split(nu, mu)
            types[triple].append((nu, mu))
    return types


def test_all(n=6):
    def inspect(w, target, good, unique=None):
        ans = []
        sort = sorted(good, key=lambda x: -unique[x]) if unique else good
        for g in sort:
            a = w[0] * g[0]
            for i in range(1, len(g)):
                a += w[i] * g[i]
            if target == a:
                ans.append(g)
        return ans

    types = get_types(n)
    success = [0, 0]
    unique = {}
    refined_unique = defaultdict(list)
    good = set()

    for (split_nu, split_mu, diagonal) in types:
        # if (split_nu, split_mu) != ((3, 1, 1, 1), ()):
        #    continue
        target = RowVector()
        w = [RowVector() for i in range(9)]

        # print(10 * '\n')
        domain = types[(split_nu, split_mu, diagonal)]

        # if len(domain) <= 2:
        #     continue

        for (nu, mu) in domain:
            nu1, mu1 = Partition.decrement_one(nu), mu
            # nu2, mu2 = Partition.decrement_one(Partition.decrement_one(nu)), mu
            # nu3, mu3 = Partition.decrement_one(Partition.decrement_one(nu)), Partition.decrement_one(mu)
            nu4, mu4 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
            nu5, mu5 = Partition.trim(nu[2:]), Partition.trim(mu[2:])
            nu6, mu6 = Partition.trim(nu[3:]), Partition.trim(mu[3:])

            m = sum(nu) - sum(mu) + 2

            v = lvector(nu, mu, m)
            target |= v

            w[0] |= lvector(nu4, mu4, m)
            w[1] |= lvector(nu4, mu4, m) >> 1
            w[2] |= lvector(nu4, mu4, m) >> 2

            w[3] |= lvector(nu5, mu5, m)
            w[4] |= lvector(nu5, mu5, m) >> 1
            w[5] |= lvector(nu5, mu5, m) >> 2

            w[6] |= lvector(nu1, mu1, m)
            w[7] |= lvector(nu1, mu1, m) >> 1
            w[8] |= lvector(nu1, mu1, m) >> 2

            # v = RowVector.elementary(1, len(v))
            # w[9] |= v

            # w[9] |= lvector(nu6, mu6, m)
            # w[10] |= lvector(nu6, mu6, m) >> 1
            # w[11] |= lvector(nu6, mu6, m) >> 2

            # w[8] |= lvector(nu2, mu2, m)
            # w[9] |= lvector(nu2, mu2, m) >> 1

            # w[10] |= lvector(nu3, mu3, m)
            # w[11] |= lvector(nu3, mu3, m) >> 1

            # w[12] |= lvector(nu6, mu6, m)
            # w[13] |= lvector(nu6, mu6, m) >> 1

            try:
                assert lvector(nu, mu, m) == rvector(nu, mu, m)
            except:
                print(Partition.printable(nu, mu, shifted=True), '\n')
                print(lvector(nu, mu, m), '==', rvector(nu, mu, m))
                input('')

        k = 1
        solved = inspect(w, k * target, good, unique)
        if not solved:
            solved = [k * x for x in RowVector.solve(w, target)]
            if solved:
                integral = [RowVector([int(a) for a in solved[0]])]
                if inspect(w, k * target, integral):
                    solved = integral

        if solved:  # and all(type(t) == int for t in solved[0]):
            solved = solved[0]
            success[0] += 1
            unique[solved] = unique.get(solved, 0) + 1
            for (nu, mu) in domain:
                x = sum(nu) - sum(mu)
            refined_unique[solved].append((split_nu, split_mu, diagonal, domain, w, target))
            # if len(domain) > 1:
            #    print(len(domain), ':', solved, diagonal)
            if all(type(i) == int for i in solved):
                good.add(solved)
            if solved[-1] == 0:
                continue
        else:
            success[1] += 1
            for (nu, mu) in domain:
                # if mu == ():
                #    continue

                nu1, mu1 = Partition.decrement_one(nu), mu
                # nu2, mu2 = Partition.decrement_one(Partition.decrement_one(nu)), mu
                # nu3, mu3 = Partition.decrement_one(Partition.decrement_one(nu)), Partition.decrement_one(mu)
                nu4, mu4 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
                nu5, mu5 = Partition.trim(nu[2:]), Partition.trim(mu[2:])
                # nu6, mu6 = Partition.decrement_one(Partition.decrement_one(nu)), mu

                m = sum(nu) - sum(mu) + 2

                print('nu =', nu, '; mu =', mu)
                print()
                print(Partition.printable(nu, mu, shifted=True))
                print()

                # print(Partition.printable(nu1, mu1, shifted=True), '\n')
                # print(Partition.printable(nu2, mu2, shifted=True), '\n')
                # print(Partition.printable(nu3, mu3, shifted=True), '\n')
                # print(Partition.printable(nu4, mu4, shifted=True), '\n')
                # print(Partition.printable(nu5, mu5, shifted=True), '\n')
                # print(Partition.printable(nu6, mu6, shifted=True), '\n')

                v = [
                    lvector(nu4, mu4, m),
                    lvector(nu4, mu4, m) >> 1,
                    lvector(nu4, mu4, m) >> 2,
                    lvector(nu5, mu5, m),
                    lvector(nu5, mu5, m) >> 1,
                    lvector(nu5, mu5, m) >> 2,
                    lvector(nu1, mu1, m),
                    lvector(nu1, mu1, m) >> 1,
                    lvector(nu1, mu1, m) >> 2,
                    # lvector(nu2, mu2, m),
                    # lvector(nu2, mu2, m) >> 1,
                    # lvector(nu3, mu3, m),
                    # lvector(nu3, mu3, m) >> 1,
                    # lvector(nu6, mu6, m),
                    # lvector(nu6, mu6, m) >> 1,
                ]
                print('target =')
                print(lvector(nu, mu, m))
                print()
                RowVector.print_matrix(v)
                print()
                print('partial solutions =', RowVector.solve(v, lvector(nu, mu, m), integral=False))
                print()
                # print('solutions so far:')
                # for u in unique:
                #     print('*', u)

        RowVector.print_matrix(w)
        RowVector.print_matrix([target])

        print(Partition.printable(nu, mu, shifted=True), '\n')
        print(Partition.printable(nu1, mu1, shifted=True), '\n')

        print()
        print('TYPE:')
        print()
        print(split_nu, split_mu, diagonal)
        print()
        print(Partition.printable(split_nu, split_mu, shifted=False))
        print()
        print()
        print(solved)
        print()
        # input('')
        print(10 * '\n')
    print()

    for u in sorted(unique, key=lambda x: -unique[x]):
        print(u, ':', unique[u])
        for (split_nu, split_mu, diagonal, domain, w, target) in refined_unique[u][:10]:
            RowVector.print_matrix(w)
            RowVector.print_matrix([target])
            print()
            print('type:', split_nu, split_mu, diagonal)
            print()
            print(Partition.printable(split_nu, split_mu, shifted=False))
            print()
            # i = 1
            # for nu, mu in domain:
            #     print(i, 'of', len(domain), '\n')
            #     print(Partition.printable(nu, mu, shifted=True))
            #     print()
            #     i += 1

    print()
    print('{')
    for u in sorted(unique, key=lambda x: unique[x]):
        print(u, ':', unique[u], ',')  # , unique[u], set(refined_unique[u]))
    print('}')
    print()
    print('unique:', len(unique))
    print()
    print(success)


def summarize_tex(n=6):
    def inspect(w, target, good, unique=None):
        ans = []
        sort = good
        for solution in sort:
            g, scalar = solution
            a = w[0] * g[0]
            for i in range(1, len(g)):
                a += w[i] * g[i]
            if scalar * target == a:
                ans.append(solution)
        return ans

    def findsolutions(w, target, good, unique):
        solved = inspect(w, target, good, unique)
        if solved:
            return solved[0]
        solved = RowVector.solve(w, target)
        for c in [1, 2, 3]:
            integral = [(RowVector([int(c * a) for a in s]), c) for s in solved]
            integral = inspect(w, target, integral)
            if integral:
                return integral[0]
        return solved[0], 1

    def get_target(domain):
        target = RowVector()
        w = [RowVector() for i in range(9)]
        for (nu, mu) in domain:
            nu1, mu1 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
            nu2, mu2 = Partition.trim(nu[2:]), Partition.trim(mu[2:])
            nu3, mu3 = Partition.decrement_one(nu), mu

            m = sum(nu) - sum(mu) + 2
            v = lvector(nu, mu, m)
            target |= v

            w[0] |= lvector(nu1, mu1, m)
            w[1] |= lvector(nu1, mu1, m) >> 1
            w[2] |= lvector(nu1, mu1, m) >> 2

            w[3] |= lvector(nu2, mu2, m)
            w[4] |= lvector(nu2, mu2, m) >> 1
            w[5] |= lvector(nu2, mu2, m) >> 2

            w[6] |= lvector(nu3, mu3, m)
            w[7] |= lvector(nu3, mu3, m) >> 1
            w[8] |= lvector(nu3, mu3, m) >> 2

            assert lvector(nu, mu, m) == rvector(nu, mu, m)
        return target, w

    good = [
        # worked out cases
        (RowVector([0, 0, 0, 0, 0, 0, 0, 1, 0, ]), 1),
        (RowVector([0, 0, 0, 0, 0, 0, 1, 1, 0, ]), 1),
        (RowVector([0, -1, 0, 0, 0, 0, 0, 0, 0, ]), 1),
        (RowVector([1, 4, 2, 0, 0, 0, 0, 0, 0, ]), 1),
        # manually added
        (RowVector([0, 0, 0, 3, 2, 0, 1, 1, 0, ]), 1),
        (RowVector([0, 0, 0, 1, 0, 0, 0, 1, 0, ]), 1),
        #
        (RowVector([0, 2, 0, 0, 0, 0, 0, 0, 0, ]), 1),
        (RowVector([0, 1, 0, 0, 1, 0, 0, 0, 0, ]), 1),
        (RowVector([0, 1, 0, 0, 0, 0, 0, 1, 0, ]), 1),
        (RowVector([1, 1, 0, 0, 0, 0, 0, 0, 0, ]), 1),
        (RowVector([0, 1, 0, 0, 0, 0, -1, 0, 0, ]), 1),
        (RowVector([0, 2, 0, 1, 0, 0, 0, 0, 0, ]), 1),
        (RowVector([3, 2, 0, 0, 0, 0, 0, 0, 0, ]), 1),
        (RowVector([1, 1, 0, 0, 0, 0, 0, 1, 0, ]), 1),
        (RowVector([0, 1, 0, 0, 0, 0, 1, 1, 0, ]), 1),
        (RowVector([3, 2, 0, 0, -2, 0, 0, 0, 0, ]), 1),
        (RowVector([-1, 1, 0, 0, 1, 0, 0, 1, 0, ]), 1),
        (RowVector([0, 0, 0, 1, 1, 0, 1, 1, 0, ]), 1),
        (RowVector([1, 1, 0, 0, 0, 0, 1, 1, 0, ]), 1),
        (RowVector([1, 1, 0, 1, 0, 0, 0, 1, 0, ]), 1),
        (RowVector([0, 1, 0, 1, 0, 0, 1, 1, 0, ]), 1),
        (RowVector([-1, 1, 0, 3, 3, 0, 0, 0, 0, ]), 1),
        (RowVector([3, 2, 0, -2, -2, 0, 0, 0, 0, ]), 1),
    ]

    types = get_types(n)
    refined_unique = defaultdict(list)
    unique = {}
    for (split_nu, split_mu, diagonal) in types:
        domain = types[(split_nu, split_mu, diagonal)]
        target, w = get_target(domain)
        solved = findsolutions(w, target, good, unique)
        unique[solved] = unique.get(solved, 0) + 1
        refined_unique[solved].append((split_nu, split_mu, diagonal, domain, w, target))
        if all(type(i) == int for i in solved[0]) and solved not in good:
            good.append(solved)
        if solved[1] > 1:
            print('*', (split_nu, split_mu, diagonal))
            print(' ', solved)
            print()
            RowVector.print_matrix(w)
            RowVector.print_matrix([target])
            # input('?')

    def sorter(x):
        index = ([i for i in range(len(good)) if x == good[i]] + [len(good)])[0]
        row, scalar = x
        a = len([i for i in row if i != 0])
        b = abs(scalar) + sum([abs(i) for i in row])
        boolean = any(type(i) != int for i in row)
        return index, boolean, scalar, a, b

    print('{')
    for u in sorted(unique, key=sorter):
        print(u, ':', unique[u], ',')
    print('}')
    print()
    print('unique:', len(unique))
    print()

    tokens = [
        "\\nu', \\mu', p",
        "\\nu', \\mu', p - 1",
        "\\nu', \\mu', p - 2",
        "\\nu'', \\mu'', p",
        "\\nu'', \\mu'', p - 1",
        "\\nu'', \\mu'', p - 2",
        "\\nu^-, \\mu, p",
        "\\nu^-, \\mu, p - 1",
        "\\nu^-, \\mu, p - 2",
        '\\nu, \\mu, p',
    ]
    tokens = ['\\Sigma(' + s + ')' for s in tokens]

    def printrec(rowvec, scalar):
        ans = ''
        lhs = (str(scalar) if scalar != 1 else '') + tokens[-1]
        for i, c in enumerate(rowvec[:-1]):
            if c != 0:
                ans += ('+ ' if c > 0 else '- ') + (str(abs(c)) if abs(c) != 1 else '') + ' ' + tokens[i]
        ans = ans[1:] if ans[:1] == '+' else ans
        return '\\begin{enumerate}\\item[] $' + lhs + ' = ' + ans + '$ \\end{enumerate}\n\n'

    def printtype(split_nu, split_mu, diagonal):
        lines = []
        for i in range(1, len(split_nu) + 1):
            b = Partition.get(split_nu, i)
            a = Partition.get(split_mu, i)
            line = ' & '.join(['\\none[\\cdot]' if i < a else '\\ ' for i in range(b)])
            lines.append(line)
        ytab = '$\\ytabsmall{\n' + ' \\\\\n'.join(reversed(lines)) + '\n}\\quad$'

        ans = '\\begin{itemize}\\item '
        ans += ytab
        # ans += '$\\mu = %s' % str(split_mu) + '$, '
        # ans += '$\\nu = %s' % str(split_nu) + '$, '
        ans += 'meets diagonal: {\\tt ' + str(bool(diagonal)) + '}.'
        ans += ' \\end{itemize}\n'
        return ans

    directory = '/Users/emarberg/Dropbox/projects/shifted-stable-grothendieck/tex/2 - hopf algebras/draft1/'
    filename = directory + '_recurrences_.tex'
    with open(filename, 'w') as f:
        f.write("""
The following recurrences are observed for pairs of strict partititions
$\\mu \\subseteq \\nu \\vdash n$ for $n \\leq %s$.
I'm only including recurrences that are observed more than once, as there is
a long tail of probability undetermined cases. Some relevant notation:
\\begin{itemize}
\\item If $\\lambda$ is a partition, then form $\\lambda'$ by omitting the first part.
\\item If $\\lambda$ is a partition, then form $\\lambda'' = (\\lambda')'$ by omitting the first and second part.
\\item If $\\lambda$ is a nonempty strict partition, then let $\\lambda^-$ be the
strict partition whose shifted Young diagram is formed from that of $\\lambda$ by
removing the last box in the last column.
\\item Let $\\Sigma(\\nu, \\mu, p)$ be either side of \\eqref{eq:desired identity}.
\\end{itemize}
Here are the recurrences:
\\begin{enumerate}""" % str(n))

        index = 1
        for row, scalar in sorted(unique, key=sorter):
            count = unique[(row, scalar)]
            if count <= 1:
                continue
            f.write('\\item[(%s)] Observed %s times:\n\n' % (index, count))
            f.write(printrec(row, scalar))
            index += 1

        f.write('\\end{enumerate}')

        index = 1
        for key in sorted(unique, key=sorter):
            row, scalar = key
            count = unique[key]
            if count <= 1:
                continue
            f.write('\n\\newpage\n\\subsection{Recurrence %s}\n Observed %s times:\n\n' % (index, count))
            f.write(printrec(row, scalar))
            f.write('\\noindent Cases that exhibit this recurrence:')
            for (split_nu, split_mu, diagonal, domain, w, target) in sorted(refined_unique[key]):
                f.write(printtype(split_nu, split_mu, diagonal))
            index += 1


def seconds(t):
    return str(int(1000 * t) / 1000.0)


# def pairs_io(n=5):
#     assert 0 <= n < 256

#     def tobytes(nu, lam, mu):
#         return bytes(nu + (0,) + lam + (0,) + mu + (0,))

#     def frombytes(b, i):
#         nu = []
#         while i < len(b) and b[i] != 0:
#             nu.append(b[i])
#             i += 1
#         i += 1
#         lam = []
#         while i < len(b) and b[i] != 0:
#             lam.append(b[i])
#             i += 1
#         i += 1
#         mu = []
#         while i < len(b) and b[i] != 0:
#             mu.append(b[i])
#             i += 1
#         i += 1
#         return tuple(nu), tuple(lam), tuple(mu), i

#     def read(filename):
#         print()
#         print('Trying to read from `%s`' % filename)
#         t0 = time.time()
#         with open(filename, 'rb') as file:
#             b = bytearray(file.read())
#         print('* Opened file in %s seconds' % seconds(time.time() - t0))
#         ans = defaultdict(list)
#         i = 0
#         while i < len(b):
#             nu, lam, mu, i = frombytes(b, i)
#             ans[(nu, mu)].append(lam)
#         print('* Succeeded in %s seconds' % seconds(time.time() - t0))
#         print()
#         return ans

#     def write(file, ans):
#         t0 = time.time()
#         b = bytearray()
#         for nu, mu in ans:
#             for lam in ans[(nu, mu)]:
#                 b += tobytes(nu, lam, mu)
#         print('* Writing to file `%s`' % file)
#         with open(file, 'wb') as f:
#             f.write(b)
#         print('* Succeeded in %s seconds' % seconds(time.time() - t0))
#         print()

#     directory = '/Users/emarberg/examples/gpq/'
#     lfile = directory + 'rims_then_strips_%s.b' % n
#     rfile = directory + 'strips_then_rims_%s.b' % n

#     try:
#         lhs = read(lfile)
#     except FileNotFoundError:
#         print('* Failed, computing instead')
#         lhs = pairs_lhs(n)
#         write(lfile, lhs)

#     try:
#         rhs = read(rfile)
#     except FileNotFoundError:
#         print('* Failed, computing instead')
#         rhs = pairs_rhs(n)
#         write(rfile, rhs)

#     return lhs, rhs


def test_rims():
    assert set(Partition.rims((), 0)) == {()}
    assert set(Partition.rims((), 1)) == {(), ((1, 1),)}
    assert set(Partition.rims((), 2)) == {(), ((1, 1),), ((1, 1), (1, 2))}
    assert set(Partition.rims((1,), 1)) == {(), ((1, 2),), ((1, 2), (2, 2))}
    assert set(Partition.rims((3, 2))) == {(), ((1, 4),), ((1, 4), (2, 4)), ((3, 3),), ((1, 4), (3, 3)), ((1, 4), (2, 4), (3, 3)),  ((1, 4), (2, 4), (3, 3), (3, 4))} # noqa


def test_KOG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KOG(0, mu)) == {Tableau()}

    assert set(Tableau.KOG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KOG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): 1,
            (2, 6): 1,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 3,
            (4, 4): 3,
        })
    }


def test_KLG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KLG(0, mu)) == {Tableau()}

    assert set(Tableau.KLG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KLG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): -1,
            (2, 6): -1,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 3,
        })
    }


def test_fast_KOG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KOG_counts_by_shape(p, mu)
            tabs = Tableau.KOG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                assert counts[nu] == len(tabs[nu])


def test_fast_KLG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KLG_counts_by_shape(p, mu)
            tabs = Tableau.KLG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                print('*', nu)
                assert counts[nu] == len(tabs[nu])


def test_GP_pieri(): # noqa
    for mu in Partition.all(6, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


def test_GQ_pieri(): # noqa
    for mu in Partition.all(5, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)


@pytest.mark.slow
def test_GP_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


@pytest.mark.slow
def test_GQ_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)
