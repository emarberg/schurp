from .polynomials import Polynomial


class Vector:

    def __init__(self, dictionary={}, printer=None, multiplier=None, sorter=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer
        self.sorter = sorter
        self.multiplier = multiplier

    def __hash__(self):
        return hash(str(self))

    @classmethod
    def print_matrix(cls, matrix):
        print()
        w = max([len(str(e)) for row in matrix for e in row])
        for row in matrix:
            print('  ', '[' + ', '.join(map(lambda e: (w - len(str(e))) * ' ' + str(e), row)) + ']')
        print()

    def is_expressable(self, *args):
        keys = set(self.dictionary)
        for a in args:
            keys |= set(a.dictionary)
        keylist = sorted(keys)
        # keydict = {k: i for i, k in enumerate(keylist)}
        m, n = len(keylist), len(args) + 1
        matrix = [[args[j][keylist[i]] for j in range(n - 1)] + [self[keylist[i]]] for i in range(m)]
        rref = self.rref(matrix)
        # self.print_matrix(rref)

        pivots = []
        for i in range(m):
            for j in range(n):
                if rref[i][j] != 0:
                    pivots.append(j)
                    break

        return n - 1 not in pivots

    @classmethod
    def rref(cls, matrix):
        assert all(type(e) == int for row in matrix for e in row)
        m = len(matrix)
        n = 0 if m == 0 else len(matrix[0])
        matrix = [[e for e in row] for row in matrix]

        def gcd(*args):
            if len(args) <= 1:
                return args[0]
            if len(args) > 2:
                return gcd(args[0], gcd(*args[1:]))
            a, b = args[0], args[1]
            while b != 0:
                a, b = b, a % b
            return a

        def swap(mat, i, j):
            mat[i], mat[j] = mat[j], mat[i]

        def scale(mat, i, v):
            for j in range(len(mat[i])):
                mat[i][j] *= v

        def iscale(mat, i, v=None):
            v = gcd(*mat[i]) if v is None else v
            for j in range(len(mat[i])):
                assert mat[i][j] % v == 0
                mat[i][j] //= v

        def replace(mat, i, j, v):
            for t in range(len(mat[j])):
                mat[j][t] += v * mat[i][t]

        def find_nonzeros_in_column(mat, rowstart, j):
            return [t for t in range(rowstart, m) if mat[t][j] != 0]

        row = 0
        for col in range(n):
            nonzeros = find_nonzeros_in_column(matrix, row, col)
            if len(nonzeros) == 0:
                continue

            i = nonzeros[0]
            iscale(matrix, i)
            if matrix[i][col] < 0:
                iscale(matrix, i, -1)

            for j in nonzeros[1:]:
                a, b = matrix[i][col], matrix[j][col]
                d = gcd(a, b)
                scale(matrix, j, a // d)
                replace(matrix, i, j, -b // d)

            if i != row:
                swap(matrix, row, i)

            row += 1

        return matrix

    @classmethod
    def base(cls, key, printer=None, multiplier=None):
        return cls({key: 1}, printer, multiplier)

    def keys(self):
        return self.dictionary.keys()

    def values(self):
        return self.dictionary.values()

    def items(self):
        return self.dictionary.items()

    def is_singleton(self):
        if len(self.dictionary) != 1:
            return False
        return list(self.dictionary.values())[0] == 1

    def is_nonnegative(self):
        return all(v > 0 for v in self.values())

    def is_positive(self):
        return not self.is_zero() and self.is_nonnegative()

    def __len__(self):
        return len(self.dictionary)

    def __eq__(self, other):
        return len((self - other).dictionary) == 0

    def __iter__(self):
        return self.dictionary.__iter__()

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def _instantiate(self, dictionary, other=None):
        other = other or self
        return self.__class__(
            dictionary,
            self.printer or other.printer,
            self.multiplier or other.multiplier,
            self.sorter or other.sorter,
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __add__(self, other):
        if other == 0:
            return self
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self._instantiate({key: self[key] + other[key] for key in keys}, other)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        if other == 0:
            return self
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self._instantiate({key: self[key] - other[key] for key in keys}, other)
        else:
            return other.__rsub__(self)

    def __mul__(self, other):
        if type(other) in [int, Polynomial]:
            return self._instantiate({key: self[key] * other for key in self.keys()})
        elif type(other) == type(self):
            ans = {}
            for a, x in self.items():
                for b, y in other.items():
                    for m, coeff in self.multiplier(a, b) if self.multiplier else [(a * b, 1)]:
                        ans[m] = ans.get(m, 0) + x * y * coeff
            return self._instantiate(ans, other)
        else:
            return self * self.base(other)

    def set_beta(self, v):
        return self._instantiate({key: value.set(0, v) if type(value) == Polynomial else val for (key, value) in self.items()})

    def __floordiv__(self, other):
        assert type(other) == int
        ans = self._instantiate({key: value // other for (key, value) in self.items()})
        assert ans * other == self
        return ans

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        return self.__mul__(-1)

    def __pow__(self, other):
        assert type(other) == int and other > 0
        if other == 1:
            return self
        x = self ** (other // 2)
        if other % 2 == 0:
            return x * x
        else:
            return x * x * self

    def is_zero(self):
        return len(self) == 0

    def _repr_coeff(self, coeff, key):
        if type(coeff) == Polynomial and len(coeff) > 1:
            return ' + (%s)*' % str(coeff)
        try:
            if key == '' and coeff > 0:
                return ' + %s' % coeff
            if key == '' and coeff < 0:
                return ' - %s' % -coeff
            if coeff == 1:
                return ' + '
            elif coeff == -1:
                return ' - '
            elif coeff > 0:
                return ' + %s*' % coeff
            else:
                return ' - %s*' % -coeff
        except TypeError:
            return ' + (%s)*' % str(coeff)

    def _print_sorted(self, sorted_items):
        base = ''.join(self._repr_coeff(value, key) + key for key, value in sorted_items)
        if base.startswith(' + '):
            return base[3:]
        elif base.startswith(' - '):
            return '-' + base[3:]
        else:
            return '0'

    def __repr__(self):
        printer = self.printer or repr
        sorter = lambda kv: ((self.sorter or printer)(kv[0]), kv[1])
        sorted_items = sorted([(key, value) for key, value in self.items()], key=sorter)
        sorted_items = [(printer(key), value) for key, value in sorted_items]
        return self._print_sorted(sorted_items)
