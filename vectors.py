import math


class Vector:

    def __init__(self, dictionary={}, printer=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer

    @classmethod
    def base(cls, key, printer=None):
        return Vector({key: 1}, printer)

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

    def is_positive(self):
        if self.is_zero():
            return False
        return all(v > 0 for v in self.values())

    def __len__(self):
        return len(self.dictionary)

    def __eq__(self, other):
        if other == 0:
            return self.is_zero()
        else:
            # assert type(other) == Vector
            return (self - other).is_zero()

    def __iter__(self):
        return self.dictionary.__iter__()

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def __add__(self, other):
        if type(other) == Vector:
            keys = self.keys() | other.keys()
            return Vector({
                key: self[key] + other[key]
                for key in keys if self[key] + other[key]
            }, self.printer or other.printer)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        if type(other) == Vector:
            keys = self.keys() | other.keys()
            return Vector({
                key: self[key] - other[key]
                for key in keys if self[key] - other[key]
            }, self.printer or other.printer)
        else:
            return other.__rsub__(self)

    def __mul__(self, other):
        if self.is_scalar(other):
            return Vector({key: self[key] * other for key in self.keys()}, self.printer)
        elif type(other) == Vector:
            ans = Vector(printer=self.printer or other.printer)
            for a, x in self.items():
                for b, y in other.items():
                    ans += (x * y) * (a * b)
            return ans
        else:
            return self * Vector.base(other)

    def __rmul__(self, other):
        if self.is_scalar(other):
            return Vector({key: self[key] * other for key in self.keys()}, self.printer)
        elif type(other) == Vector:
            ans = Vector(printer=self.printer or other.printer)
            for a, x in other.items():
                for b, y in self.items():
                    ans += (x * y) * (a * b)
            return ans
        else:
            return self * Vector.base(other)

    def __neg__(self):
        return self.__mul__(-1)

    def __truediv__(self, c):
        assert type(c) == int
        assert c != 0
        assert all(v % c == 0 for v in self.values())
        return Vector(dictionary={k: v // c for k, v in self.items()}, printer=self.printer)

    def is_scalar(self, other):
        return type(other) == int

    def is_zero(self):
        return len(self) == 0

    def _repr_coeff(self, coeff):
        if coeff == 1:
            return ' + '
        elif coeff == -1:
            return ' - '
        elif coeff > 0:
            return ' + %s*' % coeff
        else:
            return ' - %s*' % -coeff

    def __repr__(self):
        printer = self.printer or repr
        sorted_items = sorted([(printer(key), value) for key, value in self.items()])
        base = ''.join(self._repr_coeff(value) + key for key, value in sorted_items)
        if base.startswith(' + '):
            return base[3:]
        elif base.startswith(' - '):
            return '-' + base[3:]
        else:
            return '0'

    def any(self):
        assert len(self.dictionary) > 0
        return next(iter(self.dictionary))

    @classmethod
    def is_linearly_independent_subset(cls, vectors):
        return len(cls.get_linearly_independent_subset(vectors)) == len(vectors)

    @classmethod
    def reduce_linearly_independent_subset(cls, vectors, progress=None):
        progress = [] if progress is None else progress
        for v in vectors[len(progress):]:
            for m, u in progress:
                if m is not None:
                    u_coeff, v_coeff = u[m], v[m]
                    if v_coeff != 0:
                        d = math.gcd(u_coeff, v_coeff)
                        v = d * v  - ((d * v_coeff) // u_coeff) * u
            m = v.any() if not v.is_zero() else None
            progress.append((m, v))
        return progress

    @classmethod
    def get_linearly_independent_subset(cls, vectors):
        ans = []
        progress = cls.reduce_linearly_independent_subset(vectors)
        for i in range(len(progress)):
            m, u = progress[i]
            if m is not None:
                ans.append(i)
        return ans
