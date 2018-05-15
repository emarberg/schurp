
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

    def __len__(self):
        return len(self.dictionary)

    def __eq__(self, other):
        assert type(other) == Vector
        return len((self - other).dictionary) == 0

    def __iter__(self):
        return self.dictionary.__iter__()

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def __add__(self, other):
        if type(other) == Vector:
            keys = self.keys() | other.keys()
            return Vector({key: self[key] + other[key] for key in keys}, self.printer or other.printer)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        if type(other) == Vector:
            keys = self.keys() | other.keys()
            return Vector({key: self[key] - other[key] for key in keys}, self.printer or other.printer)
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
        return self.__mul__(other)

    def is_scalar(self, other):
        return type(other) == int

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
