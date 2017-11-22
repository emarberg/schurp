
class Vector:
    def __init__(self, dictionary={}, printer=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer or repr

    @classmethod
    def base(cls, key, printer=None):
        return Vector({key: 1}, printer)

    def keys(self):
        return self.dictionary.keys()

    def values(self):
        return self.dictionary.values()

    def items(self):
        return self.dictionary.items()

    def __eq__(self, other):
        assert type(other) == Vector
        return len((self - other).dictionary) == 0

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def __add__(self, other):
        assert type(other) == Vector
        keys = self.keys() | other.keys()
        return Vector({key: self[key] + other[key] for key in keys}, self.printer)

    def __sub__(self, other):
        assert type(other) == Vector
        keys = self.keys() | other.keys()
        return Vector({key: self[key] - other[key] for key in keys}, self.printer)

    def __mul__(self, other):
        return Vector({key: self[key] * other for key in self.keys()}, self.printer)

    def __rmul__(self, other):
        return Vector({key: self[key] * other for key in self.keys()}, self.printer)

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
        base = ''.join(self._repr_coeff(value) + self.printer(key) for key, value in self.items())
        if base.startswith(' + '):
            return base[3:]
        elif base.startswith(' - '):
            return '-' + base[3:]
        else:
            return '0'
