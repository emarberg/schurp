from partitions import StrictPartition
from vectors import Vector


class SchurP:

    def __init__(self, strict_partition):
        self.mu = strict_partition

    def __hash__(self):
        return hash(self.mu)

    def __eq__(self, other):
        assert type(other) == SchurP
        return self.mu == other.mu

    def __lt__(self, other):
        assert type(other) == SchurP
        return self.mu < other.mu

    def __repr__(self):
        return 'P(%s)' % ', '.join(str(i) for i in self.mu.parts)

    def __add__(self, other):
        if type(other) == Vector:
            return Vector.base(self) + other
        else:
            assert type(other) == SchurP
            return Vector.base(self) + Vector.base(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == Vector:
            return Vector.base(self) - other
        else:
            assert type(other) == SchurP
            return Vector.base(self) - Vector.base(other)

    def __rsub__(self, other):
        return other.__sub__(Vector.base(self))

    def __neg__(self):
        return -Vector.base(self)

    def __mul__(self, other):
        assert type(other) == SchurP and len(other.mu) == 1
        return Vector({SchurP(p): v for p, v in self.mu.pieri(other.mu(1)).items()})

    def __rmul__(self, i):
        self.__mul__(i)
