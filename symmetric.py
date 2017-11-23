from partitions import StrictPartition
from vectors import Vector


class SchurP:

    def __init__(self, strict_partition):
        self.mu = strict_partition

    def __repr__(self):
        return 'P%s' % str(self.mu.parts)

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
