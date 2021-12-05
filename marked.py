
class MarkedNumber:
    def __init__(self, i):
        assert type(i) == int
        self.number = i

    def weight(self):
        return abs(self.number)

    def is_primed(self):
        return self.is_marked()

    def is_marked(self):
        return self.number < 0

    def is_zero(self):
        return self.number == 0

    @classmethod
    def swap_primes(cls, a, b):
        if a.is_primed() and b.is_primed():
            return (a, b)
        if not a.is_primed() and not b.is_primed():
            return (a, b)
        return (-a, -b)

    def __abs__(self):
        return abs(self.number)

    def __neg__(self):
        return MarkedNumber(-self.number)

    def __add__(self, other):
        if type(other) == int:
            other = MarkedNumber(other)
        return MarkedNumber(self.number + other.number)

    def __sub__(self, other):
        if type(other) == int:
            other = MarkedNumber(other)
        return MarkedNumber(self.number - other.number)

    def ceil(self):
        return MarkedNumber(abs(self))

    def floor(self):
        return MarkedNumber(-abs(self))

    def increment(self):
        if self.number > 0:
            return MarkedNumber(self.number + 1)
        else:
            return MarkedNumber(self.number - 1)

    def __eq__(self, other):
        if other is None:
            return False
        assert type(other) == MarkedNumber
        return self.number == other.number

    def __lt__(self, other):
        assert type(other) == MarkedNumber
        if abs(self.number) != abs(other.number):
            return abs(self.number) < abs(other.number)
        else:
            return self.number < 0 and other.number > 0

    def __le__(self, other):
        return self == other or self < other

    def __repr__(self):
        if self.number == 0:
            return "0"  # "*"
        if self.number < 0:
            return str(-self.number) + "'"
        else:
            return str(self.number)

    def __hash__(self):
        return self.number
