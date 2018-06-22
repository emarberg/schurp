
class MarkedNumber:
    def __init__(self, i):
        assert type(i) == int
        self.number = i

    def weight(self):
        return abs(self.number)

    def is_marked(self):
        return self.number < 0

    def is_zero(self):
        return self.number == 0

    def __neg__(self):
        return MarkedNumber(-self.number)

    def increment(self):
        if self.number > 0:
            return MarkedNumber(self.number + 1)
        else:
            return MarkedNumber(self.number - 1)

    def __eq__(self, other):
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
            return "*"
        if self.number < 0:
            return str(-self.number) + "'"
        else:
            return str(self.number)

    def __hash__(self):
        return self.number
