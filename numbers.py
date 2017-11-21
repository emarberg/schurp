
class MarkedNumber:
    def __init__(self, i):
        assert type(i) == int and i != 0
        self.number = i

    def weight(self):
        return abs(self.number)

    def is_marked(self):
        return self.number < 0

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
        if self.number < 0:
            return str(-self.number) + "'"
        else:
            return str(self.number)

    def __hash__(self):
        return self.number

    def __len__(self):
        return len(str(self))
