

class SignedPermutation:

    def __init__(self, *oneline):
        self.oneline = oneline
        self.rank = len(oneline)
        assert tuple(range(1, self.rank + 1)) == tuple(abs(i) for i in self.oneline)

    def __repr__(self):
        s = []
        for i in self.oneline:
            s += [str(abs(i))]
            if i < 0:
                s += ['\u0305']
        return ''.join(s)
