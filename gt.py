from cached import cached_value
from schubert import X
from stable.partitions import Partition
import itertools

GT_CACHE = {}


class GTPattern:

    def __init__(self, array=(), compute_weight=None):
        assert all(len(row) == len(array) - i for i, row in enumerate(array))
        
        self.rows = tuple(tuple(row) for row in array)
        self._string = None
        self._length = len(self.rows)
        
        if compute_weight is None:
            w = [0]
            for row in reversed(self.rows):
                w.append(sum(row) - w[-1])
            self._weight = tuple(w[1:])
        else:
            self._weight = compute_weight(self.rows)

    def __str__(self):
        if self._string is None:
            space = '.'
            pad = max([0] + [len(str(x)) for row in self.rows for x in row])
            gap = (pad + 2) * space
            ans = []
            for i in range(len(self)):
                print('gap =', gap)
                s = i * gap + gap.join([space + (pad - len(str(x))) * space + str(x) + space for x in self.rows[i]])
                ans.append(s[1:])
            self._string = '\n' + '\n'.join(ans) + '\n'
        return self._string

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self.rows)

    def __eq__(self, other):
        return self.weight == other.weight and self.rows == other.rows

    def __len__(self):
        return self._length

    @property
    def length(self):
        return len(self)

    @property
    def weight(self):
        return self._weight

    @property
    def monomial(self):
        ans = X(0) ** 0
        for i, w in enumerate(self.weight):
            ans *= X(i + 1) ** w
        return ans

    def __getitem__(self, *args):
        (i, j) = args if len(args) == 2 else (args[0][0], args[0][1])
        return self.get(i, j)

    def get(self, i, j):
        return self.rows[i][j]

    def is_strict(self):
        return all(self[i, j - 1] < self[i, j] for i in range(len(self)) for j in range(1, len(self.rows[i])))

    @classmethod
    def from_partition(cls, mu, size=None):
        if size is not None:
            assert size >= len(mu)
            mu = list(mu) + (size - len(mu)) * [0]
        return cls.all(mu)

    @classmethod
    def all(cls, *args):
        if len(args) == 1 and type(args[0]) in [list, tuple]:
            args = args[0]
        
        assert all(type(args[i]) == int for i in range(len(args)))
        assert all(args[i] >= 0 for i in range(len(args)))
        
        args = tuple(sorted(args))
        return cls._all(args)

    @cached_value(GT_CACHE)
    def _all(cls, mu):
        if len(mu) == 0:
            return [GTPattern()]
        if len(mu) == 1:
            return [GTPattern([mu])]
        ans = []
        intervals = [range(mu[i - 1], mu[i] + 1) for i in range(1, len(mu))]
        for nu in itertools.product(*intervals):
            ans.extend([GTPattern((mu,) + x.rows) for x in cls._all(nu)])
        return ans

