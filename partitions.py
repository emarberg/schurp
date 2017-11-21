

class StrictPartition:
    def __init__(self, *args):
        self.parts = tuple(sorted(args, reverse=True))
        assert len(set(self.parts)) == len(self.parts)
        assert self.parts == tuple(args)

    def __repr__(self):
        return '\n'.join(i * '  ' + self.parts[i] * '* ' for i in range(len(self.parts)))

    @property
    def shape(self):
        return {(i + 1, i + j) for i in range(len(self.parts)) for j in range(1, self.parts[i] + 1)}

    def row(self, i):
        if not (1 <= i <= self.num_rows):
            return ()
        return tuple((i, i + j) for j in range(self.parts[i - 1]))

    @property
    def num_rows(self):
        return len(self.parts)

    def column(self, j):
        if not (1 <= j <= self.num_columns):
            return ()
        ans = []
        for i in range(1, j + 1):
            if j <= self.parts[i - 1] + i - 1:
                ans.append((i, j))
        return tuple(ans)

    @property
    def num_columns(self):
        if self.parts:
            return max(self.parts)
        else:
            return 0
