from pipedreams import BumplessPipedream
from permutations import Permutation

# TODO

def test_get_tile():
    pass


def test_is_blank():
    pass


def test_get_blank_tiles():
    pass


def test_get_minimal_blank_tile():
    pass


def test_get_pipe(n=5):
    for w in Permutation.all(n):
        for b in BumplessPipedream.from_permutation(w, n):
            for wi in range(1, n):
                j = n
                i = w.inverse()(wi)
                print(w)
                print(b)
                x = b.get_pipe(wi, j, 'H')
                print(wi, j, '-->', i, '=?=', x)
                print()
                assert x == i


def test_delta():
    pass

