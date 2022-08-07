from permutations import Permutation
from qp_utils import beissinger_rsk, rsk


def test_all_row(n=5):
    for w in Permutation.involutions(n):
        w = [w(i + 1) for i in range(n)]
        btab = beissinger_rsk(w, sgn=False)
        p, q = rsk(w)
        assert p == q
        rsktab = p
        # print(w)
        # print(btab)
        # print(rsktab)
        assert btab == rsktab 


def test_all_col(n=5):
    seen = {}
    for w in Permutation.involutions(n):
        w = [w(i + 1) for i in range(n)]
        btab = beissinger_rsk(w, sgn=True)
        rsktab = rsk(w)[0]
        print(w)
        print(btab)
        print(rsktab)
        print()
        assert btab not in seen
        seen[btab] = tuple(w)

    for tab in seen:
        print(seen[tab])
        print(tab)
        for i in range(2, n):
            print('D_%s' % i)
            print()
            alt = tab.dual_equivalence_operator(i - 1)
            print(seen[alt])
            print(alt)
            print()
        print()