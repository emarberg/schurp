from stable.tableaux import Tableau, Partition
from stable.utils import (
    GQ, GP,
    jp, decomposition_jp, 
)
from stable.polynomials import Polynomial, X


def marked_to_decomposition(tab):
    rows = tab.get_rows(unpack=False)
    newrows = []
    for row in rows:
        d = []
        for s in row:
            m = sorted(s, key=lambda x: 2 * x if x > 0 else -1 - 2 * x)[0]
            if m > 0:
                d = [[i for i in s if i > 0]] + d
                d[-1] += [-i for i in s if i < 0]
            else:
                d = d + [[-i for i in s if i < 0]]
                d[0] += [i for i in s if i > 0]
        newrows.append(d)
    return Tableau.from_rows(newrows, shifted=True)


def test_marked_to_decomposition(max_entry=5, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        for n in range(1, max_entry + 1):
            print('mu =', mu, 'n =', n)
            seen = set()
            tabs = Tableau.semistandard_shifted_marked(n, mu, diagonal_primes=False, setvalued=True)
            try:
                for t in tabs:
                    out = marked_to_decomposition(t)
                    assert out not in seen
                    assert out.is_decomposition_tableau()
                    seen.add(out)
            except:
                print('  failed')
                


def test_decomposition_jp(max_entry=5, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        for n in range(1, max_entry + 1):
            guess = decomposition_jp(n, mu)
            expected = jp(n, mu).polynomial()
            print(mu, n)
            assert guess == expected


def setvalued_decomposition_f(tab, index):
    return tab.f_operator_on_setvalued_decomposition_tableaux(index)


def setvalued_decomposition_e(tab, index):
    return tab.e_operator_on_setvalued_decomposition_tableaux(index)


def test_setvalued_decomposition_f(n=3, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        print(n, mu)
        seen = set()
        for tab in Tableau.all(n, mu, shifted=True, marked=False, setvalued=True):
            if tab.is_decomposition_tableau():
                for i in range(1, n):
                    # print(tab, i)
                    res = setvalued_decomposition_f(tab, i)
                    # print(res)
                    if not (res is None or res.is_decomposition_tableau()):
                        print(tab, i, res)
                        print(res.boxes)
                    assert res is None or res.is_decomposition_tableau()
                    if res is not None:
                        assert (i, res) not in seen
                        seen.add((i, res))
                        back = setvalued_decomposition_e(res, i) 
                        if back is None or back != tab:
                            print(res)
                            print('index =', i)
                            print(back)
                            print('expected:')
                            print(tab)
                        assert back is not None and back == tab


def GQ_decomposition_gf(n, mu):
    ans = 0
    summary = []
    for tab in Tableau.all(n, mu, shifted=True, marked=True, setvalued=True):
        if tab.is_primed_decomposition_tableau():
            #if 1 in tab.get(2, 2, unpack=False) and -1 in tab.get(2, 2, unpack=False):
            #    continue
            weight = tab.weight()
            summary += [(weight, tab)]
            weight = (sum(weight) - sum(mu),) + weight
            ans += Polynomial.from_tuple(weight)
        for weight, tab in sorted(summary):
            print(weight, tab)
    return ans


def test_GQ_decomposition_gf(max_entry=3, max_size=10):
    for mu in Partition.all(max_size, strict=True):
        for n in range(1, max_entry + 1):
            guess = GQ_decomposition_gf(n, mu)
            expected = GQ(n, mu).polynomial()
            print(mu, n)
            assert guess == expected


def GP_decomposition_gf(n, mu):
    ans = 0
    for tab in Tableau.all(n, mu, shifted=True, marked=False, setvalued=True):
        if tab.is_decomposition_tableau():
            weight = tab.weight()
            print(weight, tab)
            weight = (sum(weight) - sum(mu),) + weight
            ans += Polynomial.from_tuple(weight)
    return ans


def test_GP_decomposition_gf(max_entry=3, max_size=10):
    for mu in Partition.all(max_size, strict=True):
        for n in range(1, max_entry + 1):
            guess = GP_decomposition_gf(n, mu)
            expected = GP(n, mu).polynomial()
            print(mu, n)
            assert guess == expected