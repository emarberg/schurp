from stable.tableaux import Tableau, Partition
from stable.utils import (
    GQ, GP,
    jp, decomposition_jp, 
    schur_expansion,
    SymmetricPolynomial
)
from stable.polynomials import Polynomial, X
from stable.vst import combine_str
from crystals import AbstractGLCrystal, AbstractQCrystal
import time


def test_primed_decomposition_insert_relation():
    def span(w):
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if abs(a) <= abs(b):
                yield w[:i] + (a, -b) + w[i + 2:]
            else:
                yield w[:i] + (-a, b) + w[i + 2:]

    def generate(w):
        k = set()
        add = {w}
        while add:
            k |= add
            add = {x for v in add for x in span(v) if x not in k}
        return k

    from words import decomposition_insert
    nums = [-4, -3, -2, -1, 1, 2, 3, 4]
    seen = {}
    for a in nums:
        for b in nums:
            for c in nums:
                for d in nums:
                    t = decomposition_insert(d, c, b, a)[0]
                    if t not in seen:
                        seen[t] = []
                    seen[t].append((a, b, c, d))
    for k in seen.values():
        x = k[0]
        y = [y for y in k if tuple(abs(_) for _ in y) != tuple(abs(_) for _ in x)]
        y = y[0] if y else x
        if generate(x) | generate(y) != set(k):
            print(k)
            print(generate(k[0]))
        assert generate(x) | generate(y) == set(k)
    return seen


def test_sv_decomposition_generator(n=3, k=5):
    for mu in Partition.all(k, strict=True):
        t0 = time.time()
        expected = set(Tableau.setvalued_decomposition_tableaux_slow(n, mu))
        t1 = time.time()
        received = set(Tableau.setvalued_decomposition_tableaux(n, mu))
        t2 = time.time()
        print('n =', n, 'mu =', mu, (t1 - t0) / (t2 - t1))
        if expected != received:
            print()
            print(expected - received)
            print(received - expected)
            print()
        assert expected == received


def test_decomposition_semicrystal(n=3, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        if len(mu) > n:
            continue
        c = AbstractQCrystal.decomposition_semicrystal_from_strict_partition(mu, n)
        print(n, mu, len(c))
        nhw = c.naive_highest_weights()
        ahw = c.get_highest_weights()
        print('  ', c.is_connected(), len(ahw), len(nhw))
        assert len(ahw) == 1 or not mu
        #if len(ahw) > 1:
        #    c.draw(extended=True)
        #    input('')


def test_sv_decomposition_crystal(n=3, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        c = AbstractQCrystal.setvalued_decomposition_tableaux_from_strict_partition(mu, n)
        nhw = c.naive_highest_weights()
        ahw = c.get_highest_weights()
        print(n, mu, len(c), c.is_connected(), len(ahw), len(nhw))
        if len(ahw) > 1:
            c.draw(extended=True)
            input('')


def test_decomposition_stembridge(n=3, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        c = AbstractQCrystal.decomposition_tableaux_from_strict_partition(mu, n)
        print(n, mu, c.is_stembridge())


def test_setvalued_decomposition_stembridge(n=3, max_size=5):
    for mu in Partition.all(max_size, strict=True):
        c = AbstractQCrystal.setvalued_decomposition_tableaux_from_strict_partition(mu, n)
        b = c.is_stembridge()

        ch = sum([Polynomial.from_tuple((0,) + a) for a in c.weights.values()])
        ch = SymmetricPolynomial.from_polynomial(ch)
        u = schur_expansion(ch).dictionary
        v = {Partition.trim(mu): x for (mu, x) in c.get_highest_weight_multiplicities().items()}

        print(n, mu, b, u == v)
        if not b:
            b = c.is_stembridge(return_counterexample=True)
            print()
            print(len(b), len(c), len(b) == len(c))
            print()
            b.draw()
            input('')


def test_setvalued_semistandard_stembridge(n=3, max_size=5):
    for mu in Partition.all(max_size):
        c = AbstractGLCrystal.setvalued_semistandard_tableaux_from_partition(mu, n)
        b = c.is_stembridge()

        ch = sum([Polynomial.from_tuple((0,) + a) for a in c.weights.values()])
        ch = SymmetricPolynomial.from_polynomial(ch)
        u = schur_expansion(ch).dictionary
        v = {Partition.trim(mu): x for (mu, x) in c.get_highest_weight_multiplicities().items()}

        print(n, mu, b, u == v)


def test_primed_insert(n=5, mu=(3,2)):
    from tableaux import Tableau as T
    from words import decomposition_insert
    for tab in T.get_semistandard_decomposition(n, mu, primed=False):
        word = list(tab.row_reading_word())
        print('unprimed word =', word)
        print()
        print(tab)
        print()
        for v in range(2**len(word)):
            newword = word[:]
            for i in range(len(word)):
                newword[i] *= -1 if v % 2 else 1
                v = v // 2
            tabs = [decomposition_insert(*newword[:i + 1])[0] for i in range(len(word))]
            label = [str(newword)]
            label += (len(str(tab).split('\n')) - 1) * [' ' * len(label[-1])]
            print(combine_str('\n'.join(label), *tabs))
            print()
            tabs = [decomposition_insert(*newword[:i + 1])[1] for i in range(len(word))]
            label = [str(newword)]
            label += (len(str(tab).split('\n')) - 1) * [' ' * len(label[-1])]
            print(combine_str('\n'.join(label), *tabs))
            print()
            print()
            print()
            print()
        input('')


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


def test_setvalued_decomposition(n=3, max_size=5):
    for mu in [(3,2)]:#Partition.all(max_size, strict=True):
        print(n, mu)
        f_seen = set()
        e_seen = set()
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
                        assert (i, res) not in f_seen
                        f_seen.add((i, res))
                        back = setvalued_decomposition_e(res, i) 
                        if back is None or back != tab:
                            print(res)
                            print('index =', i)
                            print(back)
                            print('expected:')
                            print(tab)
                        assert back is not None and back == tab
                    #
                    #
                    #
                    res = setvalued_decomposition_e(tab, i)
                    if not (res is None or res.is_decomposition_tableau()):
                        print(tab, i, res)
                        print(res.boxes)
                    assert res is None or res.is_decomposition_tableau()
                    if res is not None:
                        assert (i, res) not in e_seen
                        e_seen.add((i, res))
                        back = setvalued_decomposition_f(res, i) 
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
