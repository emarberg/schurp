from permutations import Permutation
from tableaux import Tableau
from partitions import Shape
from words import Word
from marked import MarkedNumber
import random


ck = Word.coxeter_knuth_move
pairing = Word._incr_pairing


def shifted_diagram(mu):
    return {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(mu[i])}


def shifted_horizontal_strips(mu, maxdiff=None):
    m = mu[0] if mu else 0
    maxdiff = m if maxdiff is None else maxdiff
    diagram = shifted_diagram(mu)
    cells = {(i + 1, j) for (i, j) in diagram if (i + 1, j) not in diagram and i + 1 <= j}
    cells |= {(1, m + i + 1) for i in range(maxdiff)}
    return sorted(cells, key=lambda x: (-x[0], x[1]))


def shifted_vertical_strips(mu):
    diagram = shifted_diagram(mu)
    cells = {(i, j + 1) for (i, j) in diagram if (i, j + 1) not in diagram and i <= j + 1}
    q = len(mu)
    cells.add((q + 1, q + 1))
    return sorted(cells, key=lambda x: (-x[1], x[0]))


def shifted_corners(mu):
    s = shifted_vertical_strips(mu)
    return [(i, j) for (i, j) in s if (i - 1, j) not in s]


def random_strict_partition(largest_part=20):
    ans = []
    next_part = largest_part
    while next_part > 0:
        ans.append(next_part)
        next_part = random.randint(next_part // 2, next_part - 1)
    return tuple(ans)


def random_vertical_strip(mu):
    ans = []
    v = shifted_vertical_strips(mu)
    i = 0
    while i < len(v):
        j = i
        while j < len(v) and v[i][-1] == v[j][-1]:
            j += 1
        for t in range(i, random.randint(i, j)):
            ans.append(v[t])
        i = j
    return ans


def random_horizontal_strip(mu):
    ans = []
    h = shifted_horizontal_strips(mu)
    i = 0
    while i < len(h):
        j = i
        while j < len(h) and h[i][0] == h[j][0]:
            j += 1
        for t in range(i, random.randint(i, j)):
            ans.append(h[t])
        i = j
    return ans


def random_outer_corner(mu):
    return random.choice(shifted_corners(mu))


def add_strict_partition_strip(mu, s):
    diagram = shifted_diagram(mu) | set(s)
    ans = []
    for (i, j) in diagram:
        while len(ans) < i:
            ans += [0]
        ans[i - 1] += 1
    return tuple(ans)


def test_iterated_tableau_operators(n, it=1000):
    for _ in range(it):
        _test_iterated_tableau_operators(n)


def _test_iterated_tableau_operators(n):
    mu = random_strict_partition(n)
    v = random_vertical_strip(mu)
    nu = add_strict_partition_strip(mu, v)
    h = random_horizontal_strip(nu)
    ou = add_strict_partition_strip(nu, h)
    c = random_outer_corner(ou)

    tab = {}

    diagram = shifted_diagram(mu)
    for (i, j) in diagram:
        tab[(i, j)] = 0

    e = 1
    for (i, j) in v:
        tab[(i, j)] = -e
        e += 1
    for (i, j) in h:
        tab[(i, j)] = e
        e += 1
    tab[c] = random.choice([-e, e])
    tab = Tableau(tab)

    if e < 3 or tab == tab.dual_equivalence_operator(e - 2):
        print('no descent')
        return

    ans = tab
    for i in range(e - 2, 0, -1):
        # print(ans)
        ans = ans.dual_equivalence_operator(i)

    w, p = tab.shifted_crystal_word()
    ii = [i for i in range(len(w)) if abs(w[i]) == e][0]
    jj = [i for i in range(ii + 1, len(w)) if abs(w[i]) != 0][0]

    sab = tab
    for (i, j) in sab:
        v = sab[(i, j)]
        if 0 != abs(v) < e:
            sab = sab.set(i, j, -2 if v.is_marked() else 2)
        elif 0 != abs(v):
            sab = sab.set(i, j, -3 if v.is_marked() else 3)

    (i, j) = p[jj]
    v = sab[(i, j)]
    sab = sab.set(i, j, -1 if v.is_marked() else 1)

    guess = sab
    i, j = list(guess.find(3, -3))[0]
    v = guess[(i, j)]
    if v.is_marked():
        if guess[(i - 1, j)] == MarkedNumber(2):
            guess = guess.set(i, j, 2).set(i - 1, j, -2)
        elif guess[(i, j - 1)] == MarkedNumber(2):
            k = 1
            while guess[(i, j - k - 1)] == MarkedNumber(2):
                k += 1
            if guess[(i, j - k - 1)] == MarkedNumber(-2):
                guess = guess.set(i, j, 2)
            else:
                guess = guess.set(i, j, 2).set(i, j - k, -2)
        else:
            guess = guess.set(i, j, -2)
    else:
        guess = guess.set(i, j, 2)
    i, j = list(guess.find(1, -1))[0]

    if guess[(i, j - 1)] in [MarkedNumber(2), MarkedNumber(-2)]:
        guess = guess.set(i, j, -2).set(i, j - 1, 1)
    elif guess[(i - 1, j)] in [MarkedNumber(2), MarkedNumber(-2)]:
        k = 1
        while guess[(i - k - 1, j)] in [MarkedNumber(2), MarkedNumber(-2)]:
            k += 1
        guess = guess.set(i, j, -2).set(i - k, j, 1)

    bns = ans
    for (i, j) in bns:
        v = bns[(i, j)]
        if 0 != abs(v) > 1:
            bns = bns.set(i, j, -2 if v.is_marked() else 2)
        elif 0 != abs(v):
            bns = bns.set(i, j, -1 if v.is_marked() else 1)

    if guess != bns:
        print(ans)
        print(tab.shifted_reading_word())
        print()
        print(ans.shifted_reading_word())
        print()

        # print(e, tab, ans)
        print(sab, guess, bns)
        assert False
    else:
        print('success')


def descents(w):
    return [i for i in range(len(w) - 1) if w[i] > w[i + 1]]


def rightop(w, i):
    for j in range(i, len(w)):
        w = ck(w, j)
    return w


def leftop(w, i):
    for j in range(i, -1, -1):
        w = ck(w, j)
    return w


def optest(w, i):
    for j in list(range(i, -1, -1)):
        w = leftop(w, j)
    return w


def generate(n):
    w = Permutation()
    a, b = [], []
    for i in range(1, n):
        if w(i) < w(i + 1) and random.random() > 0.5:
            a += [i]
            w *= Permutation.s_i(i)
    for i in range(1, n):
        if w(i) < w(i + 1) and random.random() > 0.5:
            b += [i]
            w *= Permutation.s_i(i)
    return tuple(a + b), tuple(a), tuple(b)


def test(n=10):
    w, a, b = generate(10)
    print('w =', w, '\n')
    print(pairing(a, b), '\n')

    fw = Word.incr_crystal_f((a, b), 1)
    fw = fw[0] + fw[1] if fw else None
    print('fw =', fw, '\n')

    ew = Word.incr_crystal_e((a, b), 1)
    ew = ew[0] + ew[1] if ew else None
    print('ew =', ew, '\n')

    for i in range(len(w)):
        v = optest(w, i)
        print(i, v, v == fw, v == ew)
