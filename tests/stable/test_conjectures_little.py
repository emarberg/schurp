from stable.permutations import Permutation
from stable.insertion import InsertionAlgorithm
from stable.tableaux import Tableau
from itertools import chain, combinations
from stable.partitions import Partition
from stable.words import Word
import pytest


def powerset(iterable, maxsize=None):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    maxsize = len(s) if maxsize is None else maxsize
    return chain.from_iterable(combinations(s, r) for r in range(maxsize + 1))


@pytest.mark.slow
def test_strong_little(n=5, k=5, maxcount=20):
    cases = 0
    for b, pi in enumerate(Permutation.all(n)):
        print(b, len(pi), pi, cases)
        for c, w in enumerate(pi.get_reduced_words()):
            if c > 0 and c % 10 == 0:
                print('  ', c)
            if c > maxcount:
                break
            p, q = InsertionAlgorithm.hecke(w)
            for i in powerset(range(len(w)), k):
                v = list(w)
                for x in i:
                    v[x] += 1
                if len(Permutation.from_word(v)) == len(v):
                    cases += 1
                    pp, qq = InsertionAlgorithm.hecke(v)
                    assert q == qq
                    assert all(pp.get(i, j) - p.get(i, j) in [0, 1] for i, j in p.boxes)
    assert cases > 0


@pytest.mark.slow
def test_strong_inv_little(n=6, k=5, maxcount=20):
    cases = 0
    for b, pi in enumerate(Permutation.involutions(n)):
        print(b, pi.involution_length, pi, cases)
        for c, w in enumerate(pi.get_involution_words()):
            if c > 0 and c % 10 == 0:
                print('  ', c)
            if c > maxcount:
                break
            p, q = InsertionAlgorithm.orthogonal_hecke(w)
            for i in powerset(range(len(w)), k):
                v = list(w)
                for x in i:
                    v[x] += 1
                if Permutation.from_involution_word(v, strict=True) is not None:
                    cases += 1
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(v)
                    assert q == qq
                    assert all(pp.get(i, j) - p.get(i, j) in [0, 1] for i, j in p.boxes)
    assert cases > 0


@pytest.mark.slow
def test_strong_fpf_little(n=6, k=5, maxcount=20):
    assert n % 2 == 0
    cases = 0
    minfpf = Permutation.minfpf(n)

    def invol(w, i):
        v = [w[j] + 1 if j in i else w[j] for j in range(len(w))]
        if len(Permutation.from_word(v)) != len(v):
            return None, None
        z = minfpf
        v = []
        for j in range(len(w)):
            a = w[j] + (1 if j in i else 0)
            if z(a) == a + 1:
                a += 1
            if z(a) > z(a + 1):
                return None, None
            s = Permutation.s_i(a)
            z = s * z * s
            v += [a]
        return z, v

    for b, pi in enumerate(Permutation.fpf_involutions(n)):
        print(b, pi.fpf_involution_length, pi, cases)
        for c, w in enumerate(pi.get_fpf_involution_words()):
            if c > 0 and c % 10 == 0:
                print('  ', c)
            if c > maxcount:
                break
            p, q = InsertionAlgorithm.symplectic_hecke(w)
            for i in powerset(range(len(w)), k):
                z, v = invol(w, i)

                if z is None:
                    continue
                # if not all((v[j] > v[j + 1]) == (w[j] > w[j + 1]) for j in range(len(w) - 1)):
                #     continue

                pp, qq = InsertionAlgorithm.symplectic_hecke(v)
                if q != qq or not all(pp.get(i, j) - p.get(i, j) in [0, 1, 2] for i, j in p.boxes):
                    print()
                    print('*', w, '->', v)
                    print(p)
                    print(q)
                    print()
                    print(pp)
                    print(qq)
                    input('?')
                else:
                    cases += 1
                # assert q == qq
                # assert all(pp.get(i, j) - p.get(i, j) in [0, 1, 2] for i, j in p.boxes)
    assert cases > 0


def test_push():
    word = (1, 3, 4, 3, 2, 3, 1)
    i = 6

    new, j = (1, 3, 4, 3, 2, 3, 2), 3
    assert Permutation.little_push(word, i) == (new, j)
    word, i = new, j

    new, j = (1, 3, 4, 4, 2, 3, 2), 2
    assert Permutation.little_push(word, i) == (new, j)
    word, i = new, j

    new, j = (1, 3, 5, 4, 2, 3, 2), None
    assert Permutation.little_push(word, i) == (new, j)


def test_bump():
    word = (1, 3, 4, 3, 2, 3, 1)
    end = (1, 3, 5, 4, 2, 3, 2)
    j, k = 1, 2
    assert Permutation.little_bump(word, j, k) == end


def get_ck_type_two(word):
    for i in range(len(word) - 2):
        a, c, b = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (c, a, b) + word[i + 3:]
        c, a, b = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (a, c, b) + word[i + 3:]


def get_ck_type_one(word):
    for i in range(len(word) - 2):
        b, a, c = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (b, c, a) + word[i + 3:]
        b, c, a = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (b, a, c) + word[i + 3:]


def get_ck_type_three(word):
    for i in range(len(word) - 2):
        b, a, bb = word[i:i + 3]
        if a < b == bb:
            yield i, word[:i] + (a, b, a) + word[i + 3:]
        a, b, aa = word[i:i + 3]
        if a == aa < b:
            yield i, word[:i] + (b, a, b) + word[i + 3:]


def get_ck_type_inv(word):
    if len(word) >= 2 and word[0] != word[1]:
        yield -1, (word[1], word[0]) + word[2:]


def get_ck_type_fpf(word):
    if len(word) >= 2 and word[0] != word[1]:
        if word[0] % 2 == word[1] % 2:
            yield -1, (word[1], word[0]) + word[2:]
        elif word[1] == word[0] - 1:
            yield -1, (word[0], word[1] + 2) + word[2:]
        elif word[1] == word[0] + 1:
            yield -1, (word[0], word[1] - 2) + word[2:]


def get_inv_ck_moves(word):
    for i, v in get_ck_moves(word):
        yield i, v
    for i, v in get_ck_type_inv(word):
        yield i, v


def get_fpf_ck_moves(word):
    for i, v in get_ck_moves(word):
        yield i, v
    for i, v in get_ck_type_fpf(word):
        yield i, v


def get_ck_moves(word):
    for i, v in get_ck_type_one(word):
        yield i, v
    for i, v in get_ck_type_two(word):
        yield i, v
    for i, v in get_ck_type_three(word):
        yield i, v


def get_bruhat_covers(w):
    n = w.rank
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            t = Permutation.transposition(j, k)
            if len(w * t) == len(w) - 1:
                yield j, k


def get_inv_bruhat_covers(w):
    n = w.rank
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            t = Permutation.transposition(j, k)
            subatoms = [x * t for x in w.get_atoms() if len(x * t) == len(x) - 1]
            subatoms = [x for x in subatoms if (x.inverse() % x).involution_length == len(x)]
            if subatoms:
                yield j, k


def adjacent_unmatched_primes(q, k):
    if k in q.values() and -k - 1 in q.values():
        (x1, y1) = q.find(k)[0]
        (x2, y2) = q.find(-k - 1)[0]
        return x1 == x2 or y1 == y2
    elif -k in q.values() and k + 1 in q.values():
        (x1, y1) = q.find(-k)[0]
        (x2, y2) = q.find(k + 1)[0]
        return x1 == x2 or y1 == y2
    elif k in q.values() and k + 1 in q.values():
        (x1, y1) = q.find(k)[0]
        (x2, y2) = q.find(k + 1)[0]
        return x1 == x2 and x1 == y1
    return False


def on_diagonal(q, k):
    if k in q.values():
        i, j = q.find(k)[0]
        return i == j
    return False


def test_ck_moves():
    for w in Permutation.all(4):
        covers = list(get_bruhat_covers(w))

        for word in w.get_reduced_words():
            alphas = {}
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_type_one(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print('  type I:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_three(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print('type III:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print(' type II:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                assert q == s * qq or q == t * qq
                alphas[i] = s if q == s * qq else t

            for j, k in covers:
                bumped = Permutation.little_bump(word, j, k)
                p, q = InsertionAlgorithm.hecke(bumped)
                label = 'bump(%s, %s, %s) =' % (str(word), j, k)

                moves = []
                for i, other in get_ck_type_one(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print('  type I:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_three(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print('type III:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_two(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print(' type II:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                assert sorted(moves) == sorted(alphas)


def test_unique_transposition():
    for w in Permutation.all(3):
        for word in w.get_reduced_words():
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_moves(word):
                pp, qq = InsertionAlgorithm.hecke(other)

                print(word, '<->', other, 'at i =', i)
                print(q)
                print(qq)

                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                pi = [s, t, s * t * s]

                a = [x for x in pi if qq == x * q]
                b = [x for x in pi if qq.descent_set() == (x * q).descent_set() and (x * q).is_semistandard()]
                print('a =', a)
                print('b =', b)
                print()
                assert a == b
                assert len(a) == 1


def test_unique_transposition_invol():
    for w in Permutation.involutions(5):
        for word in w.get_involution_words():
            p, q = InsertionAlgorithm.orthogonal_hecke(word)
            for i, other in get_inv_ck_moves(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)

                if i == -1 and other == (word[1], word[0]) + word[2:]:
                    r = Permutation(1, -2)
                    pi = [r]
                else:
                    s = Permutation.s_i(i + 1)
                    t = Permutation.s_i(i + 2)

                    f = Permutation.reflection_s(i + 1)
                    g = Permutation.reflection_s(i + 2)
                    h = Permutation.reflection_s(i + 3)

                    pi = []
                    if adjacent_unmatched_primes(q, i + 1):
                        pi += [g] if on_diagonal(q, i + 1) else [f * g]
                    else:
                        pi.append(s)
                    if adjacent_unmatched_primes(q, i + 2):
                        pi += [h] if on_diagonal(q, i + 2) else [g * h]
                    else:
                        pi.append(t)

                a = set(x for x in pi if qq == x * q)
                b = set(x for x in pi if qq.descent_set() == (x * q).descent_set() and (x * q).is_semistandard())

                if len(b) != 1:
                    print(pi)
                    print(word, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('a =', a)
                    print('b =', b)
                    print()

                assert a == b
                assert len(b) == 1


def test_inv_ck_moves():
    for w in Permutation.involutions(4):
        covers = list(get_inv_bruhat_covers(w))

        for word in w.get_involution_words():
            alphas = {}
            p, q = InsertionAlgorithm.orthogonal_hecke(word)
            for i, other in get_ck_type_one(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print('  type I:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_three(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print('type III:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print(' type II:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                if adjacent_unmatched_primes(q, i + 1):
                    if on_diagonal(q, i + 1):
                        t = Permutation.reflection_s(i + 2)
                    else:
                        t = Permutation.reflection_s(i + 1) * Permutation.reflection_s(i + 2)
                else:
                    t = Permutation.s_i(i + 1)
                assert q == s * qq or q == t * qq
                alphas[i] = s if q == s * qq else t

            for i, other in get_ck_type_inv(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print(' type IV:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                r = Permutation(1, -2)
                assert q == r * qq

            for j, k in covers:
                bumped = Permutation.involution_little_bump(word, j, k)
                p, q = InsertionAlgorithm.orthogonal_hecke(bumped)
                label = 'ibump(%s, %s, %s) =' % (str(word), j, k)

                moves = []
                for i, other in get_ck_type_one(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print('  type I:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_three(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print('type III:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_two(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print(' type II:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_inv(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print(' type IV:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    s = Permutation(1, -2)
                    print('s =', s)
                    print()
                    assert q == s * qq

                assert sorted(moves) == sorted(alphas)


def test_ck_moves_type_two():
    for w in Permutation.all(3):
        for word in w.get_reduced_words():
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                a = Tableau({(x, y): v for x, y, v in q if i + 3 in v or i + 1 in v or i + 2 in v})
                aa = Tableau({(x, y): v for x, y, v in qq if i + 3 in v or i + 1 in v or i + 2 in v})
                assert q == s * qq or q == t * qq
                if q == t * qq:
                    print('type IIb:', word, '<->', other, 'at i =', i + 1)
                    print(p)
                    print(a)
                    print(aa)
                else:
                    print('type IIa:', word, '<->', other, 'at i =', i + 1)
                    print(p)
                    print(a)
                    print(aa)


def test_tab():
    for w in Permutation.involutions(3):
        if w.is_inv_grassmannian():
            for word in w.get_involution_words():
                p, q = InsertionAlgorithm.orthogonal_hecke(word)
                print(w)
                print(word)
                print(q)
                print()


def test_q_tab():
    n = 4
    for pi in Permutation.all(n):
        for w in pi.get_reduced_words():
            _, q = InsertionAlgorithm.hecke(w)
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.little_bump(w, i, j)
                    if w != bumped:
                        _, qq = InsertionAlgorithm.hecke(bumped)
                        print(pi)
                        print(i, j, ':', w, '->', bumped)
                        print(q)
                        print(qq)
                        assert q == qq


@pytest.mark.slow
def test_inv_q_tab():
    n = 5
    countmax = 5
    for pi in Permutation.involutions(n):
        if len(pi) == 0:
            continue
        for count, w in enumerate(pi.get_involution_words()):
            if count >= countmax:
                break
            _, q = InsertionAlgorithm.orthogonal_hecke(w)
            cases = 0
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.involution_little_bump(w, i, j)
                    if w != bumped:
                        cases += 1
                        _, qq = InsertionAlgorithm.orthogonal_hecke(bumped)
                        print(pi)
                        print(i, j, ':', w, '->', bumped)
                        print(q)
                        print(qq)
                        assert q == qq
            assert cases != 0


@pytest.mark.slow
def test_inv_bumping_multiplicity():
    n = 5
    countmax = 100
    for pi in Permutation.involutions(n):
        for count, w in enumerate(pi.get_involution_words()):
            if count >= countmax:
                break
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.involution_little_bump(w, i, j)
                    if w != bumped:
                        assert {bumped[e] - w[e] for e in range(len(w))}.issubset({0, 1})


@pytest.mark.slow
def test_fpf_bumping_multiplicity():
    n = 4
    countmax = 100
    for pi in Permutation.fpf_involutions(n):
        for count, w in enumerate(pi.get_fpf_involution_words()):
            if count >= countmax:
                break
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.fpf_involution_little_bump(w, i, j)
                    if w != bumped:
                        assert {bumped[e] - w[e] for e in range(len(w))}.issubset({0, 1, 2})


@pytest.mark.slow
def test_fpf_q_tab():
    n = 4
    countmax = None
    for pi in Permutation.fpf_involutions(n):
        if len(pi) == 0:
            continue
        for count, w in enumerate(pi.get_fpf_involution_words()):
            if len(w) == 0:
                continue
            if countmax and count >= countmax:
                break
            _, q = InsertionAlgorithm.symplectic_hecke(w)
            cases = 0
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.fpf_involution_little_bump(w, i, j)
                    if w != bumped:
                        cases += 1
                        _, qq = InsertionAlgorithm.symplectic_hecke(bumped)
                        if q != qq:
                            print(pi)
                            print(i, j, ':', w, '->', bumped)
                            print(q)
                            print(qq)
                        assert q == qq
            assert cases != 0


@pytest.mark.slow
def test_dual_equiv():
    for w in Permutation.involutions(7):
        print()
        print('w =', w)
        words = w.get_involution_words()
        tab = {
            v: InsertionAlgorithm.orthogonal_hecke(v)[1]
            for v in words
        }
        for v in words:
            for i, u in get_inv_ck_moves(v):

                a = tab[v].shifted_reading_word()
                b = tab[u].shifted_reading_word()

                loc = {x - 1: i for i, x in enumerate(a)}

                print(tab[v])
                print(tab[u])
                print('i =', i + 1)
                print()
                print(' ' + i * '   ' + '*')
                print(v, '          ', a)
                print(u, '          ', b)
                print()

                if i == -1:
                    print('Case IV')
                    t = Permutation.reflection_s(2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                    continue

                x, y, z = loc[i], loc[i + 1], loc[i + 2]
                if x < z < y or y < z < x:
                    print('Case I')
                    t = Permutation.s_i(i + 1) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 1) * Permutation.reflection_s(i + 2) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                elif y < x < z or z < x < y:
                    print('Case II')
                    t = Permutation.s_i(i + 2) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                elif x < y < z or z < y < x:
                    print('Case III')
                    assert False
                else:
                    raise Exception
    print('done')


@pytest.mark.slow
def test_fpf_dual_equiv():
    for w in Permutation.fpf_involutions(6):
        print()
        print('w =', w)
        words = w.get_fpf_involution_words()
        tab = {
            v: InsertionAlgorithm.symplectic_hecke(v)[1]
            for v in words
        }
        for v in words:
            for i, u in get_fpf_ck_moves(v):

                a = tab[v].shifted_reading_word()
                b = tab[u].shifted_reading_word()

                loc = {x - 1: i for i, x in enumerate(a)}

                print(tab[v])
                print(tab[u])
                print('i =', i + 1)
                print()
                print(' ' + i * '   ' + '*')
                print(v, '          ', a)
                print(u, '          ', b)
                print()

                if i == -1:
                    print('Case IV')
                    t = Permutation.reflection_s(2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                    continue

                x, y, z = loc[i], loc[i + 1], loc[i + 2]
                if x < z < y or y < z < x:
                    print('Case I')
                    t = Permutation.s_i(i + 1) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 1) * Permutation.reflection_s(i + 2) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                elif y < x < z or z < x < y:
                    print('Case II')
                    t = Permutation.s_i(i + 2) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3) * tab[v]
                    if not t.is_semistandard():
                        t = Permutation.reflection_s(i + 2) * tab[v]
                    if tab[u] != t:
                        print(t)
                    assert tab[u] == t
                elif x < y < z or z < y < x:
                    print('Case III')
                    assert False
                else:
                    raise Exception
    print('done')


def ck(i, w):
    if i == -1:
        assert len(w) >= 2
        a, b = w[:2]
        if a % 2 == b % 2 == 0:
            return (b, a) + w[2:]
        if a % 2 == 0 and b in [a - 1, a + 1]:
            pre = (a, a + 1) if b < a else (a, a - 1)
            return pre + w[2:]
        raise Exception
    else:
        assert i + 2 < len(w)
        pre, mid, pos = w[:i], w[i:i + 3], w[i + 3:]

        while True:
            c, a, b = mid
            if a < b < c:
                mid = (a, c, b)
                break
            a, c, b = mid
            if a < b < c:
                mid = (c, a, b)
                break

            b, c, a = mid
            if a < b < c:
                mid = (b, a, c)
                break
            b, a, c = mid
            if a < b < c:
                mid = (b, c, a)
                break

            a, b, c = mid
            if a == c:
                mid = (b, a, b)
                break

            raise Exception

        return pre + mid + pos


def descent(i, w):
    assert 0 <= i < len(w) - 1
    return w[i] > w[i + 1]


def ck_noop(i, w, vee=True):
    try:
        if i >= 0 and not vee:
            assert not (w[i] > w[i + 1] < w[i + 2])
        return ck(i, w)
    except:
        return w


def get_rows(word):
    rows = []
    for x in word:
        if len(rows) == 0 or rows[-1][-1] > x:
            rows += [(x,)]
        else:
            rows[-1] += (x,)
    return rows


def get_columns(word):
    cols = []
    for x in word:
        if len(cols) == 0 or cols[-1][-1] < x:
            cols += [(x,)]
        else:
            cols[-1] += (x,)
    return cols


def to_row_reading(word):
    columns = get_columns(word)
    t = Tableau()
    for j, col in enumerate(columns):
        for i, v in enumerate(reversed(col)):
            t = t.add(i + 1, j + 1, v)
    row = t.row_reading_word()
    ans, col = to_column_reading(row)
    assert col == word
    return list(reversed(ans)), row


def to_column_reading(word):
    rows = get_rows(word)
    if len(rows) <= 1:
        return [], word

    starts = []
    for r in rows[:-1]:
        starts += [len(r) + (starts[-1] if starts else -1)]
    m = len(rows[-1])

    starts = [starts[i] for i in range(len(rows) - 1) if len(rows[i]) + len(rows) - i - 1 == m]
    m = len(starts) + 1

    ans = []
    for i, a in enumerate(starts):
        for b in range(i + 1):
            c = (starts[i + 1] - 1) if i + 1 < len(starts) else len(word) - 2
            for x in range(a - b, c - b):
                ans += [x]
                word = ck(x, word)
    subword = word[:-m]
    bns, newword = to_column_reading(subword)

    ans += bns
    word = newword + word[-m:]
    return ans, word


def ck_compute_functional(w, a):

    def des(i, w):
        return w[i - 1] > w[i]

    def compose(*args):
        def f(x):
            for fun in args:
                x = fun(x)
            return x
        return f

    def convert_to_col(mu):
        if len(mu) <= 1:
            return []
        starts = []
        for r in list(reversed(mu[1:])):
            starts += [r + (starts[-1] if starts else -1)]
        starts = [starts[i] for i in range(len(mu) - 1) if mu[len(mu) - 1 - i] + len(mu) - i - 1 == mu[0]]

        ans = []
        for i, a in enumerate(starts):
            for b in range(i + 1):
                c = (starts[i + 1] - 1) if i + 1 < len(starts) else (sum(mu) - 2)
                ans += list(range(a - b + 1, c - b + 1))

        nu = tuple((mu[i] - 1) if mu[i] + i == mu[0] else mu[i] for i in range(len(mu)))
        while nu and nu[-1] == 0:
            nu = nu[:-1]

        return ans + convert_to_col(nu)

    def convert_to_row(mu):
        return list(reversed(convert_to_col(mu)))

    if len(w) == 0:
        return compose(), Tableau()
    if len(w) == 1:
        return compose(), Tableau().add(1, 1, 1)

    op, Q = ck_compute_functional(w[:-1], a)

    mu = Q.shape()
    n = len(w) - 1
    r = len(mu)

    d = {i: n - sum(mu[:i - 1]) for i in range(1, r + 2)}

    def rho(i):
        f = op
        for j in range(n - 1, d[i], -1):
            f = compose(f, a(j))
        return f

    for i in range(1, r + 1):
        if not des(d[i], rho(i)(w)):
            return rho(i), Q.add(i, i - 1 + mu[i - 1] + 1, n + 1)

    f = rho(r)
    for j in range(mu[r - 1] - 1):
        f = compose(f, a(j))

    if not des(d[r], f(w)):
        return rho(r + 1), Q.add(r + 1, r + 1, n + 1)

    def psi(i):
        return compose(*[a(j) for j in range(d[i] + i - 2, 0, -1)])

    def phi(i):
        delta = (i + 1) * i // 2
        factors = i * [compose(*[a(j) for j in range(delta + r - 1, delta, -1)])]
        return compose(*factors)

    Psi = compose(*[psi(i) for i in range(1, r + 1)])

    factors = r * [compose(*[a(j) for j in range(r)])]
    Theta = compose(*factors)

    Phi = compose(*[phi(i) for i in range(1, r)])

    nu = tuple(mu[i] - 1 for i in range(len(mu)))
    nu = nu if nu[-1] != 0 else nu[:-1]
    indices = convert_to_col(nu)
    Gamma = compose(*[a(r + 1 + j) for j in indices])

    q = mu[0]
    h = Q.transpose().shape() + (0,)
    e = {i: sum(h[:i - 1]) for i in range(1, q + 2)}

    def gamma(i):
        f = compose(op, Psi, Theta, Gamma, Phi)
        delta = r * (r + 1) // 2
        for j in range(delta + 1, e[i]):
            f = compose(f, a(j))
        return f

    i = None
    for index in range(r + 1, q + 1):
        if des(e[index] + 1, gamma(index)(w)):
            i = index
            break
    if i is None:
        i = q + 1
    f = gamma(i)

    newQ = Q.add(h[i - 1] + 1, i, -n - 1)
    nu = newQ.shape()
    for j in convert_to_row(nu):
        f = compose(f, a(j))
    return f, newQ


def sp_ck_compute_functional(w):

    def a(i):
        return lambda x: ck_noop(i - 1, x)

    return ck_compute_functional(w, a)


def o_ck_compute_functional(w):

    def a(i):
        def f(x):
            if i == 0:
                if len(x) >= 2:
                    return tuple(reversed(x[:2])) + x[2:]
                else:
                    return x
            else:
                return ck_noop(i - 1, x)
        return f

    return ck_compute_functional(w, a)


def sp_ck_compute_improved(tab, letter):
    read = tab.row_reading_word()
    rows = get_rows(read)
    columns = get_columns(tab.column_reading_word())

    n = len(read)
    r = len(rows)
    word = read + (letter,)

    _print()
    _print('rows:', rows, '<-', letter)
    _print()

    if n == 0:
        return 1, 1, True, (letter,)

    s = [-1]
    for row in rows:
        s = [s[0] + len(row)] + s
    assert s[0] == n - 1

    def rowinsert(i, w):
        for a in range(n - 2, s[i], -1):
            w = ck_noop(a, w)
        return w

    h = [r * (r + 1) // 2]
    for c in columns[r:]:
        h += [h[-1] + len(c)]
    q = len(h)

    def colinsert(i, w):
        for a in range(h[0], h[i] - 1):
            w = ck_noop(a, w)
        return w

    def collectdiag(w):
        for i in range(1, r):
            for a in range(s[i] + i - 1, -1, -1):
                w = ck_noop(a, w)
        for _ in range(r):
            for a in range(-1, r - 1):
                w = ck_noop(a, w)
        return w

    def diagtest(w):
        for a in range(-1, s[-2] - 1):
            w = ck(a, w)
        return not descent(s[-2], w)

    def integrate(w):
        start = r - 1
        for a in range(r - 1):
            for _ in range(a + 1):
                for b in range(start, start - r + 1 + a, -1):
                    w = ck(b, w)
                start += 1
        return w

    _print('* word')
    _print('*', word)
    _print('*', s)
    _print()

    for i in range(r):
        test = rowinsert(i, word)
        if not descent(s[i], test):
            return i + 1, len(rows[-i - 1]) + i + 1, True, test

    test = rowinsert(r - 1, word)
    word = rowinsert(r, word)
    _print('* diagonal bump?')
    _print('*', test)
    _print()

    if diagtest(test):
        return r + 1, r + 1, True, word

    _print('* collecting diagonal')
    _print('*', word)
    _print()

    word = collectdiag(word)

    _print('* left to column word')
    _print('*', word)
    _print()

    _, left = to_column_reading(word[r + 1:])
    word = word[:r + 1] + left

    _print('* integration')
    _print('*', word)
    _print()

    word = integrate(word)

    left = word[r * (r + 1) // 2 + 1:]

    _print('* column insertion')
    _print('*', word, left)
    _print('*', h)
    _print()

    for i in range(q - 1):
        test = colinsert(i, word)
        if descent(h[i], test):
            return len(columns[r + i]) + 1, r + i + 1, False, test

    test = colinsert(q - 1, word)
    return 1, q + r, False, test


def sp_ck_compute(tab, letter):
    ans = []
    rows = get_rows(tab.row_reading_word())

    _print('rows:', rows, '<-', letter)
    _print()

    if len(rows) == 0:
        return 1, 1, True, ans

    left = ()

    for i in range(1, len(rows) + 1):
        o = sum([len(r) for r in rows[:-i]])
        row = rows[-i]

        _print('*', i)
        _print('*', row, ':', letter, ':', left)
        _print()

        working = row + (letter,)
        if not descent(len(row) - 1, working):
            j = len(row) + i
            return i, j, True, ans

        if i == len(rows):
            testing = working
            for a in range(-1, len(row) - 2):
                testing = ck(a, testing)
            if not descent(len(row) - 1, testing):
                ans += [a + o for a in range(len(row) - 2, -1, -1)]
                return i + 1, i + 1, True, ans

        for a in range(len(row) - 2, -1, -1):
            ans += [a + o]
            working = ck(a, working)

        letter = working[0]
        left = working[1:] + left

    working = (letter,) + left

    _print('* collecting diagonal')
    _print('*', working, '|', ans)
    _print()

    for i in range(1, len(rows)):
        for b in range(i):
            s = sum([len(r) for r in rows[-i:]])
            n = len(working) - s - 2 + b
            m = len(working) - s - len(rows[-i - 1]) - 1 + b
            for a in range(n, m, -1):
                ans += [a]
                working = ck(a, working)

    bit = working[:len(rows) + 1]
    left = working[len(rows) + 1:]

    _print('* reversing bit')
    _print('*', bit, ':', left, '|', ans)
    _print()

    for b in range(len(bit) - 2, -1, -1):
        for a in range(-1, b):
            ans += [a]
            bit = ck(a, bit)

    _print('* left to column word')
    _print('*', bit, ':', left, '|', ans)
    _print()

    seq, left = to_column_reading(left)
    ans += [a + len(bit) for a in seq]

    _print('* integration')
    _print('*', bit, ':', left, '|', ans)
    _print()

    m = len(rows) - 1
    working = bit + left
    start = m
    for a in range(m):
        for _ in range(a + 1):
            for b in range(start, start - m + a, -1):
                ans += [b]
                working = ck(b, working)
            start += 1

    i = m * (m + 1) // 2 + len(rows)
    bit, letter, left = working[:i], working[i], working[i + 1:]

    _print('* column insertion')
    _print('*', bit, ':', letter, ':', left, '|', ans)
    _print()

    offset = len(rows)
    columns = get_columns(left)

    if len(columns) == 0:
        ans += to_row_reading(bit + (letter,))[0]
        return 1, offset + 1, False, ans

    for j in range(1, 1 + len(columns)):
        o = i + sum([len(c) for c in columns[:j - 1]])
        col = columns[j - 1]
        working = (letter,) + col

        rest = ()
        for k in range(j, len(columns)):
            rest += columns[k]

        _print('* column', j + offset, j, len(columns))
        _print('*', working)
        _print()

        if descent(0, working):
            i = len(col) + 1
            ans += to_row_reading(bit + working + rest)[0]
            return i, j + offset, False, ans

        if j == len(columns):
            for a in range(len(working) - 2):
                ans += [a + o]
                working = ck(a, working)
            ans += to_row_reading(bit + working + rest)[0]
            return 1, j + offset + 1, False, ans

        for a in range(len(working) - 2):
            ans += [a + o]
            working = ck(a, working)

        letter = working[-1]
        bit += working[:-1]

    raise Exception


def _print(*args):
    pass  # print(*args)


def test_fpf(n=6, maxcount=0):
    for a, pi in enumerate(Permutation.fpf_involutions(n)):
        print(a, pi.fpf_involution_length, pi)
        for b, w in enumerate(pi.get_fpf_involution_words()):
            if len(w) == 0:
                continue
            if b > maxcount > 0:
                break

            v = w[:-1]
            p, q = InsertionAlgorithm.symplectic_hecke(v)
            t, r = InsertionAlgorithm.symplectic_hecke(w)

            print()
            print('word =', w)
            print(p)
            _print(q)
            _print(r)

            newboxes = [b for b in r.boxes if b not in q.boxes]
            assert len(newboxes) == 1

            i, j, sgn, seq = sp_ck_compute(p, w[-1])
            assert (i, j) == newboxes[0]
            assert sgn == (r.get(i, j) > 0)

            _print(seq)

            word = p.row_reading_word() + (w[-1],)
            for a in seq:
                word = ck(a, word)
            assert word == t.row_reading_word()

            a, b, s, word = sp_ck_compute_improved(p, w[-1])
            assert (a, b, s) == (i, j, sgn)
            assert (s and word == t.row_reading_word()) or (not s and word == t.column_reading_word())

            op, Q = sp_ck_compute_functional(w)
            print(Q)
            print(r)
            print(t)
            print(op(w))
            assert Q == r
            assert t.row_reading_word() == op(w)


def test_inv(n=5, maxcount=0):
    for a, pi in enumerate(Permutation.involutions(n)):
        print(a, pi.involution_length, pi)
        for b, w in enumerate(pi.get_involution_words()):
            if len(w) == 0:
                continue
            if b > maxcount > 0:
                break

            v = w[:-1]
            p, q = InsertionAlgorithm.orthogonal_hecke(v)
            t, r = InsertionAlgorithm.orthogonal_hecke(w)

            print()
            print('word =', w)
            _print(p)
            _print(q)
            _print(r)

            op, Q = o_ck_compute_functional(w)
            print(Q)
            print(r)
            print(t)
            print(op(w))
            assert Q == r
            assert t.row_reading_word() == op(w)


def test_standard():
    for mu in Partition.all(6):
        tab = Tableau.standard(mu)
        print(mu)
        for t in tab:
            print(t)
        words = set(Permutation.get_grassmannian(*mu).get_reduced_words())
        assert len(tab) == len(words)
        print()


def test_shifted_standard():
    for mu in Partition.all(6, strict=True):
        tab = Tableau.standard_shifted_marked(mu)
        print(mu)
        for t in tab:
            print(t)
        words = set(Permutation.get_inv_grassmannian(*mu).get_involution_words())
        assert len(tab) == len(words)
        print()


def test_descents():
    for mu in Partition.all(6, strict=True):
        tab = Tableau.standard_shifted_marked(mu)
        print(mu)
        for t in tab:
            w = {a: i for i, a in enumerate(t.shifted_reading_word())}
            d = t.descent_set()
            print(t)
            print(t.shifted_reading_word())
            print(t.descent_set())
            print()
            assert d == {i for i in w if i < sum(mu) and w[i] > w[i + 1]}
        print()


def test_grassmannian():
    w = Permutation.get_grassmannian(4, 4, 3, 2, 2, 2, 1)
    assert w.shape() == (4, 4, 3, 2, 2, 2, 1)

    shapes = set(Partition.subpartitions([4, 3, 2, 1]))
    gr = list(Permutation.grassmannians(5))
    assert len(gr) == len(shapes)
    assert {w.shape() for w in gr} == shapes


def test_inv_grassmannian():
    w = Permutation.get_inv_grassmannian(4, 3, 1)
    t = Permutation.transposition
    assert w == t(3, 7) * t(2, 6) * t(1, 4)
    assert w.involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([5, 3, 1], strict=True))
    gr = list(Permutation.inv_grassmannians(6))
    assert len(gr) == len(shapes)
    assert {w.involution_shape() for w in gr} == shapes


def test_fpf_grassmannian():
    w = Permutation.get_fpf_grassmannian(4, 3, 1)
    t = Permutation.transposition

    assert w == t(1, 5) * t(2, 7) * t(3, 8) * t(4, 6)
    assert w.fpf_involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([6, 4, 2], strict=True))
    gr = list(Permutation.fpf_grassmannians(8))
    assert len(gr) == len(shapes)
    assert {w.fpf_involution_shape() for w in gr} == shapes


def test_eg_insertion():
    rank = 4
    for w in Permutation.grassmannians(rank):
        print(w.shape(), '=', 'shape(', w, ') = shape(', w.inverse(), '^-1 )')
        print()

        n = 0 if len(w) == 0 else list(w.left_descent_set)[0]
        a = tuple(w.inverse()(i) for i in range(n, 0, -1))
        b = tuple(w.inverse()(i) for i in range(n + 1, w.rank + 1))
        mapping = {(x, y): (i + 1, j + 1) for i, x in enumerate(a) for j, y in enumerate(b)}
        print(mapping)
        for word in w.get_reduced_words():

            tab = Tableau()
            partial = [Permutation()]
            for i in reversed(word):
                partial = [partial[0] * Permutation.s_i(i)] + partial
            for i, e in enumerate(word):
                x, y = partial[i](e), partial[i](e + 1)
                x, y = mapping[(x,y)]
                tab = tab.add(x, y, i + 1)

            print(w)
            print(word)
            print(Word.wiring_diagram(word))
            p, q = InsertionAlgorithm.hecke(word)
            print(p)
            print(q)
            print()

            assert q == tab
        print()
        print()
        print()


def double_inv_word(word):
    w = Permutation()
    ans = []
    for i in word:
        s = Permutation.s_i(i)
        if w * s == s * w:
            ans = [None] + ans + [i]
            w = w * s
        else:
            ans = [i] + ans + [i]
            w = s * w * s
    return tuple(ans)


def test_inv_eg_insertion():
    rank = 4
    for w in Permutation.inv_grassmannians(rank):
        w = w.star()
        print(w.involution_shape(), '=', 'shape(', w, ')')
        print()

#        n = 0 if len(w) == 0 else list(w.visible_descent_set)[0]
#        a = tuple(w.inverse()(i) for i in range(n, 0, -1))
#        b = tuple(w.inverse()(i) for i in range(n + 1, w.rank + 1))
#        mapping = {(x, y): (i + 1, j + 1) for i, x in enumerate(a) for j, y in enumerate(b)}

        for word in w.get_involution_words():
            p, q = InsertionAlgorithm.orthogonal_hecke(word)

            # if not q.is_shifted_column_superstandard():
            #     continue

            # tab = Tableau()
            # partial = [Permutation()]
            # for i in reversed(word):
            #     partial = [partial[0] * Permutation.s_i(i)] + partial
            # for i, e in enumerate(word):
            #     x, y = partial[i](e), partial[i](e + 1)
            #     x, y = mapping[(x,y)]
            #     tab = tab.add(x, y, i + 1)

            print()
            print()
            print()
            print(w.cycle_repr())
            print()
            print(w.get_min_atom())
            print()
            print(word)

            labels = {}
            for i in range(1, len(word) + 1):
                j = len(word) + i
                k = len(word) + 1 - i
                if q.find(i):
                    a = q.find(i)[0]
                else:
                    a = tuple(reversed(q.find(-i)[0]))
                labels[j] = a
                labels[k] = ''

            double = double_inv_word(word)
            #print(double)
            #print()
            print(q)
            print(Word.wiring_diagram(double, labels))
            print()
            print(p)
            #print()

        print()
        print()
        print()


def standardize(t):
    seen = set()
    offset = 0
    offsets = {}
    for i, j in sorted(t.boxes, key=lambda box: (box[1], box[0])):
        v = t.get(i, j)
        while v + offset in seen:
            offset += 1
        seen.add(v + offset)
        offsets[(i, j)] = offset
    for i, j in offsets:
        t = t.replace(i, j, t.get(i, j) + offsets[(i, j)])
    return t


def rectify_print(v, printw=True):
    p, q = InsertionAlgorithm.orthogonal_hecke(v)
    t = standardize(p)
    u, _ = InsertionAlgorithm.inverse_orthogonal_hecke(t, q)

    def st(word):
        w = len(word) * [0]
        for i, e in enumerate(sorted(enumerate(word), key=lambda p: (p[1], p[0]))):
            w[e[0]] = i + 1
        return tuple(w)

    if printw:
        labels = get_crossings(v)
        print(Word.involution_wiring_diagram(v, [''] + labels))
        # print(Word.wiring_diagram(u))
        print()
        # print(', '.join(tuple(str(i) + (' ' if i < 10 else '') for i in v)))
        # print(', '.join(tuple(str(i) + (' ' if i < 10 else '') for i in u)))
        # print()
    return u


def fpf_minimal_word(mu):
    return tuple(i + 1 for i in minimal_word(mu))


def minimal_word(mu):
    if len(mu) == 0:
        return ()
    columns = [sum([1 for i in range(len(mu)) if mu[i] + i >= j]) for j in range(len(mu) + 1, mu[0] + 1)]
    word = [i for j in range(len(mu)) for i in range(2 * j + 1, j, -1)]
    m = max(word) + 1
    for c in columns:
        word += [i for i in range(m, m - c, -1)]
        m += 1
    return tuple(word)


def get_crossings(word):
    w = Permutation.from_involution_word(word)

    def endpoint(arrows, index, i):
        key = (index, i)
        while key in arrows:
            key = arrows[key]
        _, a = key
        if w(a) < a:
            a = w(a)
        return a

    n = (max(word) + 1) if word else 0
    arrows = {}
    for index, i in enumerate(word):
        arrows[(index, i)] = (index + 1, i + 1)
        arrows[(index, i + 1)] = (index + 1, i)
        for j in range(1, n + 1):
            if j not in [i, i + 1]:
                arrows[(index, j)] = (index + 1, j)
    ans = []
    for index, i in enumerate(word):
        a = endpoint(arrows, index, i)
        b = endpoint(arrows, index, i + 1)
        ans += [(a, b)]
    return ans


def predict(word):
    z = Permutation.from_involution_word(word)
    minword = minimal_word(z.involution_shape())
    mincrossings = get_crossings(minword)
    values = rectify_print(minword, False)
    base = {}
    for i in range(len(minword)):
        base[mincrossings[i]] = values[i]
        base[tuple(reversed(mincrossings[i]))] = values[i]

    n = z.rank
    fpts = [i for i in range(1, n + 1) if z(i) == i]
    apts = [i for i in range(1, n + 1) if i < z(i)]

    crossings = get_crossings(word)
    keys = sorted({tuple(sorted(key)) for key in crossings})
    index = {}
    for i, key in enumerate(keys):
        index[key] = i + 1
        index[tuple(reversed(key))] = i + 1

    sigma = Permutation()
    for f in reversed(fpts):
        factor = Permutation()
        for j in [ _ for _ in reversed(apts) if f < z(_)]:
            if (f, j) in crossings:
                factor *= Permutation.transposition(index[(f, j)], index[(j, j)])
            elif (j, f) in crossings:
                pass
            else:
                raise Exception
            for i in [_ for _ in reversed(apts) if _ < j and f < z(_)]:
                order = [tuple(sorted(k)) for k in crossings if k in [(i, f), (f, i), (j, f), (f, j)]]
                if order == [(i, f), (j, f)]:
                    factor *= Permutation.cycle(index[(i, j)], index[(i, f)], index[(j, f)])
                elif order == [(j, f), (i, f)]:
                    pass
                else:
                    raise Exception
        sigma = factor * sigma
    return tuple(base[keys[sigma(index[c]) - 1]] for c in crossings)


def columns_subwords(p):
    columns = []
    for i, j in sorted(p.boxes, key=lambda x: (x[1], -x[0])):
        if j > len(columns):
            columns += [[]]
        columns[-1].append(p.get(i, j))
    return columns


@pytest.mark.slow
def test_inv_predict(rank=8, bound=0):
    delta = tuple(range(rank - 1, 0, -2))
    partitions = list(Partition.subpartitions(delta, strict=True))
    for count, mu in enumerate(partitions):
        w = Permutation.get_inv_grassmannian(*mu).star()
        p, _ = InsertionAlgorithm.orthogonal_hecke(w.get_involution_word())
        columns = columns_subwords(p)
        stcolumns = columns_subwords(standardize(p))
        print(p)
        print(columns)
        print(stcolumns)
        print()

        nodes = []
        for i in range(len(columns) - 1, -1, -1):
            while columns[i] != stcolumns[i]:
                word = [a for j, col in enumerate(columns) for a in (col[:-1] if i == j else col)]
                nodes.append(Permutation.from_involution_word(word))
                columns[i] = [a + 1 for a in columns[i]]
        assert all(y is not None for y in nodes)

        for index, word in enumerate(w.get_involution_words()):
            if index > bound > 0:
                break
            testpr = word
            for y in nodes:
                testpr = Permutation.involution_little_bump(testpr, y)
            pr = predict(word)
            print(len(partitions) - count, ':', w, word, ':', index)
            print(' test =', testpr)
            print(' real =', pr)
            print()
            assert testpr == pr


@pytest.mark.slow
def test_inv_grassmannian_braids(rank=7):
    delta = tuple(range(rank - 1, 0, -2))
    rect = {}
    for mu in Partition.subpartitions(delta, strict=True):
        # w = (1, 3, 2, 7, 9, 8, 4, 6, 5, 6, 4, 3, 7, 4, 5, 6, 5, 4)
        # z = Permutation.from_involution_word(w)
        # mu = z.involution_shape()

        minword = minimal_word(mu)
        rect[minword] = rectify_print(minword, False)

        def getmap(word):
            return {
                tuple(sorted(x)): rect[word][i]
                for i, x in enumerate(get_crossings(word))
            }

        mmap = getmap(minword)

        w = minword
        # w = (1, 3, 2, 7, 9, 8, 4, 6, 5, 6, 4, 3, 7, 4, 5, 6, 5, 4)

        seen = set()
        add = {w}
        while True:
            nextadd = set()
            for v in add:
                rect[v] = rectify_print(v, False)
                vmap = getmap(v)

                for i, u in get_inv_ck_moves(v):
                    r = rect[v]
                    if i >= 0 and v[i] == v[i + 2] < v[i + 1]:
                        assert r[i + 2] < r[i] < r[i + 1]
                    if i >= 0 and v[i] == v[i + 2] > v[i + 1]:
                        assert r[i + 1] < r[i] < r[i + 2]
                    if i >= 0 and v[i + 1] < v[i] < v[i + 2]:
                        assert r[i + 1] < r[i] < r[i + 2]
                    if i >= 0 and v[i + 2] < v[i] < v[i + 1]:
                        assert r[i + 2] < r[i] < r[i + 1]
                    if i >= 0 and v[i] < v[i + 2] < v[i + 1]:
                        assert r[i] < r[i + 2] < r[i + 1]
                    if i >= 0 and v[i + 1] < v[i + 2] < v[i]:
                        assert r[i + 1] < r[i + 2] < r[i]

                    if u not in seen:
                        rect[u] = rectify_print(u, False)
                        umap = getmap(u)
                        nextadd.add(u)

                        if len([key for key in mmap if vmap[key] != umap[key]]) >= 0:
                            rectify_print(minword, True)
                            rectify_print(v, True)
                            rectify_print(u, True)

                            print(' ' + i * '   ' + '*')
                            print(v, '          ', rect[v])
                            print(u, '          ', rect[u])
                            print()

                            for key in mmap:
                                print(key, ':', mmap[key], '=>', vmap[key], '=>', umap[key], '*' if vmap[key] != umap[key] else '')

                            print(10 * '\n')

                    pv = predict(v)
                    rectify_print(minword, True)
                    rectify_print(v, True)
                    print(v)
                    print(rect[v])
                    print(pv, 'predicted')
                    print()
                    assert rect[v] == pv

            if len(nextadd) == 0:
                break
            seen |= add
            add = nextadd

    print('\ndone')


def fpf_rectify_print(v, printw=True):
    p, q = InsertionAlgorithm.symplectic_hecke(v)
    t = standardize(p).apply(lambda x: x - 1).double()
    u, _ = InsertionAlgorithm.inverse_symplectic_hecke(t, q)
    assert all(i % 2 == 0 for i in u)
    u = tuple(i // 2 for i in u)
    if printw:
        print(Word.fpf_involution_wiring_diagram(v))
        print(Word.involution_wiring_diagram(u))
        print()
    return u


def fpf_crossings(word):
    n = (max(word) + 1) if word else 0
    n = n + 1 if n % 2 != 0 else n
    w = Permutation.from_fpf_involution_word(word)

    v = Permutation()
    for i in range(1, n + 1):
        if i < w(i) and not any(j < w(j) for j in range(i + 1, w(i))):
            v *= Permutation.transposition(i, w(i))
    w = w * v

    def endpoint(arrows, index, i):
        key = (index, i)
        while key in arrows:
            key = arrows[key]
        _, a = key
        if w(a) < a:
            a = w(a)
        return a

    arrows = {}
    for index, i in enumerate(word):
        arrows[(index, i)] = (index + 1, i + 1)
        arrows[(index, i + 1)] = (index + 1, i)
        for j in range(1, n + 1):
            if j not in [i, i + 1]:
                arrows[(index, j)] = (index + 1, j)
    ans = []
    for index, i in enumerate(word):
        a = endpoint(arrows, index, i)
        b = endpoint(arrows, index, i + 1)
        ans += [(a, b)]
    return ans


def fpf_predict(word):
    z = Permutation.from_fpf_involution_word(word)
    minword = fpf_minimal_word(z.fpf_involution_shape())
    mincrossings = fpf_crossings(minword)
    values = fpf_rectify_print(minword, False)
    base = {}
    for i in range(len(minword)):
        base[mincrossings[i]] = values[i]
        base[tuple(reversed(mincrossings[i]))] = values[i]

    n = (max(word) + 1) if word else 0
    n = n + 1 if n % 2 != 0 else n
    v = Permutation()
    for i in range(1, n + 1):
        if i < z(i) and not any(j < z(j) for j in range(i + 1, z(i))):
            v *= Permutation.transposition(i, z(i))
    z = z * v

    fpts = [i for i in range(1, n + 1) if z(i) == i]
    apts = [i for i in range(1, n + 1) if i < z(i)]

    crossings = fpf_crossings(word)
    keys = sorted({tuple(sorted(key)) for key in crossings})
    index = {}
    for i, key in enumerate(keys):
        index[key] = i + 1
        index[tuple(reversed(key))] = i + 1

    sigma = Permutation()
    for f in reversed(fpts):
        factor = Permutation()
        for j in [_ for _ in reversed(apts) if f < z(_)]:
            for i in [_ for _ in reversed(apts) if _ < j and f < z(_)]:
                order = [tuple(sorted(k)) for k in crossings if k in [(i, f), (f, i), (j, f), (f, j)]]
                if order == [(i, f), (j, f)]:
                    factor *= Permutation.cycle(index[(i, j)], index[(j, f)], index[(i, f)])
                elif order == [(j, f), (i, f)]:
                    pass
                else:
                    raise Exception
        sigma = factor * sigma
    return tuple(base[keys[sigma(index[c]) - 1]] for c in crossings)


@pytest.mark.slow
def test_fpf_grassmannian_braids(rank=8):
    delta = tuple(range(rank - 2, 0, -2))
    rect = {}
    for mu in Partition.subpartitions(delta, strict=True):

        minword = fpf_minimal_word(mu)
        rect[minword] = fpf_rectify_print(minword, False)
        w = minword

        def getmap(word):
            return {
                tuple(sorted(x)): rect[word][i]
                for i, x in enumerate(fpf_crossings(word))
            }

        mmap = getmap(minword)

        p, q = InsertionAlgorithm.symplectic_hecke(w)

        seen = set()
        add = {w}
        while True:
            nextadd = set()
            for v in add:
                rect[v] = fpf_rectify_print(v, False)
                vmap = getmap(v)

                for i, u in get_fpf_ck_moves(v):
                    if u not in seen:
                        rect[u] = fpf_rectify_print(u, False)
                        umap = getmap(u)
                        nextadd.add(u)

                        # if (i == -1 and abs(u[0] - u[1]) == 1) or
                        # if (i >= 0 and u[i] == u[i + 2]):
                        if any(vmap[key] != umap[key] for key in vmap) and not (i >= 0 and u[i] == u[i + 2]):
                            labels = [''] + fpf_crossings(minword)
                            print(Word.fpf_involution_wiring_diagram(minword, labels))
                            print(rect[minword])
                            print()
                            labels = [''] + fpf_crossings(v)
                            print(Word.fpf_involution_wiring_diagram(v, labels))
                            print(2 * '\n')
                            labels = [''] + fpf_crossings(u)
                            print(Word.fpf_involution_wiring_diagram(u, labels))

                            # vv = rectify_print(rect[v], False)
                            # uu = rectify_print(rect[u], False)
                            # if rect[v] == vv and rect[u] == uu:
                            #    continue

                            print(' ' + i * '   ' + '*')
                            print(v, '          ', rect[v])
                            print(u, '          ', rect[u],)
                            print()

                            for key in mmap:
                                print(key, ':', mmap[key], '=>', vmap[key], '=>', umap[key], '*' if vmap[key] != umap[key] else '')

                            print(2 * '\n')
                            # input(str(p))

                        r = rect[v]
                        s = rect[u]

                        if i >= 0:
                            assert r[:i] == s[:i]
                            assert r[i + 3:] == s[i + 3:]
                        elif i == -1:
                            assert r[2:] == s[2:]

                        if i >= 0 and v[i] == v[i + 2] < v[i + 1]:
                            assert r[i] < r[i + 2] < r[i + 1]
                        if i >= 0 and v[i] == v[i + 2] > v[i + 1]:
                            assert r[i + 1] < r[i + 2] < r[i]

                        if i >= 0 and v[i + 1] < v[i] < v[i + 2]:
                            assert r[i + 1] < r[i] < r[i + 2]
                        if i >= 0 and v[i + 2] < v[i] < v[i + 1]:
                            assert r[i + 2] < r[i] < r[i + 1]

                        if i >= 0 and v[i] < v[i + 2] < v[i + 1]:
                            assert r[i] < r[i + 2] < r[i + 1]

                        if i >= 0 and v[i + 1] < v[i + 2] < v[i]:
                            assert r[i + 1] < r[i + 2] < r[i]

                        if i == -1 and v[0] < v[1]:
                            assert r[0] < r[1]
                        if i == -1 and v[1] < v[0]:
                            assert r[1] < r[0]

                pv = fpf_predict(v)
                if rect[v] != pv:
                    labels = [''] + fpf_crossings(minword)
                    print(Word.fpf_involution_wiring_diagram(minword, labels))
                    labels = [''] + fpf_crossings(v)
                    print(Word.fpf_involution_wiring_diagram(v, labels))
                    print(v)
                    print(rect[v])
                    print(pv, 'predicted')
                    print()
                    for key in mmap:
                        print(key, ':', mmap[key], '=>', vmap[key], '*' if vmap[key] != mmap[key] else '')
                    print()
                assert rect[v] == pv

            if len(nextadd) == 0:
                # input('?')
                break
            seen |= add
            add = nextadd

    print('\ndone')


def test_inv_conversion_to_permutation(rank=5):
    def pi(mu):
        r = len(mu)
        ans = Permutation()
        for i in range(1, r + 1):
            ans *= Permutation.transposition(i, r + mu[-i])
        return ans

    def bumping_sequence(mu):
        assert len(mu) == 0 or mu[-1] > 0
        r = len(mu)
        if r == 0:
            return tuple()
        q = mu[0]

        h = 1
        while h + 1 <= r and mu[h] + h == q:
            h += 1

        nu = list(mu)
        nu[h - 1] -= 1
        nu = tuple(nu)
        while nu and nu[-1] == 0:
            nu = nu[:-1]

        n = sum(mu)
        s = Permutation.s_i(2 * n)
        ans = [s * w for w in bumping_sequence(nu)]
        bns = (2 * n - q - r + 1) * [pi(nu)]
        return tuple(ans + bns)

    for w in Permutation.inv_grassmannians(rank):
        mu = w.involution_shape()
        print(mu, '=', 'shape(', w, ')')
        print()

        seq = bumping_sequence(mu)
        print('  bumping_sequence:')
        for v in seq:
            print('    ', v)
        print()

        print()
        word = minimal_word(mu)
        for sigma in reversed(seq):
            print('  ', sigma.get_involution_word(), 'bumps', word)
            word = Permutation.involution_little_bump(word, sigma)
        print('  result:', word)
        print()

        word = w.get_involution_word()
        for sigma in reversed(seq):
            word = Permutation.involution_little_bump(word, sigma)
        print('  result:', word)
        print()
        assert set(word) == {2 * i for i in range(1, sum(mu) + 1)}


def test_fpf_conversion_to_permutation(rank=4):
    def pi(mu):
        r = len(mu)
        ans = Permutation()
        for i in range(1, r + 1):
            ans *= Permutation.transposition(i, r + 1 + mu[-i])
        n = ans.rank
        f = [i for i in range(1, n + 1) if ans(i) == i]
        if len(f) % 2 != 0:
            f += [n + 1]
        for i in range(0, len(f), 2):
            ans *= Permutation.transposition(f[i], f[i + 1])
        return ans

    def bumping_sequence(mu):
        assert len(mu) == 0 or mu[-1] > 0
        r = len(mu)
        if r == 0:
            return tuple()
        q = mu[0]

        h = 1
        while h + 1 <= r and mu[h] + h == q:
            h += 1

        nu = list(mu)
        nu[h - 1] -= 1
        nu = tuple(nu)
        while nu and nu[-1] == 0:
            nu = nu[:-1]

        e = 1 if (q + r) % 2 == 0 else 0
        n = sum(mu)
        s = Permutation.s_i(2 * n)
        ans = []
        for w in bumping_sequence(nu):
            i = w.rank + 1
            assert i % 2 != 0
            while i <= 2 * n + 1:
                w *= Permutation.s_i(i)
                i += 2
            ans += [s * w * s]
        bns = (n - (q + r + e + 1) // 2 + 1) * [pi(nu)]
        return tuple(ans + bns)

    for w in Permutation.fpf_grassmannians(rank):
        w = w.star()
        mu = w.fpf_involution_shape()
        print(mu, '=', 'shape(', w, ')')
        print()

        seq = bumping_sequence(mu)
        print('  bumping_sequence:')
        for v in seq:
            print('    ', v)
        print()

        print()
        word = fpf_minimal_word(mu)
        for sigma in reversed(seq):
            newword = Permutation.fpf_involution_little_bump(word, sigma)
            print('  ', sigma, 'bumps', word, 'to', newword)
            word = newword
        print()
        print('  result:', word)
        print()

        word = w.get_fpf_involution_word()
        for sigma in reversed(seq):
            word = Permutation.fpf_involution_little_bump(word, sigma)
        print('  result:', word)
        print()
        assert set(word) == {2 * i for i in range(1, sum(mu) + 1)}
