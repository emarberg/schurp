from permutations import Permutation
from tableaux import Tableau
from marked import MarkedNumber
from words import (
    Word,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
    involution_insert,
    primed_sw_insert,
    sagan_worley_insert
)
import pytest
import random

sw_cache = {}
pq_cache = {}
commutations_cache = {}
insert_cache = {}
full_cseq_cache = {}
tau_cache = {}


def primed_cached_insert(a):
    p, q = cached_insert(a)
    for i in range(1, 1 + p.max_row):
        if p.entry(i, i).is_primed():
            p = p.set(i, i, -p.entry(i, i))
            q = q.set(i, i, -q.entry(i, i))
    return p, q


def cached_insert(a):
    if a not in pq_cache:
        pq_cache[a] = Word(*a).involution_insert()
    return pq_cache[a]


def cached_sw_insert(a):
    if a not in sw_cache:
        sw_cache[a] = Word(*a).primed_sw_insert()
    return sw_cache[a]


def matrix_to_biword(mat):
    ans = []
    for i in range(len(mat)):
        word = []
        for j in range(len(mat[i])):
            if mat[i][j] < 0:
                word += [-(j + 1)] + (-1 - mat[i][j]) * [j + 1]
            elif mat[i][j] > 0:
                word += mat[i][j] * [j + 1]
        ans += [Word(*word)]
    return ans


def integer_matrices(m, n, p):
    r = 2 * p + 1
    for v in range(r ** (m * n)):
        ans = [n * [0] for _ in range(m)]
        for i in range(m):
            for j in range(n):
                ans[i][j] = (v % r) - p
                v = v // r
        yield ans


def random_integer_matrix(m, n, p):
    ans = [n * [0] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            ans[i][j] = random.randint(-p, p)
    return ans


def sw_permutation(biword):
    ans = Permutation()
    word = tuple(a for w in biword for a in w)
    for i in range(len(word)):
        u = cached_sw_insert(word[:i])[0].get_rows()
        v = cached_sw_insert(word[:i + 1])[0].get_rows()
        if len(u) == len(v):
            q = len(u)
            t = Permutation()
            for i in range(q):
                a, b = u[i][0], v[i][0]
                if a != b:
                    a, b = abs(a), abs(b)
                    assert a != b
                    assert t == Permutation()
                    t = Permutation.transposition(a, b)
            ans = ans * t

    signs = {}
    for a in word:
        if abs(a) not in signs:
            signs[abs(a)] = -1 if a < 0 else 1
    return ans, signs, word


def help_test_primed_sw_insertion(biword, mapping={}):
    pi, signs, word = sw_permutation(biword)

    def func(i):
        return -i if signs[abs(i)] != signs[pi(abs(i))] else i

    def transform(tab, rec):
        tab = Tableau({
            box: -val if box[0] == box[1] and rec.entry(box[0], box[1]).is_primed() else val
            for box, val in tab.mapping.items()
        })

        return Tableau({
            box: val if not tab.is_togglable(box) else func(val)
            for box, val in tab.mapping.items()
        })

    try:
        p, q = primed_sw_insert(*biword)
        pp, qq = sagan_worley_insert(*biword)
        assert p.is_shifted_semistandard()
        assert q.is_shifted_semistandard(False)
        assert (p, q) not in mapping or mapping[(p, q)] == biword
        assert q.unprime_diagonal() == qq
        assert p.unprime_togglable() == pp.unprime_togglable()
        if pp != transform(p, q):
            print(biword)
            print(word)
            print(signs)
            print(pi)
            print(p)
            print(pp)
            print(transform(p, q))
            input('\n\n\n')
        mapping[(p, q)] = biword

        # print()
        # print('\n'.join([str(row) for row in mat]))
        # print()
        # print(biword)
        # print()
        # print(p)
        # print(q)
        # print()
        # print(pp)
        # print(qq)
    except:
        # print('\n'.join([str(row) for row in mat]))
        # print()
        print(biword)
        print()
        print(p)
        print(q)
        print(mapping.get((p, q), None))
        print()
        print(pp)
        print(qq)
        assert False


def test_primed_sw_insertion_simple():
    biword = [Word(3,), Word(-2,)]
    help_test_primed_sw_insertion(biword)

    biword = [Word(1, 3), Word(2, 6), Word(-5), Word(4), Word(3)]
    help_test_primed_sw_insertion(biword)

    biword = [Word(2, 3, 4), Word(-1, 4, 4, 5), Word(2, 2, -4, 4, 5, 5), Word(-1, 2, -4, 5, 5), Word(2, 3, -4, 4)]
    help_test_primed_sw_insertion(biword)


def test_random_primed_sw_insertion(m=5, n=5, bound=5):
    total = 1000
    mapping = {}
    iteration = 0
    for _ in range(total):
        mat = random_integer_matrix(m, n, bound)
        biword = matrix_to_biword(mat)
        print(iteration, 'of', total)
        help_test_primed_sw_insertion(biword, mapping)
        iteration += 1


def test_primed_sw_insertion(m=2, n=3, bound=2):
    total = (2 * bound + 1) ** (m * n)
    mapping = {}
    iteration = 0
    for mat in integer_matrices(m, n, bound):
        if iteration % 100 == 0:
            print(iteration, 'of', total)
        biword = matrix_to_biword(mat)
        help_test_primed_sw_insertion(biword, mapping)
        iteration += 1


def test_standard_primed_sw_insertion(n=3, bound=2):
    total = (2 * bound) ** n
    mapping = {}
    iteration = 0
    for v in range(total):
        if iteration % 100 == 0:
            print(iteration, 'of', total)
        try:
            word = []
            for _ in range(n):
                x = v % (2 * bound)
                word += [x + 1 if x < bound else -(x + 1 - bound)]
                v = v // (2 * bound)

            p, q = Word(*word).primed_sw_insert()
            assert p.is_shifted_semistandard()
            assert q.is_shifted_standard(False)
            assert (p, q) not in mapping
            mapping[(p, q)] = word

            for i in range(len(word) - 2):
                a, c, b = word[i:i + 3]
                ma, mc, mb = MarkedNumber(a), MarkedNumber(c), MarkedNumber(b)
                if ma.ceil() <= mb <= mc.floor():
                    vord = word[:i] + [c, a, b] + word[i + 3:]
                    pp = Word(*vord).primed_sw_insert()[0]
                    # print('* ', i, ':', a, c, b)
                    assert pp == p
                b, c, a = word[i:i + 3]
                ma, mc, mb = MarkedNumber(a), MarkedNumber(c), MarkedNumber(b)
                if ma <= mb.floor() and mb.ceil() <= mc:
                    vord = word[:i] + [b, a, c] + word[i + 3:]
                    pp = Word(*vord).primed_sw_insert()[0]
                    # print('* ', i, ':', a, c, b)
                    assert pp == p
            if len(word) >= 2:
                a, b = word[:2]
                if a * b > 0:
                    vord = [b, a] + word[2:]
                else:
                    vord = [-b, -a] + word[2:]
                pp = Word(*vord).primed_sw_insert()[0]
                assert pp == p

            # if iteration % 1234 == 0:
            #     print()
            #     print(word)
            #     print()
            #     print(p)
            #     print(q)
        except:
            print(word)
            print()
            print(p)
            print(q)
            print(mapping.get((p, q), None))
            print(pp)
            print(vord)
            assert False
        iteration += 1


def entries(tab):
    assert type(tab) == Tableau
    return {_.number for _ in tab.entries()}


def entries_diagonal(tab):
    assert type(tab) == Tableau
    return {tab.get(x, y).number for (x, y) in tab if x == y}


def commutations(iword):
    iword = tuple(iword)
    if iword not in commutations_cache:
        ans = {}
        w = Permutation()
        for index, i in enumerate(iword):
            i = abs(i)
            s = Permutation.s_i(i)
            if w(i) == i and w(i + 1) == i + 1:
                ans[index] = (i, i + 1)
                w = w * s
            else:
                ans = {k: (s(v[0]), s(v[1])) for (k, v) in ans.items()}
                w = s * w * s
        commutations_cache[iword] = ans
    return commutations_cache[iword]


def marked_cycles(iword):
    ans = {}
    w = Permutation()
    for i in iword:
        sgn = -1 if i < 0 else 1
        i = abs(i)
        s = Permutation.s_i(i)
        if w(i) == i and w(i + 1) == i + 1:
            ans[(i, i + 1)] = sgn
            w = w * s
        else:
            assert sgn == 1
            ans = {(s(a), s(b)): ans[(a, b)] for (a, b) in ans}
            w = s * w * s
    return ans


def shifted_rows_to_columns(rows):
    return [[rows[i][j - i] for i in range(len(rows)) if 0 <= j - i < len(rows[i])] for j in range(len(rows[0]) if rows else 0)]


def columns_to_shifted_rows(cols):
    m = max([0] + [len(c) for c in cols])
    return [[cols[j][i] for j in range(len(cols)) if i < len(cols[j])] for i in range(m)]


def extract_cycle(a, rows, rest):
    p = {a, a + 1}
    for row in reversed(rows):
        for i in row:
            s = Permutation.s_i(i)
            p = {s(x) for x in p}
    for i in rest:
        s = Permutation.s_i(abs(i))
        p = {s(x) for x in p}
    return tuple(sorted(p))


def extract_commutation_positions(rows):
    ans = {}
    w = Permutation()
    for i in range(len(rows) - 1, -1, -1):
        for j in range(len(rows[i])):
            a = rows[i][j]
            s = Permutation.s_i(a)
            if w(a) == a and w(a + 1) == a + 1:
                w = w * s
                ans[(a, a + 1)] = (i, j)
            else:
                ans = {(s(a), s(b)): ans[(a, b)] for (a, b) in ans}
                w = s * w * s
    return ans


def insert_old_version(iword):
    def apply(s, pair):
        return tuple(sorted([s(i) for i in pair]))

    def update_cycles(a, i, j, rows, rest, permuted_cycles):
        if a + 1 < rows[i][j]:
            p = extract_cycle(rows[i][j], rows[:i + 1], rest)
            q = extract_cycle(a, rows[:i], rest)
            permuted_cycles[p], permuted_cycles[q] = permuted_cycles[q], permuted_cycles[p]
        if a == rows[i][j]:
            p = extract_cycle(a, rows[:i] + [rows[i] + [a]], rest)
            q = extract_cycle(a + 2, rows[:i] + [rows[i] + [a]] + [rows[i + 1]], rest)
            permuted_cycles[p], permuted_cycles[q] = permuted_cycles[q], permuted_cycles[p]

    rows = []
    permuted_cycles = {p: p for p in marked_cycles(iword)}
    for n, a in enumerate(iword):
        rest = iword[n + 1:]
        transposed = False
        a = abs(a)
        i = 0
        while i < len(rows):
            j = [t for t in range(len(rows[i])) if a <= rows[i][t]]
            if j:
                j = j[0]
                ####
                if j == 0 and not transposed:
                    update_cycles(a, i, j, rows, rest, permuted_cycles)
                ####
                if rows[i][j] == a:
                    a = a + 1
                else:
                    a, rows[i][j] = rows[i][j], a
                if j == 0 and not transposed:
                    rows = shifted_rows_to_columns(rows)
                    transposed = True
            else:
                rows[i] += [a]
                a = None
                break
            i += 1
        if a is not None:
            rows += [[a]]
        if transposed:
            rows = columns_to_shifted_rows(rows)
    diagonal_cycles = {extract_cycle(rows[i][0], rows[:i + 1], []) for i in range(len(rows))}
    cycle_map = extract_commutation_positions(rows)
    return rows, diagonal_cycles, permuted_cycles, cycle_map


def primed_insert(a):
    rows, diagonal_cycles, permuted_cycles, cycle_map, tabs, paths = insert(a)
    marks = marked_cycles(a)
    for box in permuted_cycles:
        sgn = marks[permuted_cycles[box]]
        i, j = cycle_map[box]
        rows[i][j] *= sgn
    for box in diagonal_cycles:
        i, j = cycle_map[box]
        rows[i][j] = abs(rows[i][j])
    return Tableau.shifted_from_rows(rows), tabs, paths


def test_primed_insertion(n=5):
    for w in Permutation.involutions(n):
        for a in w.get_primed_involution_words():
            tab, inter, paths = primed_insert(a)

            v = tuple(Word(_) for _ in a)
            p, q = involution_insert(*v)

            # rows, diagonal_cycles, permuted_cycles, cycle_map = insert(a)
            # print('a =', a)
            # print(marked_cycles(a))
            # print(diagonal_cycles, permuted_cycles, cycle_map)
            # print()
            # print(tab)
            # print(p)
            # print('\n\n')
            assert p == tab


def batch_insert_sequential(tab, a):
    tab = tab if type(tab) == Tableau else Tableau.shifted_from_rows(tab)
    word = tab.column_reading_word() + tuple(a)
    ans, _, _ = primed_insert(word)
    return ans


def row_insert(row, a):
    b = []
    row = row[:]
    for i in a:
        if not row or abs(row[-1]) < abs(i):
            row.append(i)
        else:
            t = [j for j in range(len(row)) if abs(i) <= abs(row[j])][0]
            j = row[t]
            if -i == j < 0:
                row[t + 1] *= -1
                row[t] *= -1
                b += [i + 1]
            elif -i == j > 0:
                b += [i - 1]
            elif i == j:
                assert i > 0
                b += [i + 1]
            else:
                row[t] = i
                b += [j]
    return b, row


def batch_row_insert(tab, a):
    rows = tab.get_rows()
    i = 0
    while i < len(rows) and a and min(map(abs, a)) > min(map(abs, rows[i])):
        a, rows[i] = row_insert(rows[i], a)
        i += 1
    top = Tableau.shifted_from_rows(rows[i:])
    bot = Tableau.shifted_from_rows(rows[:i])
    return top, a, bot


def batch_combine(top, bot):
    top = top.get_rows()
    bot = bot.get_rows()
    r = len(bot)
    top = shifted_rows_to_columns(top)
    bot = shifted_rows_to_columns(bot)
    ans, bot = bot[:r], bot[r:]
    q = max(len(top), len(bot)) + 1
    for i in range(q):
        tcol = top[i] if i < len(top) else []
        bcol = bot[i] if i < len(bot) else []
        if len(tcol) == len(bcol) == 0:
            continue
        elif len(tcol) == 0:
            ans += [bcol]
        elif len(bcol) == 0:
            ans += [tcol]
        else:
            for j in range(len(tcol)):
                t = i
                # print('*', i, j, top, bot)
                while t < len(top) and abs(top[t][j]) <= abs(bcol[-1]):
                    t += 1
                a = [top[x][j] for x in range(i, t)]
                a = list(reversed(a))
                if t < len(top):
                    bcol = bcol + [top[t][j]]
                    a, bcol = row_insert(bcol, a)
                    a = list(reversed(a))
                    for x in range(i + 1, t + 1):
                        top[x][j] = a[x - i - 1]
                else:
                    a, bcol = row_insert(bcol, a)
                    a = list(reversed(a))
                    for x in range(i + 1, t):
                        top[x][j] = a[x - i - 1]
                    top += [[a[-1]]]

            ans += [bcol]
    ans = columns_to_shifted_rows(ans)
    return Tableau.shifted_from_rows(ans)


def batch_insert_parallel(tab, a):
    tab = tab if type(tab) == Tableau else Tableau.shifted_from_rows(tab)
    top, a, bot = batch_row_insert(tab, a)
    top = batch_insert_sequential(top, a)
    return batch_combine(top, bot)


# def _test_batch_insert(n=5):
#     for w in Permutation.involutions(n):
#         for word in w.get_primed_involution_words():
#             print(w, word)
#             p, _, _ = primed_insert(word)
#             for k in range(len(word)):
#                 tab, _, _ = primed_insert(word[:k])
#                 a = word[k:]
#                 q = batch_insert_sequential(tab, a)
#                 assert p == q
#                 r = batch_insert_parallel(tab, a)
#                 assert q == r


def extract_all_cycles(rows, rest):
    ans = []
    for i in range(len(rows)):
        ans.append(extract_cycle(rows[i][0], rows[:i + 1], rest))
    return ans


def insert(iword):
    iword = tuple(iword)
    if iword not in insert_cache:
        rows = []
        cycle_sequence = [[]]
        tabs = []
        paths = []
        for n, a in enumerate(iword):
            rest = iword[n + 1:]
            transposed = False
            a = abs(a)
            i = 0
            path = []
            while i < len(rows):
                j = [t for t in range(len(rows[i])) if a <= rows[i][t]]
                if j:
                    j = j[0]
                    b = rows[i][j] + int(rows[i][j] == a)
                    ###
                    if transposed:
                        x, y = j + 1, i + 1
                        x2, y2 = (x, y) if rows[i][j] != a else (x + 1, y)
                        path += [(x, y, x2, y2, a, b)]
                    else:
                        x, y = i + 1, i + j + 1
                        x2, y2 = (x, y) if rows[i][j] != a else (x, y + 1)
                        path += [(x, y, x2, y2, a, b)]
                    ###
                    if rows[i][j] == a:
                        a = a + 1
                    else:
                        a, rows[i][j] = rows[i][j], a
                    if j == 0 and not transposed:
                        rows = shifted_rows_to_columns(rows)
                        transposed = True
                else:
                    rows[i] += [a]
                    ###
                    if transposed:
                        x, y = len(rows[i]), i + 1
                        path += [(x, y, x, y, a, None)]
                    else:
                        x, y = i + 1, i + len(rows[i])
                        path += [(x, y, x, y, a, None)]
                    ###
                    a = None
                    break
                i += 1
            if a is not None:
                rows += [[a]]
                ###
                if transposed:
                    x, y = 1, len(rows)
                    path += [(x, y, x, y, a, None)]
                else:
                    x, y = len(rows), len(rows)
                    path += [(x, y, x, y, a, None)]
                ###
            if transposed:
                rows = columns_to_shifted_rows(rows)
            if a is not None or transposed:
                cycle_sequence.append(extract_all_cycles(rows, rest))
            paths.append(path)
            tabs.append(Tableau.shifted_from_rows(rows))

        permuted_cycles = {p: p for p in marked_cycles(iword)}
        for i in range(1, len(cycle_sequence)):
            a, b = cycle_sequence[i - 1], cycle_sequence[i]
            rng = [t for t in range(len(a)) if a[t] != b[t]]
            pq = {a[j] for j in rng} | {b[j] for j in rng}
            if pq:
                assert len(pq) == 2
                p, q = tuple(pq)
                permuted_cycles[p], permuted_cycles[q] = permuted_cycles[q], permuted_cycles[p]

        diagonal_cycles = {extract_cycle(rows[i][0], rows[:i + 1], []) for i in range(len(rows))}
        cycle_map = extract_commutation_positions(rows)
        insert_cache[iword] = rows, diagonal_cycles, permuted_cycles, cycle_map, tabs, paths
    return insert_cache[iword]


def _test_bumping_path(w):
    def mid(path):
        m = [i + 1 for i, a in enumerate(path) if a[0] == a[1]]
        return m[0] if m else len(path)

    def is_row_bumped(path, index):
        return index + 1 <= mid(path)

    _, _, _, _, tab, pat = insert(w)
    for i in range(len(w)):
        p = pat[i]
        j = mid(p)
        try:
            for t in range(len(p)):
                x, y, xx, yy, _, _ = p[t]
                if (x, y) != (xx, yy) and xx == yy:
                    assert t == j
                if (x, y) == (xx, yy) and xx == yy:
                    assert t == j - 1

            for t in range(j - 1):
                x1, y1, xx1, yy1, _, _ = p[t]
                x2, y2, xx2, yy2, _, _ = p[t + 1]
                assert x1 == xx1 == t + 1 and x2 == xx2 == t + 2
                assert y1 >= y2 >= j and yy1 >= yy2 >= j
            for t in range(j, len(p) - 1):
                x1, y1, xx1, yy1, _, _ = p[t]
                x2, y2, xx2, yy2, _, _ = p[t + 1]
                assert j >= x1 >= x2 and j + 1 >= xx1 >= xx2
                assert y1 == yy1 == t + 1 and y2 == yy2 == t + 2
            for s in range(j):
                for t in range(j, len(p)):
                    rx, ry, rxx, ryy, _, _ = p[s]
                    cx, cy, cxx, cyy, _, _ = p[t]

                    assert rx < cx or ry < cy
                    assert rxx < cxx or ryy < cyy

                    if s != j - 1 or t != j or ry == ryy:
                        assert rxx < cx or ryy < cy
                    else:
                        assert ry == rx == rxx == cx == j
                        assert cyy == cxx == cy == ryy == j + 1

        except:
            print(tab[i - 1])
            print(tab[i])
            print(p, mid(p))
            assert False

    for i in range(len(w) - 1):
        a, b = w[i:i + 2]
        p, q = pat[i:i + 2]
        m = min(len(p), len(q))
        try:
            if a < b:
                for j in range(m):
                    if is_row_bumped(p, j):
                        px, py, pxx, pyy, _, _ = p[j]
                        qx, qy, qxx, qyy, _, _ = q[j]
                        assert is_row_bumped(q, j)
                        assert px == qx == pxx == qxx == j + 1
                        assert py < qy and pyy < qyy

                if is_row_bumped(p, len(p) - 1):
                    assert is_row_bumped(q, len(q) - 1)
                    px, py, pxx, pyy, _, _ = p[-1]
                    qx, qy, qxx, qyy, _, _ = q[-1]
                    assert py < qy and pyy < qyy
                    assert px >= qx and pxx >= qxx

                for j in range(len(q)):
                    if not is_row_bumped(q, j):
                        px, py, pxx, pyy, _, _ = p[j]
                        qx, qy, qxx, qyy, _, _ = q[j]
                        assert not is_row_bumped(p, j)
                        assert py == qy == pyy == qyy == j + 1
                        assert px < qx and pxx < qxx

                if mid(p) < len(p) and mid(q) < len(q):
                    assert mid(p) < mid(q)

                if not is_row_bumped(q, len(q) - 1):
                    assert not is_row_bumped(p, len(p) - 1)
                    px, py, pxx, pyy, _, _ = p[-1]
                    qx, qy, qxx, qyy, _, _ = q[-1]
                    assert py >= qy and pyy >= qyy
                    assert px < qx and pxx < qxx

            if a > b:
                for j in range(m):
                    if is_row_bumped(q, j):
                        px, py, pxx, pyy, _, _ = p[j]
                        qx, qy, qxx, qyy, _, _ = q[j]
                        assert is_row_bumped(p, j)
                        assert px == qx == pxx == qxx == j + 1
                        assert py >= qy and pyy >= qyy

                if not is_row_bumped(p, len(p) - 1):
                    assert not is_row_bumped(q, len(q) - 1)
                    px, py, pxx, pyy, _, _ = p[-1]
                    qx, qy, qxx, qyy, _, _ = q[-1]
                    assert py < qy and pyy < qyy
                    assert px >= qx and pxx >= qxx

                for j in range(len(p)):
                    if not is_row_bumped(p, j):
                        px, py, pxx, pyy, _, _ = p[j]
                        qx, qy, qxx, qyy, _, _ = q[j]
                        assert not is_row_bumped(q, j)
                        assert py == qy == pyy == qyy == j + 1
                        assert px >= qx and pxx >= qxx

                if is_row_bumped(q, len(q) - 1):
                    assert is_row_bumped(p, len(p) - 1)
                    px, py, pxx, pyy, _, _ = p[-1]
                    qx, qy, qxx, qyy, _, _ = q[-1]
                    assert py >= qy and pyy >= qyy
                    assert px < qx and pxx < qxx
        except:
            print('a =', a, 'b =', b)
            print(tab[i - 1])
            print(tab[i])
            print(tab[i + 1])
            print(p)
            print(q)
            assert False


def test_random_bumping_path(bound=10):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        _test_bumping_path(w)


def test_bumping_path(bound=7):
    for n in range(bound):
        for pi in Permutation.involutions(n):
            for w in pi.get_involution_words():
                _test_bumping_path(w)


def partial_insert(word, x):
    rows, _, _, _, _, pat = insert(word + (x,))
    weak_path = [(t[0], t[1]) for t in pat[-1]]
    strict_path = [(t[2], t[3]) for t in pat[-1]]

    j = [j for j in range(len(weak_path)) if weak_path[j][0] == weak_path[j][1]]
    if j:
        j = j[0] + 1
    else:
        j = len(weak_path)

    weak_row = weak_path[:j]
    weak_col = weak_path[j:]

    # j = [j for j in range(len(strict_path)) if strict_path[j][0] == strict_path[j][1]]
    # if j:
    #     j = j[0] + 1
    # else:
    #     j = len(strict_path)
    strict_row = strict_path[:j]
    strict_col = strict_path[j:]

    return rows, weak_row, weak_col, strict_row, strict_col


def help_test_complete_ddiag(a, ptab, qtab, q):
    try:
        assert abs(qtab.get(q - 1, q - 1)) + 1 == abs(qtab.get(q - 1, q)) == abs(qtab.get(q, q)) - 1
        if not qtab.get(q - 1, q).is_primed():
            i = qtab.get(q - 1, q).number - 1
            assert simplify_map(tau_permutation(a, i - 1)) == {}
            assert simplify_map(tau_permutation(a, i)) == {}
            assert simplify_map(tau_permutation(a, i + 1)) == {}

            x, y, z = a[i - 1:i + 2]
            if y < x < z or z < x < y:
                mid = (x, z, y)
            elif x == z:
                mid = (y, x, y)
            elif x < z < y or y < z < x:
                mid = (y, x, z)
            b = a[:i - 1] + mid + a[i + 2:]

            assert simplify_map(tau_permutation(b, i - 1)) == {}
            assert simplify_map(tau_permutation(b, i + 1)) == {}

            t1, weak_row1, weak_col1, strict_row1, strict_col1 = partial_insert(b[:i - 1], b[i - 1])
            t2, weak_row2, weak_col2, strict_row2, strict_col2 = partial_insert(b[:i], b[i])
            t3, weak_row3, weak_col3, strict_row3, strict_col3 = partial_insert(b[:i + 1], b[i + 1])

            assert weak_row1[-1] == (q - 1, q - 1) and len(weak_col1) == 0
            assert weak_row2[-1] == (q - 1, q - 1) and weak_col2[0] == (q - 1, q) and len(weak_col2) == 1
            assert weak_row3[-2] == (q - 1, q) and weak_row3[-1] == (q, q) and len(weak_col3) == 0

            g = full_cseq(a, i + 2)[0]
            assert full_cseq(b, i - 1)[0] == g[:-2]
            assert full_cseq(b, i)[0] == g[:-2] + (g[-1],)
            assert full_cseq(b, i + 1)[0] == g[:-1]
            assert full_cseq(b, i + 2)[0] == g

            assert simplify_map(tau_permutation(b, i)) == {g[-2]: g[-1], g[-1]: g[-2]}
            return [1]
    except:
        print(ptab)
        print(qtab)
        print('q =', q)
        print('a =', a)
        print('b =', b)
        assert False


def test_random_complete_ddiag(bound=15):
    count = None
    for n in range(bound):
        a = Permutation.random_involution_word(n)
        p, q = primed_cached_insert(a)
        for i in range(2, p.max_row + 1):
            if q.get(i - 1, i - 1).number + 2 == q.get(i, i).number:
                incr = help_test_complete_ddiag(a, p, q, i)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                if count and count[0] % 100 == 0:
                    print(count)


@pytest.mark.slow
def test_complete_ddiag(bound=7):
    count = None
    for pi in Permutation.involutions(bound):
        for a in pi.get_involution_words():
            p, q = primed_cached_insert(a)
            for i in range(2, p.max_row + 1):
                if q.get(i - 1, i - 1).number + 2 == q.get(i, i).number:
                    incr = help_test_complete_ddiag(a, p, q, i)
                    count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                    if count and count[0] % 100 == 0:
                        print(count)


def bump_differential(word, i):
    tab = insert(word[:i])[0]
    rows, weak_row, weak_col, strict_row, strict_col = partial_insert(word[:i], word[i])
    gamma = gamma_map(tab, word[i:])
    ans = []
    for j in range(len(weak_row)):
        if j == 0:
            d = word[i]
            eta = commutations(word).get(i, None)
        else:
            x = j
            y = ans[-1][0]
            z = ans[-1][1]
            d = tab[x - 1][z - x]
            if y == z and (x != y or ans[-1][-1] is not None):
                eta = gamma[(x, y)]
            else:
                eta = ans[-1][-1]
        ans += [(weak_row[j][1], strict_row[j][1], d, eta)]
    return ans


def help_test_bump_differential(word, seen):
    for i in range(len(word)):
        if (word[:i], word[i]) not in seen:
            seen.add((word[:i], word[i]))
            diff = bump_differential(word, i)
            cseq_0 = full_cseq(word, i)
            cseq_1 = full_cseq(word, i + 1)

            try:
                y, z, d, eta = diff[-1]
                p = len(diff)
                q = len(cseq_0[0])

                if p < y:
                    assert cseq_0 == cseq_1

                elif q < p:
                    assert p == y == z == q + 1
                    assert cseq_1 == (cseq_0[0] + (eta,), cseq_0[1] + (d,))

                elif p == y == z <= q:
                    c = cseq_0[1][p - 1]
                    theta = cseq_0[0][p - 1] if d + 1 == c else eta
                    cseq = [list(cseq_0[0][:]), list(cseq_0[1][:])]
                    cseq[0][p - 1] = theta
                    cseq[1][p - 1] = d
                    cseq = tuple(tuple(_) for _ in cseq)
                    assert cseq == cseq_1

                elif p == y < z:
                    assert p + 1 == z <= q
                    cseq = [list(cseq_0[0][:]), list(cseq_0[1][:])]
                    cseq[0][p - 1], cseq[0][p] = cseq[0][p], cseq[0][p - 1]
                    cseq = tuple(tuple(_) for _ in cseq)
                    assert cseq == cseq_1

                else:
                    assert False
            except:
                print(word[:i], word[i], word[i + 1:])
                print(diff)
                print()
                print('cseq_%s(a)' % str(i + 1))
                print(cseq_0[0])
                print(cseq_0[1])
                print()
                print('cseq_%s(a)' % str(i + 2))
                print(cseq_1[0])
                print(cseq_1[1])
                print()
                assert False


def test_random_bump_differential(bound=15):
    seen = set()
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_bump_differential(w, seen)


def test_bump_differential(bound=7):
    seen = set()
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            help_test_bump_differential(w, seen)


def compose(map1, map2):
    return {key: map1[map2[key]] for key in map2}


def simplify_map(m):
    return {k: v for (k, v) in m.items() if k != v}


def test_conservation_lemma(bound=6):
    reference = {}
    seen = {}
    for pi in Permutation.involutions(bound):
        for a in pi.get_involution_words():
            p, q = cached_insert(a)
            for i in range(len(a) - 2):
                key = (
                    full_cseq(a, i),
                    full_cseq(a, i + 2),
                    len({i + 1, i + 2} & entries_diagonal(q)),
                    len({i + 1, i + 2} & entries(q)),
                )
                val = simplify_map(compose(tau_permutation(a, i), tau_permutation(a, i + 1)))
                if a[i] < a[i + 1]:
                    if key in reference:
                        assert val == reference[key]
                    reference[key] = val
                    seen[key] = [(a, i)]
    print('*', len(reference), ' references')
    count = 0
    for pi in Permutation.involutions(bound):
        for a in pi.get_involution_words():
            p, q = cached_insert(a)
            for i in range(len(a) - 2):
                key = (
                    full_cseq(a, i),
                    full_cseq(a, i + 2),
                    len({i + 1, i + 2} & entries_diagonal(q)),
                    len({i + 1, i + 2} & entries(q)),
                )
                val = simplify_map(compose(tau_permutation(a, i), tau_permutation(a, i + 1)))
                try:
                    if key in reference:
                        assert val == reference[key]
                        count += 1
                except:
                    print('\na =', a)
                    print('i =', i)
                    print(p)
                    print(q)
                    print(key)
                    print(val)
                    print(reference[key])
                    print(seen[key])
                    assert False
    print('*', count, 'cross-referenced')


def test_tau_works(bound=5):
    for pi in Permutation.involutions(bound):
        for a in pi.get_primed_involution_words():
            tab = Word(*a).involution_insert()[0]
            gamma = gamma_map(tab)
            tau = tau_permutation(a)
            marked = {k for (k, v) in marked_cycles(a).items() if v != 1}
            # print('\na =', a)
            # print(tab)
            # print(gamma)
            # print(tau)
            # print(marked)
            for (i, j) in tab:
                entry = tab.get(i, j).number
                if entry < 0:
                    assert i != j
                    assert gamma[(i, j)] is not None
                    assert tau[gamma[(i, j)]] in marked
                if i != j and gamma[(i, j)] is not None and tau[gamma[(i, j)]] in marked:
                    assert entry < 0


def help_test_bumping_prop_two(a, i, j):
    try:
        assert i < j
        if a[i] < a[j] and all(a[t] < a[i] < a[j] for t in range(i + 1, j)):
            t1, weak_row1, weak_col1, strict_row1, strict_col1 = partial_insert(a[:i], a[i])
            t2, weak_row2, weak_col2, strict_row2, strict_col2 = partial_insert(a[:j], a[j])

            t1 = Tableau.shifted_from_rows(t1)
            t2 = Tableau.shifted_from_rows(t2)

            m = min(len(weak_row1), len(weak_row2))
            assert all(weak_row1[t] < weak_row2[t] for t in range(m))

            if j == i + 1 and weak_row2[-1][0] == weak_row2[-1][1]:
                assert len(weak_col1) > 0
                return [0, 1, 0, 0]

            return [1, 0, 0, 0]
        elif a[j] < a[i] and all(a[t] < a[i] < a[j] for t in range(i + 1, j)):
            t1, weak_row1, weak_col1, strict_row1, strict_col1 = partial_insert(a[:i], a[i])
            t2, weak_row2, weak_col2, strict_row2, strict_col2 = partial_insert(a[:j], a[j])

            t1 = Tableau.shifted_from_rows(t1)
            t2 = Tableau.shifted_from_rows(t2)

            m = min(len(weak_row1), len(weak_row2))
            assert all(weak_row1[t] >= weak_row2[t] for t in range(m))

            if weak_row1[-1][0] == weak_row1[-1][1]:
                assert len(weak_col2) > 0
                return [0, 0, 0, 1]

            return [0, 0, 1, 0]
    except:
        print('a =', a)
        print('i =', i)
        print('j =', j)
        print()
        for x in [t1, weak_row1, weak_col1, strict_row1, strict_col1]:
            print(x)
        print()
        for x in [t2, weak_row2, weak_col2, strict_row2, strict_col2]:
            print(x)
        print()
        assert False


def test_random_bumping_prop_two(bound=15):
    count = None
    for n in range(bound):
        a = Permutation.random_involution_word(n)
        for i in range(len(a) - 1):
            for j in range(i + 1, len(a)):
                incr = help_test_bumping_prop_two(a, i, j)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                if count and count[0] % 10 == 0:
                    print(count)


@pytest.mark.slow
def test_bumping_prop_two(bound=7):
    count = None
    pi = Permutation.longest_element(bound)
    for a in pi.get_involution_words():
        for i in range(len(a) - 1):
            for j in range(i + 1, len(a)):
                incr = help_test_bumping_prop_two(a, i, j)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                if count and count[0] % 10 == 0:
                    print(count)


def help_test_complete_acb(a, i):
    def reindexed_cseq(b, j):
        return full_cseq(b, j + 1)

    def find_in_row(path, row):
        return [(x, y) for (x, y) in path if x == row][0]

    def gam(b, j):
        return gamma_map(insert(b[:j + 1])[0], b[j + 1:])

    count = 13 * [0]
    j = None
    k = None
    try:
        a1, a2, a3 = a[i:i + 3]
        b = (a2, a1, a2) if a1 == a3 else (a2, a1, a3)
        b = a[:i] + b + a[i + 3:]

        t0 = insert(a[:i])[0]
        t3 = insert(a[:i + 3])[0]
        if len(t3) == len(t0) + 2:
            return
        tab = t0
        t0 = Tableau.shifted_from_rows(t0)
        t3 = Tableau.shifted_from_rows(t3)

        t1, weak_row1, weak_col1, strict_row1, strict_col1 = partial_insert(a[:i], a1)
        t2, weak_row2, weak_col2, strict_row2, strict_col2 = partial_insert(b[:i], a2)

        ta1 = Tableau.shifted_from_rows(t1)
        tb1 = Tableau.shifted_from_rows(t2)

        ta2, wweak_row1, wweak_col1, sstrict_row1, sstrict_col1 = partial_insert(a[:i + 1], a2)
        tb2, wweak_row2, wweak_col2, sstrict_row2, sstrict_col2 = partial_insert(b[:i + 1], a1)

        ta3, wwweak_row1, wwweak_col1, ssstrict_row1, ssstrict_col1 = partial_insert(a[:i + 2], a3)
        tb3, wwweak_row2, wwweak_col2, ssstrict_row2, ssstrict_col2 = partial_insert(b[:i + 2], a3 if a1 < a3 else a2)

        ta2 = Tableau.shifted_from_rows(ta2)
        tb2 = Tableau.shifted_from_rows(tb2)

        assert ta3 == tb3
        ta3 = Tableau.shifted_from_rows(ta3)
        tb3 = ta3

        if any(x != y for (x, y) in set(strict_row1) & set(strict_row2)):
            assert reindexed_cseq(a, i) == reindexed_cseq(b, i)
            assert tau_permutation(a, i) == tau_permutation(b, i)
            assert compose(tau_permutation(a, i + 1), tau_permutation(a, i + 2)) == compose(tau_permutation(b, i + 1), tau_permutation(b, i + 2))
            count[0] += 1
            return count
        if a1 < a3 < a2 and len(set(weak_row1) & set(weak_row2)) == 0:
            assert reindexed_cseq(a, i + 1) == reindexed_cseq(b, i + 1)
            assert tau_permutation(a, i + 2) == tau_permutation(b, i + 2)
            assert compose(tau_permutation(a, i), tau_permutation(a, i + 1)) == compose(tau_permutation(b, i), tau_permutation(b, i + 1))
            count[1] += 1
            return count

        if a1 == a3:
            j = 0
            u = a1
        else:
            j, jcol = sorted(set(weak_row1) & set(weak_row2))[0]
            assert (j, jcol) in t0
            u = t0.get(j, jcol).number

        cseq_start = reindexed_cseq(a, i - 1)
        cseq_final = reindexed_cseq(a, i + 2)
        g = lambda t: cseq_start[0][t - 1] # noqa
        e = lambda t: cseq_final[0][t - 1] # noqa

        if (j, j) not in weak_row1:
            k = weak_row1[-1][0]
            assert j < k
            exceptions = set()
            # (A1)
            for t in range(0 if j > 0 else 1, k - j):
                row = tab[j + t - 1]
                assert u + t in row and u + t + 1 in row

                (x, y) = find_in_row(strict_row1, j + t)
                assert t0.get(x, y).number == u + t

                (x, y) = find_in_row(weak_row2, j + t)
                assert t0.get(x, y).number == u + t

                (x, y) = find_in_row(strict_row2, j + t)
                assert t0.get(x, y).number == u + t + 1

                (x, y) = find_in_row(weak_row1, j + t)
                if u + t - 1 in row:
                    assert t0.get(x, y).number == u + t - 1
                else:
                    exceptions.add(j + t - 1)
                    assert t0.get(x, y).number == u + t
                    assert t0.get(x, y + 1).number == u + t + 1
                    assert (x, y + 1) in wweak_row1[:k - 1]

                if j > 0 == t:
                    assert u - 1 not in row
            # (A2)
            assert (k, k) in weak_row1
            # (A3)
            assert strict_row1[:k - 1] == sstrict_row2[:k - 1]
            assert strict_row2[:k - 1] == sstrict_row1[:k - 1]
            assert weak_row1[:k - 1] == wweak_row2[:k - 1]
            # (A4)
            assert all(strict_row2[t] == sstrict_row1[t] for t in range(k - 1) if t not in exceptions)
            # (A5)
            assert ssstrict_row1[:j] == ssstrict_row2[:j]
            assert wwweak_row1[:j] == wwweak_row2[:j]
            if j > 0:
                (x, y) = find_in_row(ssstrict_row1, j)
                assert t0.get(x, y).number == u + 1
                assert (x, y) in ssstrict_row2
                assert (x, y) in wwweak_row1
                assert (x, y) in wwweak_row2
            # (A6)
            for t in range(1, k - j):
                row = tab[j + t - 1]
                if u + t - 1 in row:
                    (x, y) = find_in_row(wwweak_row1, j + t)
                    assert t0.get(x, y).number == u + t - 1

                    (x, y) = find_in_row(ssstrict_row1, j + t)
                    assert t0.get(x, y).number == u + t

                    (x, y) = find_in_row(wwweak_row2, j + t)
                    assert t0.get(x, y).number == u + t

                    (x, y) = find_in_row(ssstrict_row2, j + t)
                    assert t0.get(x, y).number == u + t + 1
                else:
                    (x, y) = find_in_row(wwweak_row1, j + t)
                    assert t0.get(x, y).number == u + t

                    (x, y) = find_in_row(ssstrict_row1, j + t)
                    assert t0.get(x, y).number == u + t + 1

                    (x, y) = find_in_row(wwweak_row2, j + t)
                    assert t0.get(x, y).number == u + t + 1

                    (x, y) = find_in_row(ssstrict_row2, j + t)
                    assert t0.get(x, y).number == u + t + 1
            # (A7)
            v = u + k - j - 1
            if k > 1:
                (x, y) = find_in_row(strict_row1, k - 1)
                assert t0.get(x, y).number == v
                (x, y) = find_in_row(sstrict_row1, k - 1)
                assert ta1.get(x, y).number == v + 1
                (x, y) = find_in_row(ssstrict_row1, k - 1)
                assert ta2.get(x, y).number == v

                (x, y) = find_in_row(strict_row2, k - 1)
                assert t0.get(x, y).number == v + 1
                (x, y) = find_in_row(sstrict_row2, k - 1)
                assert tb1.get(x, y).number == v
                (x, y) = find_in_row(ssstrict_row2, k - 1)
                assert tb2.get(x, y).number == v + 1
            # (A8)
            assert t0.get(k, k).number in [v, v + 1, v + 2]
            count[2] += 1

            if t0.get(k, k).number == v:
                for t in [t0, ta1, ta2, ta3, tb1, tb2]:
                    assert t.get(k, k).number == v
                    assert t.get(k, k + 1).number == v + 1
                    assert t.get(k, k + 2).number == v + 2
                    assert t.get(k + 1, k + 1).number == v + 2
                    assert t.get(k + 1, k + 2).number == v + 3
                    assert t.get(k + 2, k + 2).number == v + 4

                assert simplify_map(tau_permutation(a, i)) == simplify_map(tau_permutation(b, i + 2)) == {g(k): g(k + 1), g(k + 1): g(k)}
                assert simplify_map(tau_permutation(a, i + 1)) == simplify_map(tau_permutation(b, i + 1)) == {g(k): g(k + 2), g(k + 2): g(k)}
                assert simplify_map(tau_permutation(a, i + 2)) == simplify_map(tau_permutation(b, i)) == {g(k + 1): g(k + 2), g(k + 2): g(k + 1)}
                count[3] += 1
            elif t0.get(k, k).number == v + 1:
                assert t0.get(k, k + 1).number == v + 2
                assert t0.get(k + 1, k + 1).number == v + 3
                if k > 1:
                    assert t0.get(k - 1, k + 1).number <= v + 1
                if k > 1 and t0.get(k - 1, k + 1).number == v + 1:
                    assert t0.get(k - 1, k).number == v

                    # assert ta1.get(k - 1, k).number == v - 1 #
                    assert ta1.get(k - 1, k + 1).number == v + 1
                    assert ta1.get(k, k).number == v
                    assert ta1.get(k, k + 1).number == v + 2
                    assert ta1.get(k + 1, k + 1).number == v + 3

                    # assert ta2.get(k - 1, k).number == v - 1 #
                    assert ta2.get(k - 1, k + 1).number == v
                    assert ta2.get(k, k).number == v
                    assert ta2.get(k, k + 1).number == v + 1
                    assert ta2.get(k + 1, k + 1).number == v + 2

                    # assert ta3.get(k - 1, k).number == v - 1 #
                    # assert ta3.get(k - 1, k + 1).number == v #
                    assert ta3.get(k, k).number == v
                    assert ta3.get(k, k + 1).number == v + 1
                    assert ta3.get(k + 1, k + 1).number == v + 2

                    assert tb1.get(k - 1, k).number == v
                    assert tb1.get(k - 1, k + 1).number == v + 1
                    assert tb1.get(k, k).number == v + 1
                    assert tb1.get(k, k + 1).number == v + 2
                    assert tb1.get(k + 1, k + 1).number == v + 3

                    # assert tb2.get(k - 1, k).number == v - 1 #
                    assert tb2.get(k - 1, k + 1).number == v + 1
                    assert tb2.get(k, k).number == v
                    assert tb2.get(k, k + 1).number == v + 2
                    assert tb2.get(k + 1, k + 1).number == v + 3
                    count[4] += 1
                else:
                    assert t0.get(k, k).number == v + 1
                    assert t0.get(k, k + 1).number == v + 2
                    assert t0.get(k + 1, k + 1).number == v + 3

                    assert ta1.get(k, k).number == v
                    assert ta1.get(k, k + 1).number == v + 1
                    assert ta1.get(k, k + 2).number == v + 2
                    assert ta1.get(k + 1, k + 1).number == v + 3

                    assert ta2.get(k, k).number == v
                    assert ta2.get(k, k + 1).number == v + 1
                    assert ta2.get(k, k + 2).number == v + 2
                    assert ta2.get(k + 1, k + 1).number == v + 2
                    assert ta2.get(k + 1, k + 2).number == v + 3

                    assert ta3.get(k, k).number == v
                    assert ta3.get(k, k + 1).number == v + 1
                    assert ta3.get(k, k + 2).number == v + 2
                    assert ta3.get(k + 1, k + 1).number == v + 2
                    assert ta3.get(k + 1, k + 2).number == v + 3

                    assert tb1.get(k, k).number == v + 1
                    assert tb1.get(k, k + 1).number == v + 2
                    assert tb1.get(k + 1, k + 1).number == v + 3

                    assert tb2.get(k, k).number == v
                    assert tb2.get(k, k + 1).number == v + 1
                    assert tb2.get(k, k + 2).number == v + 2
                    assert tb2.get(k + 1, k + 1).number == v + 3
                    count[5] += 1

                gg = gam(a, i - 1)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 1)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 2)
                assert gg[(k, k)] == g(k + 1)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(b, i)
                assert gg[(k, k)] == g(k + 1)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(b, i + 1)
                assert gg[(k, k)] == g(k + 1)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

            elif t0.get(k, k).number == v + 2:
                if k > 1:
                    assert t0.get(k - 1, k + 1).number < v + 2
                if k > 1 and ((k - 1, k + 2) not in t0 or t0.get(k - 1, k + 2).number >= v + 2):
                    assert t0.get(k - 1, k).number == v
                    assert t0.get(k - 1, k + 1).number == v + 1

                    # assert ta1.get(k - 1, k).number == v - 1 #
                    assert ta1.get(k - 1, k + 1).number == v + 1
                    assert ta1.get(k, k).number == v
                    assert ta1.get(k, k + 1).number == v + 2

                    # assert ta2.get(k - 1, k).number == v - 1 #
                    assert ta2.get(k - 1, k + 1).number == v
                    assert ta2.get(k, k).number == v
                    assert ta2.get(k, k + 1).number == v + 1
                    assert ta2.get(k + 1, k + 1).number == v + 2

                    # assert ta3.get(k - 1, k).number == v - 1 #
                    # assert ta3.get(k - 1, k + 1).number == v #
                    assert ta3.get(k, k).number == v
                    assert ta3.get(k, k + 1).number == v + 1
                    assert ta3.get(k + 1, k + 1).number == v + 2

                    assert tb1.get(k - 1, k).number == v
                    assert tb1.get(k - 1, k + 1).number == v + 1
                    assert tb1.get(k, k).number == v + 1
                    assert tb1.get(k, k + 1).number == v + 2

                    # assert tb2.get(k - 1, k).number == v - 1 #
                    assert tb2.get(k - 1, k + 1).number == v + 1
                    assert tb2.get(k, k).number == v
                    assert tb2.get(k, k + 1).number == v + 2
                    count[6] += 1
                else:
                    assert t0.get(k, k).number == v + 2

                    assert ta1.get(k, k).number == v
                    assert ta1.get(k, k + 1).number == v + 2

                    assert ta2.get(k, k).number == v
                    assert ta2.get(k, k + 1).number == v + 1
                    assert ta2.get(k + 1, k + 1).number == v + 2

                    assert ta3.get(k, k).number == v
                    assert ta3.get(k, k + 1).number == v + 1
                    assert ta3.get(k, k + 2).number == v + 2
                    assert ta3.get(k + 1, k + 1).number == v + 2

                    assert tb1.get(k, k).number == v + 1
                    assert tb1.get(k, k + 1).number == v + 2

                    assert tb2.get(k, k).number == v
                    assert tb2.get(k, k + 1).number == v + 1
                    assert tb2.get(k, k + 2).number == v + 2
                    count[7] += 1

                assert e(k) == g(k)

                gg = gam(a, i - 1)
                assert gg[(k, k)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i)
                assert gg[(k, k)] == e(k + 1)
                assert gg[(k, k + 1)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 1)
                assert gg[(k, k)] == e(k + 1)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(a, i + 2)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == e(k + 1)

                gg = gam(b, i)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(b, i + 1)
                assert gg[(k, k)] == g(k)
                if k > 1 and ((k - 1, k + 2) not in t0 or t0.get(k - 1, k + 2).number >= v + 2):
                    assert gg[(k, k + 1)] == e(k + 1)
                else:
                    assert gg[(k, k + 1)] is None
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)
            else:
                raise Exception
        else:
            assert a1 < a3 < a2
            k = j
            # (B1)
            assert all(strict_row1[t] == sstrict_row2[t] for t in range(k - 1))
            assert all(strict_row2[t] == sstrict_row1[t] for t in range(k - 1))
            assert all(strict_row1[t][0] != strict_row1[t][1] and strict_row1[t][1] < strict_row2[t][1] for t in range(k - 1))
            assert all(weak_row1[t] == wweak_row2[t] for t in range(k - 1))
            assert all(weak_row2[t] == wweak_row1[t] for t in range(k - 1))
            assert all(weak_row1[t][0] != weak_row1[t][1] and weak_row1[t][1] < weak_row2[t][1] for t in range(k - 1))
            # (B2)
            assert all(ssstrict_row1[t] == ssstrict_row2[t] for t in range(k - 1))
            assert all(ssstrict_row1[t][1] > strict_row1[t][1] and ssstrict_row1[t][1] <= strict_row2[t][1] for t in range(k - 1))
            assert all(wwweak_row1[t] == wwweak_row2[t] for t in range(k - 1))
            assert all(wwweak_row1[t][1] > weak_row1[t][1] and wwweak_row1[t][1] <= weak_row2[t][1] for t in range(k - 1))

            if k == 1:
                u = a1
                v = a3
                w = a2
            else:
                (x, y) = find_in_row(strict_row1, k - 1)
                u = t0.get(x, y).number

                (x, y) = find_in_row(ssstrict_row1, k - 1)
                v = ta2.get(x, y).number

                (x, y) = find_in_row(sstrict_row1, k - 1)
                w = ta1.get(x, y).number

            # (B3)
            assert u < v < w
            if k > 1:
                (x, y) = find_in_row(sstrict_row2, k - 1)
                assert u == tb1.get(x, y).number

                (x, y) = find_in_row(ssstrict_row2, k - 1)
                assert v == tb2.get(x, y).number

                (x, y) = find_in_row(strict_row2, k - 1)
                assert w == t0.get(x, y).number

            # (B4)
            assert (k, k) in t0
            assert t0.get(k, k).number >= w
            count[8] += 1

            if t0.get(k, k).number == w:
                if k > 1:
                    assert t0.get(k - 1, k + 1).number <= w
                if k > 1 and t0.get(k - 1, k + 1).number == w:
                    assert t0.get(k - 1, k).number == u
                    assert t0.get(k - 1, k + 1).number == w
                    assert t0.get(k, k).number == w
                    assert t0.get(k, k + 1).number == w + 1
                    assert t0.get(k + 1, k + 1).number == w + 2

                    # assert ta1.get(k - 1, k).number == ?
                    assert ta1.get(k - 1, k + 1).number == w
                    assert ta1.get(k, k).number == u
                    assert ta1.get(k, k + 1).number == w + 1
                    assert ta1.get(k + 1, k + 1).number == w + 2

                    # assert ta2.get(k - 1, k).number == ?
                    assert ta2.get(k - 1, k + 1).number == v
                    assert ta2.get(k, k).number == u
                    assert ta2.get(k, k + 1).number == w
                    assert ta2.get(k + 1, k + 1).number == w + 1

                    # assert ta3.get(k - 1, k).number == ?
                    # assert ta3.get(k - 1, k + 1).number == ?
                    assert ta3.get(k, k).number == u
                    assert ta3.get(k, k + 1).number == v
                    assert ta3.get(k + 1, k + 1).number == w

                    assert tb1.get(k - 1, k).number == u
                    assert tb1.get(k - 1, k + 1).number == v
                    assert tb1.get(k, k).number == w
                    assert tb1.get(k, k + 1).number == w + 1
                    assert tb1.get(k + 1, k + 1).number == w + 2

                    # assert tb2.get(k - 1, k).number == ?
                    assert tb2.get(k - 1, k + 1).number == v
                    assert tb2.get(k, k).number == u
                    assert tb2.get(k, k + 1).number == w
                    assert tb2.get(k + 1, k + 1).number == w + 2
                    count[9] += 1
                else:
                    assert t0.get(k, k).number == w
                    assert t0.get(k, k + 1).number == w + 1
                    assert t0.get(k + 1, k + 1).number == w + 2

                    assert ta1.get(k, k).number == u
                    assert ta1.get(k, k + 1).number == w
                    assert ta1.get(k, k + 2).number == w + 1
                    assert ta1.get(k + 1, k + 1).number == w + 2

                    assert ta2.get(k, k).number == u
                    assert ta2.get(k, k + 1).number == w
                    assert ta2.get(k, k + 2).number == w + 1
                    assert ta2.get(k + 1, k + 1).number == w + 1
                    assert ta2.get(k + 1, k + 2).number == w + 2

                    assert ta3.get(k, k).number == u
                    assert ta3.get(k, k + 1).number == v
                    assert ta3.get(k, k + 2).number == w + 1
                    assert ta3.get(k + 1, k + 1).number == w
                    assert ta3.get(k + 1, k + 2).number == w + 2

                    assert tb1.get(k, k).number == w
                    assert tb1.get(k, k + 1).number == w + 1
                    assert tb1.get(k, k + 2).number == w + 2
                    assert tb1.get(k + 1, k + 1).number == w + 2

                    assert tb2.get(k, k).number == u
                    assert tb2.get(k, k + 1).number == w
                    assert tb2.get(k, k + 2).number == w + 1
                    assert tb2.get(k + 1, k + 1).number == w + 2
                    count[10] += 1

                gg = gam(a, i - 1)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i)
                assert gg[(k, k)] == e(k)
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 2)
                assert gg[(k, k)] == e(k)
                assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(b, i)
                assert gg[(k, k)] == g(k + 1)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(b, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == g(k + 1)
                assert gg[(k + 1, k + 1)] == g(k)

            elif t0.get(k, k).number == w + 1:
                if k > 1:
                    assert t0.get(k - 1, k + 1).number <= w
                assert t0.get(k, k).number == w + 1

                assert ta1.get(k, k).number == u
                assert ta1.get(k, k + 1).number == w + 1

                assert ta2.get(k, k).number == u
                assert ta2.get(k, k + 1).number == w
                assert ta2.get(k + 1, k + 1).number == w + 1

                assert ta3.get(k, k).number == u
                assert ta3.get(k, k + 1).number == v
                assert ta3.get(k + 1, k + 1).number == w

                assert tb1.get(k, k).number == w
                assert tb1.get(k, k + 1).number == w + 1

                assert tb2.get(k, k).number == u
                assert tb2.get(k, k + 1).number == w

                gg = gam(a, i - 1)
                assert gg[(k, k)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] is None
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(a, i + 2)
                assert gg[(k, k)] == e(k)
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(b, i)
                assert gg[(k, k)] == g(k)
                assert gg[(k, k + 1)] is None
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(b, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                count[11] += 1
            elif t0.get(k, k).number > w + 1:
                x = t0.get(k, k).number
                if k > 1:
                    assert t0.get(k - 1, k + 1).number <= w
                assert t0.get(k, k).number == x

                assert ta1.get(k, k).number == u
                assert ta1.get(k, k + 1).number == x

                assert ta2.get(k, k).number == u
                assert ta2.get(k, k + 1).number == w
                assert ta2.get(k + 1, k + 1).number == x

                assert ta3.get(k, k).number == u
                assert ta3.get(k, k + 1).number == v
                assert ta3.get(k + 1, k + 1).number == w

                assert tb1.get(k, k).number == w
                assert tb1.get(k, k + 1).number == x

                assert tb2.get(k, k).number == u
                assert tb2.get(k, k + 1).number == w

                gg = gam(a, i - 1)
                assert gg[(k, k)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(a, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == e(k + 1)
                assert gg[(k + 1, k + 1)] == g(k)

                gg = gam(a, i + 2)
                assert gg[(k, k)] == e(k)
                assert gg[(k + 1, k + 1)] == e(k + 1)

                gg = gam(b, i)
                assert gg[(k, k)] == e(k + 1)
                assert gg[(k, k + 1)] == g(k)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                gg = gam(b, i + 1)
                assert gg[(k, k)] == e(k)
                assert gg[(k, k + 1)] == e(k + 1)
                if len(cseq_start[0]) > k:
                    assert gg[(k + 1, k + 1)] == g(k + 1)

                count[12] += 1
            else:
                raise Exception
        return count
    except:
        print('a =', a)
        print('i =', i)
        print()
        print('b =', b)
        print()
        print(a1, a2, a3)
        print()
        print(t0)
        print(t3)
        print(reindexed_cseq(a, i - 1))
        print()
        print(reindexed_cseq(a, i))
        print(reindexed_cseq(a, i + 1))
        print(reindexed_cseq(a, i + 2))
        print()
        print(reindexed_cseq(b, i))
        print(reindexed_cseq(b, i + 1))
        print(reindexed_cseq(b, i + 2))
        print()
        print(tau_permutation(a, i))
        print(tau_permutation(a, i + 1))
        print(tau_permutation(a, i + 2))
        print()
        print(tau_permutation(b, i))
        print(tau_permutation(b, i + 1))
        print(tau_permutation(b, i + 2))
        print()
        print(weak_row1, strict_row1)
        print(wweak_row1, sstrict_row1)
        print()
        print(weak_row2, strict_row2)
        print(wweak_row2, sstrict_row2)
        print()
        print('j =', j, 'k =', k)
        print('u =', u, 'v =', v)
        assert False


def test_simple_complete_acb():
    count = None

    def process(count, incr):
        return incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count

    a = (4, 13, 16, 15, 17, 6, 10, 8, 12, 5, 18, 20, 17, 3, 9, 4, 7, 14, 15, 1, 11, 16, 8, 13, 3, 2, 7, 12, 10, 11, 6, 9, 1, 3, 8, 19, 10, 13, 5, 18, 14, 11, 7, 8, 4, 12, 17, 20, 18, 15, 11, 16, 6, 19, 14, 13, 3, 2, 12, 5, 11, 14, 17, 9, 13, 6, 10, 7, 18, 15, 11, 4, 8, 6, 20, 17, 16, 9, 14, 19, 12, 13, 12, 10, 17, 18, 14, 15, 11, 20, 3, 12, 10, 13, 19, 14, 15, 13, 11, 16, 7, 15, 17, 8, 14, 9, 10, 16, 15, 1)
    i = 5
    count = process(count, help_test_complete_acb(a, i))

    a = (1, 5, 2, 3, 4, 2, 3, 2, 5)
    i = 5
    count = process(count, help_test_complete_acb(a, i))

    a = (5, 1, 3, 2, 4, 3, 1, 2, 1)
    i = 6
    count = process(count, help_test_complete_acb(a, i))

    a = (2, 4, 3, 1, 2, 1)
    i = 3
    count = process(count, help_test_complete_acb(a, i))

    a = (3, 1, 2, 1)
    i = 1
    count = process(count, help_test_complete_acb(a, i))

    a = (1, 3, 5, 4, 2, 3, 2, 4, 5)
    i = 4
    count = process(count, help_test_complete_acb(a, i))

    a = (5, 3, 6, 4, 5, 1, 2, 3, 2, 4, 5, 6)
    i = 6
    count = process(count, help_test_complete_acb(a, i))

    a = (6, 4, 1, 5, 2, 4, 3, 4, 2, 3, 5, 6)
    i = 4
    count = process(count, help_test_complete_acb(a, i))

    a = (6, 1, 2, 3, 4, 5, 2, 4, 3, 6, 4, 5)
    i = 6
    count = process(count, help_test_complete_acb(a, i))

    a = (5, 3, 4, 1, 3, 2, 3, 4, 5)
    i = 3
    count = process(count, help_test_complete_acb(a, i))

    a = (4, 1, 3, 2, 3, 4)
    i = 1
    count = process(count, help_test_complete_acb(a, i))

    a = (6, 1, 5, 3, 2, 3, 4, 5, 3, 4, 6, 5)
    i = 1
    count = process(count, help_test_complete_acb(a, i))

    a = (5, 1, 3, 2, 4, 5, 3, 4, 5)
    i = 1
    count = process(count, help_test_complete_acb(a, i))

    a = (7, 5, 6, 1, 3, 5, 2, 4, 3)
    i = 6
    count = process(count, help_test_complete_acb(a, i))

    a = (6, 5, 7, 6)
    i = 1
    count = process(count, help_test_complete_acb(a, i))

    a = (2, 5, 7, 4, 6, 5, 3, 4, 3)
    i = 6
    count = process(count, help_test_complete_acb(a, i))

    print()
    print(count)
    print()
    assert all(c > 0 for c in count)


def test_random_complete_acb(bound=15):
    count = None
    for n in range(bound):
        a = Permutation.random_involution_word(n)
        for i in range(len(a) - 2):
            a1, a2, a3 = a[i:i + 3]
            if a1 <= a3 < a2:
                incr = help_test_complete_acb(a, i)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                if count and count[0] % 10 == 0:
                    print(count)


@pytest.mark.slow
def test_complete_acb(bound=7):
    count = None
    for pi in Permutation.involutions(bound):
        for a in pi.get_involution_words():
            for i in range(len(a) - 2):
                a1, a2, a3 = a[i:i + 3]
                if a1 <= a3 < a2:
                    incr = help_test_complete_acb(a, i)
                    count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                    if count and count[0] % 10 == 0:
                        print(count)


def help_test_complete_disjoint(a, u, v):
    t0 = Tableau.shifted_from_rows(insert(a)[0])
    gamma = gamma_map(t0, (u, v))
    assert gamma == gamma_map(t0, (v, u))

    t1, weak_row1, weak_col1, strict_row1, strict_col1 = partial_insert(a, u)
    t2, weak_row2, weak_col2, strict_row2, strict_col2 = partial_insert(a, v)

    if set(weak_row1) & set(weak_row2):
        return

    n = len(a)
    auv = a + (u, v)
    avu = a + (v, u)
    assert bump_differential(auv, n) == bump_differential(avu, n + 1)

    uu = [u] + [t0.get(x, y).number for (x, y) in (strict_row1 + strict_col1)[:-1]]
    vv = [v] + [t0.get(x, y).number for (x, y) in (strict_row2 + strict_col2)[:-1]]

    t3, weak_row3, weak_col3, strict_row3, strict_col3 = partial_insert(a + (u,), v)

    def print_tab(tab, weak, strict):
        tab = Tableau.shifted_from_rows(tab) if type(tab) != Tableau else tab
        tab = Tableau({b: v for b, v in tab.mapping.items() if b in strict or b in weak})
        print(tab)

    cap = set(strict_row2) & set(strict_col1)
    wcap = set(strict_row2) & (set(weak_col1) - set(strict_col1))

    b1 = bump_differential(avu, n)
    b2 = bump_differential(auv, n + 1)

    try:
        assert len(cap) <= 1
        assert len(wcap) <= 1
        if len(wcap) == 1:
            assert len(cap) == 1
            assert list(cap)[0][0] == list(wcap)[0][0] + 1
            assert list(cap)[0][1] == list(wcap)[0][1]

        for i in range(min(len(weak_row1), len(weak_row2))):
            x1, y1 = weak_row1[i]
            x2, y2 = strict_row1[i]
            x3, y3 = weak_row2[i]
            x4, y4 = strict_row2[i]
            assert i + 1 == x1 == x2 == x3 == x4
            assert y1 <= y2 < y3 <= y4

        diag1 = [(x, y) for (x, y) in weak_row1 if x == y]
        diag2 = [(x, y) for (x, y) in weak_row2 if x == y]
        if diag1 and diag2:
            x1, y1 = diag1[0]
            x2, y2 = diag2[0]
            if (x1, y1) not in strict_row1:
                assert x1 + 1 < x2
            else:
                assert x1 < x2

        alt = set(weak_row2) & set(strict_col1)
        if alt:
            assert len(alt) == 1
            (i, j) = next(iter(alt))
            assert (i, j) in cap or (i, j + 1) in cap
        alt = set(strict_row2) & (set(weak_col1) - set(strict_col1))
        if alt:
            assert len(alt) == 1
            (i, j) = next(iter(alt))
            assert (i + 1, j) in cap
        alt = (set(weak_row2) - set(strict_row2)) & (set(weak_col1) - set(strict_col1))
        if alt:
            assert len(alt) == 1
            (i, j) = next(iter(alt))
            assert (i + 1, j + 1) in cap

        count = 10 * [0]
        if len(cap) == 0:
            assert b1 == b2
            count[0] += 1
        if len(cap) == 1:
            (i, j) = next(iter(cap))
            _j = j
            _x = t1[i - 1][j - i]
            _y = b1[i - 1][2]
            t1 = Tableau.shifted_from_rows(t1)

            j = min([y for (x, y) in strict_col1 if x == i]) - 1
            l = len([(x, y) for (x, y) in weak_col1 if x == i - 1 and (x + 1, y) in strict_col1])
            k = len([(x, y) for (x, y) in strict_col1 if x == i]) - l
            delta = _j - j

            if (i, j + k + l) in t0:
                assert i > 1

                assert all(t0.get(i, j + t).number == uu[j + t] for t in range(1, k + l + 1))
                assert all(t0.get(i, j + k + t).number == uu[j + k] + t for t in range(1, l + 1))
                assert all(t0.get(i - 1, j + k + t).number == uu[j + k] + t - 1 for t in range(1, l + 1))

                assert all(t1.get(i, j + t).number == uu[j + t - 1] for t in range(1, k + 1))
                assert all(t1.get(i, j + k + t).number == uu[j + k] + t for t in range(1, l + 1))
                assert all(t1.get(i - 1, j + k + t).number == uu[j + k] + t - 1 for t in range(1, l + 1))

                if (i - 1, j + k + l + 1) in t0:
                    assert t0.get(i - 1, j + k + l + 1).number > uu[j + k] + l
                if (i, j + k + l + 1) in t0:
                    assert t0.get(i, j + k + l + 1).number > uu[j + k] + l + 1

                assert (i - 1, j + k + l) not in set(weak_row2) - set(strict_row2)
                assert (i, j + k + l) not in set(weak_row2) - set(strict_row2)

                assert i <= j + 1
                assert i < j + 1 or (k == 0 and t0.get(i - 1, i - 1).number + 1 == t0.get(i - 1, i).number == t0.get(i, i).number - 1 == uu[j])
                assert i != j or t0.get(j, j).number == uu[j]
                if i < j:
                    col = []
                    x = i + 1
                    while (x, j) in t0:
                        col.append(t0.get(x, j).number)
                        x += 1
                    assert uu[j] in col

                if k < delta <= l:
                    count[1] += 1
                    assert (i - 1, j + delta) in strict_row2
                if k + 1 < delta <= l:
                    count[2] += 1
                    assert b1 == b2
                if k > 0 and delta == k + 1:
                    count[3] += 1
                    assert l > 0
                    assert (i - 1, j + k + 1) in strict_row2
                    assert (i, j + k + 1) in strict_row2
                    assert t0.get(i - 1, j + k).number != uu[j + k] - 1
                    assert (i - 1, j + k + 1) in weak_row2

                    tup1, tup2 = b1[i - 1], b1[i]
                    theta = gamma[(i - 1, j + k + 1)]
                    eta = gamma[(i, j + k + 1)]
                    y, yy = tup2[:2]

                    assert tup1 == (j + k, j + k + 1, uu[j + k], theta)
                    assert tup2 == (y, yy, uu[j + k] + 1, theta)

                    tup1, tup2 = b2[i - 1], b2[i]
                    assert tup1 == (j + k + 1, j + k + 1, uu[j + k], eta)
                    assert tup2 == (y, yy, uu[j + k] + 1, theta)

                    assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1, i])
                if k == 0 and delta == 1:
                    count[4] += 1
                    assert l > 0
                    assert t0.get(i, j).number <= uu[j] - 1
                    assert t0.get(i + 1, j).number <= uu[j]

                    theta = gamma[(i - 1, j + 1)]
                    eta = gamma[(i, j + 1)]
                    tup1, tup2 = b1[i - 1], b1[i]
                    assert tup1 == (j + 1, j + 1, uu[j], theta)
                    assert tup2 == (j + 1, j + 1, uu[j] + 1, eta)

                    tup1, tup2 = b2[i - 1], b2[i]
                    assert tup1 == (j + 1, j + 1, uu[j], eta)
                    assert tup2 == (j + 1, j + 1, uu[j] + 1, theta)

                    assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1, i])
                    assert len(b1) > i + 1 or i < j

                assert not (k > 1 and delta == k)

                if (k > 0 and 1 < delta < k) or (k > 0 and delta == 1 and uu[j] < vv[i - 1]):
                    count[5] += 1
                    assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1])
                    y, yy, d, eta = b1[i - 1]
                    assert b2[i - 1] == (1 + y, 1 + yy, d, eta)
                    assert len(b1) > i

                if k > 0 and delta == 1 and vv[i - 1] <= uu[j]:
                    count[6] += 1
                    assert t0.get(i, j).number < vv[i - 1]

                    assert t0.get(i + 1, j).number <= uu[j]
                    col = []
                    x = i + 1
                    while (x, j) in t0:
                        col.append(t0.get(x, j).number)
                        x += 1
                    assert uu[j] in col

                    assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1, i])
                    assert len(b1) > i + 1 or i < j

                    tup1, tup2 = b1[i - 1], b1[i]
                    theta, eta = tup1[-1], tup2[-1]
                    assert tup1 == (j + 1, j + 1, vv[i - 1], theta)
                    assert tup2 == (j + 1, j + 1, uu[j + 1], eta)

                    tup1, tup2 = b2[i - 1], b2[i]
                    phi = tup2[-1]
                    if vv[i - 1] == uu[j]:
                        assert tup1 == (j + 1, j + 2, vv[i - 1], theta)
                        assert tup2 == (j + 1, j + 1, uu[j + 1], theta)
                    else:
                        assert tup1 == (j + 1, j + 1, vv[i - 1], theta)
                        assert tup2 == (j + 1, j + 1, uu[j], phi)

            else:
                assert l == 0
                assert all(t0.get(i, j + t).number == uu[j + t] for t in range(1, k))
                assert all(t1.get(i, j + t).number == uu[j + t - 1] for t in range(1, k + 1))
                if uu[j] < vv[i - 1]:
                    count[7] += 1
                    y, yy, d, eta = b1[i - 1]
                    assert b2[i - 1] == (1 + y, 1 + yy, d, eta)
                    assert len(b1) > i or y == yy == j + k
                    assert k > 0
                    assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1])
                else:
                    assert t0.get(i, j).number < vv[i - 1]

                    assert t0.get(i + 1, j).number <= uu[j]
                    col = []
                    x = i + 1
                    while (x, j) in t0:
                        col.append(t0.get(x, j).number)
                        x += 1
                    assert uu[j] in col

                    if delta < k:
                        count[8] += 1
                        assert len(b1) == len(b2) and all(b1[t] == b2[t] for t in range(len(b1)) if t not in [i - 1, i])
                        assert len(b1) > i + 1 or i < j

                        tup1, tup2 = b1[i - 1], b1[i]
                        theta, eta = tup1[-1], tup2[-1]
                        assert tup1 == (j + 1, j + 1, vv[i - 1], theta)
                        assert tup2 == (j + 1, j + 1, uu[j + 1], eta)

                        tup1, tup2 = b2[i - 1], b2[i]
                        phi = tup2[-1]
                        if vv[i - 1] == uu[j]:
                            assert tup1 == (j + 1, j + 2, vv[i - 1], theta)
                            assert tup2 == (j + 1, j + 1, uu[j + 1], theta)
                        else:
                            assert tup1 == (j + 1, j + 1, vv[i - 1], theta)
                            assert tup2 == (j + 1, j + 1, uu[j], phi)
                    else:
                        count[9] += 1
                        assert k == 1
                        assert all(b1[t] == b2[t] for t in range(len(b1) - 1))
                        assert len(b1) == i
                        assert len(b2) == i + 1

                        tup1 = b1[i - 1]
                        theta = tup1[-1]
                        assert tup1 == (j + 1, j + 1, vv[i - 1], theta)

                        tup1, tup2 = b2[i - 1], b2[i]
                        phi = tup2[-1]
                        assert vv[i - 1] < uu[j]
                        assert tup1 == (j + 1, j + 1, vv[i - 1], theta)
                        assert tup2 == (j + 1, j + 1, uu[j], phi)

            # a0
            #    b1 b2 ... bk bk+1 bk+2 ... bk+l
            #    a1 a2 ... ak bk   bk+1 ... bk+l-1 c
            #
        #     j = _j
        #     x = _x
        #     y = _y
        #     if (i, j) in t0:
        #         assert weak_row2[i:] == weak_row3[i:]
        #         assert strict_row2[i:] == strict_row3[i:]
        #         assert len(b1) == len(b2) > i

        #     # intersect at b1
        #     if (i, j - 1) not in strict_col1 and (i, j) in weak_col1 and (i, j) in t0:
        #         if x < y:
        #             c1, c2, d, eta = b1[i - 1]
        #             assert (c1 + 1, c2 + 1, d, eta) == b2[i - 1]
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 != i)
        #         elif x == y:
        #             c1, c2, d, eta = b1[i - 1]
        #             assert (c1, c2 + 1, d, eta) == b2[i - 1]
        #             c1, c2, d, _ = b1[i]
        #             assert (c1, c2, d, eta) == b2[i]
        #             assert len(b1) == len(b2) > i + 1 or i + 1 < c1
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 not in [i, i + 1])
        #         else:
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 != i + 1)
        #             c1, c2, d, eta = b1[i]
        #             assert b2[i][:2] == (c1, c2)
        #             assert len(b1) == len(b2) > i + 1 or i + 1 < c1

        #         assert (i, j) in strict_row3 or (i, j + 1) in strict_row3
        #         assert (i + 1, j) in strict_row3
        #         assert (i + 1, j) in strict_row2
        #         assert (i + 1, j) in weak_row2
        #         count[0] += 1
        #     elif (i, j) not in t0:
        #         assert x != y
        #         if x < y:
        #             c1, c2, d, eta = b1[i - 1]
        #             assert c1 == c2 and (c1 + 1, c2 + 1, d, eta) == b2[i - 1]
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 != i)
        #         else:
        #             assert b1 == b2[:-1] and i == len(b1) == len(b2) - 1
        #             c1, c2, d, eta = b2[-1]
        #             assert d == x and i + 1 < c1 == c2

        #         assert len(strict_col2) == len(strict_col3) == 0
        #         assert weak_row3[-1] in [(i, j + 1), (i + 1, j)]
        #         assert len(set(weak_row3[-1])) != 1
        #         assert weak_row2[-1] == (i, j)
        #         count[1] += 1

        #     # intersect at bk+?
        #     elif (i, j) not in weak_col1:
        #         if (i - 1, j - 1) in weak_col1:
        #             assert b1 == b2
        #         elif (i, j - 1) in strict_col1:
        #             assert i > 1
        #             c1, c2, _, _ = b1[i - 2]
        #             assert c1 == c2 == j
        #             c1, c2, d, _ = b1[i - 1]
        #             assert c1 + 1 == c2 and b2[i - 1][:-1] == (c1 + 1, c2, d)
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 != i)

        # #  a0    [a0+2]
        # # [a0-1]  a0+1  a0+2 ... a0+l
        # # [a0-2]  a0    a0+1 ... a0+l-1 c
        # #
        #         # intersect at a0+1
        #         else:
        #             assert 1 < i < len(b1) == len(b2)
        #             c1, c2, d, _ = b1[i - 1]
        #             assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 not in [i, i + 1])
        #             assert c1 == c2 and (c1, c2, d + 1) == b1[i][:-1] == b2[i][:-1]

        #         count[2] += 1
        #     elif (i, j - 1) in weak_col1 and (i, j) in weak_col1:
        #         c1, c2, d, eta = b1[i - 1]
        #         assert (c1 + 1, c2 + 1, d, eta) == b2[i - 1]
        #         assert all(b1[k] == b2[k] for k in range(len(b1)) if k + 1 != i)

        #         assert (i, j + 1) in strict_row3
        #         count[3] += 1
        #     else:
        #         assert False

    except:
        print(a, u, v, '\n')
        print(b1)
        print(b2)
        print()
        print('x =', _x, 'y =', _y)
        print()
        print_tab(t0, [(p, q) for (p, q) in strict_col1 if p == i], [(p, q) for (p, q) in weak_col1 if p == i - 1 and (p + 1, q) in strict_col1])
        print_tab(t1, [(p, q) for (p, q) in strict_col1 if p == i], [(p, q) for (p, q) in weak_col1 if p == i - 1 and (p + 1, q) in strict_col1])
        print()
        print('i =', i, 'j =', j, 'k =', k, 'l =', l, 'delta =', delta)
        print('\ntab =\n', t0)
        print('a =', a)
        print('u =', u, uu)
        print('v =', v, vv)
        print(strict_col1, strict_row2, strict_row3)
        assert False
    return count


def test_complete_disjoint_simple():
    a = (8, 10, 9, 5, 4, 13, 14, 6, 2, 10, 7, 6, 12, 11, 3, 4, 5, 8, 12)
    u = 4
    v = 6
    help_test_complete_disjoint(a, u, v)

    a = (5, 8, 7, 14, 11, 10, 3, 9, 4, 1, 12, 13, 16, 8, 6, 14, 2, 11, 15, 5, 3, 7, 6, 16, 10, 8, 12, 5, 14, 9, 2, 13, 8, 7)
    u = 6
    v = 8
    help_test_complete_disjoint(a, u, v)

    a = (1, 14, 12, 11, 13, 8, 10, 4, 7, 2, 14, 11, 6, 3, 9, 12, 7, 8, 5)
    u = 2
    v = 4
    help_test_complete_disjoint(a, u, v)

    a = (5, 1, 4, 10, 3, 12, 7, 11, 9, 6, 10, 2)
    u = 5
    v = 7
    help_test_complete_disjoint(a, u, v)

    a = (3, 9, 7)
    u = 4
    v = 8
    help_test_complete_disjoint(a, u, v)

    a = (3, 5, 6, 1, 4, 5, 6, 2, 3)
    u = 1
    v = 4
    help_test_complete_disjoint(a, u, v)


def test_random_complete_disjoint(bound=15):
    count = None
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        for i in range(len(w) - 1):
            a, u, v = w[:i], w[i], w[i + 1]
            if u + 1 < v:
                incr = help_test_complete_disjoint(a, u, v)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                print(count)


@pytest.mark.slow
def test_complete_disjoint(bound=7):
    wseen = set()
    seen = set()
    count = None
    pi = Permutation.longest_element(bound)
    for w in pi.get_involution_words():
        for i in range(len(w) - 1):
            a, u, v = w[:i], w[i], w[i + 1]
            if u + 1 >= v:
                continue
            if (a, u, v) in wseen:
                continue
            else:
                wseen.add((a, u, v))
            tab = Tableau.from_rows(insert(a)[0])
            if u + 1 < v and (tab, u, v) not in seen:
                seen.add((tab, u, v))
                incr = help_test_complete_disjoint(a, u, v)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                print(count)
    assert all(c > 0 for c in count)


def help_test_disjoint(a, u, v):
    try:
        t1, w1, cw1, s1, cs1 = partial_insert(a, u)
        t2, w2, cw2, s2, cs2 = partial_insert(a, v)
    except:
        print(a, u, v)
        assert False
    intersect = set(w1) & set(w2)
    strict = set(s1) & set(s2)

    count = [0, 0, 0]

    def printout(a, u, v, t1, t2, w1, s1, w2, s2, intersect, strict):
        print('a =', a)
        print('u =', u)
        print('v =', v)
        t0 = insert(a)[0]
        print(Tableau.shifted_from_rows(t0))
        print(Tableau.shifted_from_rows(t1))
        print(Tableau.shifted_from_rows(t2))
        print(Tableau.shifted_from_rows(insert(a + (u, v))[0]))
        print(Tableau.shifted_from_rows(insert(a + (v, u))[0]))
        print()
        print(w1, s1)
        print()
        print(w2, s2)
        print()
        print(intersect, ':', strict)
        print()
        print(cseq(t0, (u, v)))
        print(cseq(t1, (v,)))
        print()
        print(cseq(t0, (u, v)))
        print(cseq(t2, (u,)))
        print('\n\n\n')

    if any(x != y for (x, y) in strict):
        count[0] += 1

        q1 = involution_insert(*tuple(Word(_) for _ in a + (u,)))[1]
        q2 = involution_insert(*tuple(Word(_) for _ in a + (v,)))[1]
        n = len(a)
        e1 = entries(q1)
        e2 = entries(q2)
        diag1 = {q1[box].number for box in q1 if box[0] == box[1]}
        diag2 = {q2[box].number for box in q2 if box[0] == box[1]}

        (x, y) = sorted(strict)[0]

        try:
            if u + 1 == v:
                assert cseq(t1, (v, u)) == cseq(t2, (u, v))
            else:
                assert cseq(t1, (v,)) == cseq(t2, (u,))

            if -n - 1 in e1:
                assert -n - 1 in e2

            if n + 1 in diag1:
                assert n + 1 in diag2

            assert cs1 == cs2
            assert cw1 == cw2
            assert s1[x - 1:] == s2[x - 1:]
            assert w1[x - 1:] == w1[x - 1:]
            assert (x, y) in intersect
        except:
            printout(a, u, v, t1, t2, w1, s1, w2, s2, intersect, strict)
            # print('(1) cseq:', cseq(t1, (v,)), '==', cseq(t2, (u,)))
            assert False

    if u + 1 == v:
        return count

    if len(intersect) == 0:
        count[1] += 1
        tt1 = insert(a + (u, v))[0]
        tt2 = insert(a + (v, u))[0]

        q1 = involution_insert(*tuple(Word(_) for _ in a + (u, v)))[1]
        q2 = involution_insert(*tuple(Word(_) for _ in a + (v, u)))[1]
        n = len(a)
        e1 = {_.number for _ in q1.entries()}
        e2 = {_.number for _ in q2.entries()}
        diag1 = {q1[box].number for box in q1 if box[0] == box[1]}
        diag2 = {q2[box].number for box in q2 if box[0] == box[1]}

        try:
            assert len(strict) == 0
            assert len(set(s1) & set(w2)) == 0
            assert cseq(tt1) == cseq(tt2)

            if -n - 1 in e1:
                assert -n - 2 in e2
            if -n - 2 in e1:
                assert -n - 1 in e2

            if n + 1 in diag1:
                assert n + 2 in diag2
            if n + 2 in diag1:
                assert n + 1 in diag2

        except:
            printout(a, u, v, t1, t2, w1, s1, w2, s2, intersect, strict)
            print('(2) cseq:', cseq(tt1), '==', cseq(tt2))
            print(q1)
            print(q2)
            assert False

    if not any(x != y for (x, y) in strict):
        t3, w3, _, s3, _ = partial_insert(a + (v,), u)
        t4, w4, _, s4, _ = partial_insert(a + (u,), v)
        if s1 != s3 or w1 != w3:
            printout(a, u, v, t1, t2, w1, s1, w2, s2, intersect, strict)
            print(w3, s3, '!=', s1)
            print(w4, s4, '!=', s2)
            assert False

    if not any(x != y for (x, y) in strict) and any(x != y for (x, y) in intersect):
        try:
            x = min([x for (x, y) in intersect if x != y])
            c = w1[x - 1][1] - x + 1
            y = w1[-1][0]
            if x < y:
                print(x, c, y)
                count[2] += 1
                rows = insert(a)[0]
                val = rows[x - 1][c - 1]
                for i in range(y - x):
                    assert val + i in rows[x - 1 + i]
                    assert val + i + 1 in rows[x - 1 + i]
                assert (val + y - x + 2) in rows[y - 1] or (val + y - x + 1) in rows[y - 1]
        except:
            printout(a, u, v, t1, t2, w1, s1, w2, s2, intersect, strict)
            assert False

    return count


def test_disjoint_simple():
    a = (6, 4, 1, 5)
    u = 2
    v = 4
    help_test_disjoint(a, u, v)

    a = (3, 5, 4)
    u = 1
    v = 3
    help_test_disjoint(a, u, v)

    a = (4, 1, 6, 5)
    u = 2
    v = 4
    help_test_disjoint(a, u, v)

    a = (11, 8, 10, 3, 15, 7, 19, 20, 9, 4, 14, 17, 12, 2, 6, 11, 5, 1, 18, 13, 8, 4, 14, 3, 10, 22, 7, 8, 16, 2, 19, 12, 15, 11, 13, 17, 6, 14, 18, 1, 12, 16, 15, 20, 9, 17, 19, 7, 8, 9, 4, 18, 5, 10, 4, 9, 13, 14, 6, 19, 21, 22, 20, 11, 3, 16, 7, 6, 8, 12, 21, 19, 11, 7, 20, 18, 13)
    u = 15
    v = 17
    help_test_disjoint(a, u, v)

    a = (5, 1, 13, 3, 19, 15, 7, 11, 17, 21, 9, 10, 11, 4, 14, 12, 2, 8, 13, 11, 6, 16, 18, 9, 14, 15, 20, 5, 21, 3, 19, 14, 17, 18, 4, 16, 7, 17, 6, 12, 5, 13, 11, 10, 11, 14, 3, 2, 7, 18, 8, 12, 9, 20, 15, 19, 21, 11, 13, 6, 4)
    u = 7
    v = 10
    help_test_disjoint(a, u, v)


def test_random_disjoint(bound=15):
    c1, c2, c3 = 0, 0, 0
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        for i in range(len(w) - 1):
            a = w[:i]
            u = w[i]
            v = w[i + 1]
            if u + 1 < v or (i + 2 < len(w) and u + 1 == v == w[i + 2] + 1):
                count = help_test_disjoint(a, u, v)
                c1 += count[0]
                c2 += count[1]
                c3 += count[2]
                print(c1, c2, c3)


@pytest.mark.slow
def test_disjoint(bound=7):
    wseen = set()
    seen = set()
    c1, c2, c3 = 0, 0, 0
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            for i in range(len(w) - 1):
                a = w[:i]
                u = w[i]
                v = w[i + 1]

                check = u + 1 < v or (i + 2 < len(w) and u + 1 == v == w[i + 2] + 1)

                if not check:
                    continue
                if (a, u, v) in wseen:
                    continue
                else:
                    wseen.add((a, u, v))
                t = Tableau.from_rows(insert(a)[0])
                if (t, u, v) in seen:
                    continue
                if check:
                    seen.add((t, u, v))
                    count = help_test_disjoint(a, u, v)
                    c1 += count[0]
                    c2 += count[1]
                    c3 += count[2]
                    print(c1, c2, c3)


def gamma_map(tab, b=()):
    rows = tab.get_rows() if type(tab) == Tableau else tab
    word = ()
    for row in reversed(rows):
        word += tuple(row)
    word += tuple(b)
    comm = commutations(word)
    ans = {}
    for i in range(len(rows)):
        for j in range(len(rows[i])):
            index = j + sum([len(r) for r in rows[i + 1:]])
            ans[(i + 1, i + j + 1)] = comm.get(index, None)
    return ans


def full_cseq(word, j):
    if (word, j) not in full_cseq_cache:
        rows = insert(word[:j])[0]
        b = word[j:]
        gamma = gamma_map(rows, b)
        full_cseq_cache[(word, j)] = (
            tuple(gamma[(i, i)] for i in range(1, 1 + len(rows))),
            tuple(rows[i - 1][0] for i in range(1, 1 + len(rows))),
        )
    return full_cseq_cache[(word, j)]


def tau_permutation(word, j=None):
    if (word, j) not in tau_cache:
        if j is None and len(word) == 0:
            ans = {}
        elif j is None and len(word) > 0:
            ans = tau_permutation(word, 0)
            for k in range(1, len(word)):
                ans = compose(ans, tau_permutation(word, k))
        else:
            comm = commutations(word)
            ans = {p: p for p in comm.values()}
            a = full_cseq(word, j)[0]
            b = full_cseq(word, j + 1)[0]
            rng = [t for t in range(len(a)) if a[t] != b[t]]
            pq = {a[j] for j in rng} | {b[j] for j in rng}
            if pq:
                assert len(pq) == 2
                p, q = tuple(pq)
                ans[p], ans[q] = ans[q], ans[p]
        tau_cache[(word, j)] = ans
    return tau_cache[(word, j)]


def cseq(tab, b=()):
    rows = tab.get_rows() if type(tab) == Tableau else tab
    gamma = gamma_map(rows, b)
    return [gamma[(i, i)] for i in range(1, 1 + len(rows))]


def help_test_gamma(w):
    _, _, _, _, tabs, paths = insert(w)
    tabs = [Tableau()] + tabs
    for i in range(len(w)):
        g1 = gamma_map(tabs[i], w[i:])
        g2 = gamma_map(tabs[i + 1], w[i + 1:])
        path = paths[i]

        # print()
        # print(w, i)
        # print(tabs[i])
        # print(tabs[i + 1])
        # print(' first gamma map:', g1)
        # print('second gamma map:', g2)
        # print('path =', path)

        row = tabs[i].row_reading_word()
        gamma = [commutations(row + w[i:]).get(len(row), None)]
        ww = [w[i]] + [tabs[i].get(x2, y2).number for (_, _, x2, y2, _, _) in path[:-1]]
        for j in range(len(path) - 1):
            (x1, y1, x2, y2, _, _) = path[j]
            if x1 == y1:
                assert gamma[-1] is None or ww[j] < ww[j + 1]
            if (x1, y1) == (x2, y2) and (x1 != y1 or gamma[-1] is not None):
                gamma.append(g1[(x1, y1)])
            else:
                gamma.append(gamma[-1])

        # print('gamma list =', gamma)
        seen = set()
        index = 0
        for (x1, y1, x2, y2, _, _) in path:
            # print('* index:', index, x1 == y1, gamma[index])
            if (x1, y1) == (x2, y2):
                if x1 != y1 or gamma[index] is not None:
                    assert g2[(x1, y1)] == gamma[index]
                else:
                    assert g2[(x1, y1)] == g1[(x1, y1)]
                seen.add((x1, y1))
            elif x1 != y1 and x2 != y2:
                assert g2[(x1, y1)] == g1[(x2, y2)]
                assert g1[(x1, y1)] == g2[(x2, y2)]
                seen.add((x1, y1))
                seen.add((x2, y2))
            elif x1 == y1:
                x = x1
                assert g2[(x, x)] == g1[(x + 1, x + 1)] and g2[(x, x)] is not None
                assert g2[(x + 1, x + 1)] == g1[(x, x)] and g2[(x + 1, x + 1)] is not None
                assert g1[(x, x + 1)] is None and g2[(x, x + 1)] is None
                seen.add((x, x))
                seen.add((x, x + 1))
                seen.add((x + 1, x + 1))
            index += 1
        for box in g2:
            if box not in seen:
                assert g2[box] == g1[box]


def test_random_gamma(bound=15):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_gamma(w)


@pytest.mark.slow
def test_gamma(bound=7):
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            help_test_gamma(w)


def help_test_bac(w):
    for i in range(len(w) - 2):
        b, a, c = w[i: i + 3]
        if a < b < c:
            u = w[:i + 3]
            v = w[:i] + (b, c, a)
            urows, _, utau, _, _, _ = insert(u)
            vrows, _, vtau, _, _, _ = insert(v)
            assert urows == vrows
            assert utau == vtau


def test_random_bac(bound=15):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_bac(w)


@pytest.mark.slow
def test_bac(bound=7):
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            help_test_bac(w)


def help_test_acb(w):
    ans = [0, 0]
    for i in range(len(w) - 2):
        a, c, b = w[i: i + 3]
        if a <= b < c:
            ans[int(a == b)] += 1
            u = w[:i + 3]
            v = w[:i] + ((c, a, b) if a < b else (c, a, c))
            trows, _, ttau, _, _, _ = insert(w[:i])
            urows, _, utau, _, _, _ = insert(u)
            vrows, _, vtau, _, _, _ = insert(v)
            if len(urows) != len(trows) + 2 and utau != vtau:
                print(u, '<->', v)
                print()
                print(ttau)
                print(Tableau.shifted_from_rows(trows), cseq(trows, w[i: i + 3]))
                # print(Tableau(gamma_map(trows, w[i: i + 3])))
                trows, _, _, _, _, _ = insert(w[:i + 1])
                print(Tableau.shifted_from_rows(trows), cseq(trows, w[i + 1: i + 3]))
                # print(Tableau(gamma_map(trows, w[i + 1: i + 3])))
                trows, _, _, _, _, _ = insert(w[:i + 2])
                print(Tableau.shifted_from_rows(trows), cseq(trows, w[i + 2: i + 3]))
                # print(Tableau(gamma_map(trows, w[i + 2: i + 3])))
                trows, _, _, _, _, _ = insert(w[:i + 3])
                print(Tableau.shifted_from_rows(trows), cseq(trows))
                # print(Tableau(gamma_map(trows)))
                print()
                print('***')
                print()
                trows, _, _, _, _, _ = insert(v[:i])
                print(Tableau.shifted_from_rows(trows), cseq(trows, v[i: i + 3]))
                # print(Tableau(gamma_map(trows, v[i: i + 3])))
                trows, _, _, _, _, _ = insert(v[:i + 1])
                print(Tableau.shifted_from_rows(trows), cseq(trows, v[i + 1: i + 3]))
                # print(Tableau(gamma_map(trows, v[i + 1: i + 3])))
                trows, _, _, _, _, _ = insert(v[:i + 2])
                print(Tableau.shifted_from_rows(trows), cseq(trows, v[i + 2: i + 3]))
                # print(Tableau(gamma_map(trows, v[i + 2: i + 3])))
                trows, _, _, _, _, _ = insert(v[:i + 3])
                print(Tableau.shifted_from_rows(trows), cseq(trows))
                # print(Tableau(gamma_map(trows)))
                print()
                print(utau)
                print(vtau)
                print()
                assert urows == vrows
                assert utau == vtau
                # input('?')
    return ans


def test_random_acb(bound=15):
    count = None
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        incr = help_test_acb(w)
        count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
        print(count)


@pytest.mark.slow
def test_acb(bound=7):
    count = None
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            incr = help_test_acb(w)
            count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
            print(count)
