from permutations import Permutation
from tableaux import Tableau
from words import (
    Word,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
    involution_insert
)


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


def test_batch_insert(n=5):
    for w in Permutation.involutions(n):
        for word in w.get_primed_involution_words():
            print(w, word)
            p, _, _ = primed_insert(word)
            for k in range(len(word)):
                tab, _, _ = primed_insert(word[:k])
                a = word[k:]
                q = batch_insert_sequential(tab, a)
                assert p == q
                r = batch_insert_parallel(tab, a)
                assert q == r


def extract_all_cycles(rows, rest):
    ans = []
    for i in range(len(rows)):
        ans.append(extract_cycle(rows[i][0], rows[:i + 1], rest))
    return ans


def insert(iword):
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
        for j in range(len(a)):
            p, q = a[j], b[j]
            permuted_cycles[p], permuted_cycles[q] = permuted_cycles[q], permuted_cycles[p]

    diagonal_cycles = {extract_cycle(rows[i][0], rows[:i + 1], []) for i in range(len(rows))}
    cycle_map = extract_commutation_positions(rows)
    return rows, diagonal_cycles, permuted_cycles, cycle_map, tabs, paths


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


def test_bumping_path(bound=5):
    for n in range(bound):
        for pi in Permutation.involutions(n):
            for w in pi.get_involution_words():
                _test_bumping_path(w)
