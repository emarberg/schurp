from permutations import Permutation
from tableaux import Tableau
from words import (
    Word,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
    involution_insert
)


def commutations(iword):
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
    return ans


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
        rng = [t for t in range(len(a)) if a[t] != b[t]]
        pq = {a[j] for j in rng} | {b[j] for j in rng}
        if pq:
            assert len(pq) == 2
            p, q = tuple(pq)
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
                    assert cseq_1 == [cseq_0[0] + [eta], cseq_0[1] + [d]]

                elif p == y == z <= q:
                    c = cseq_0[1][p - 1]
                    theta = cseq_0[0][p - 1] if d + 1 == c else eta
                    cseq = [cseq_0[0][:], cseq_0[1][:]]
                    cseq[0][p - 1] = theta
                    cseq[1][p - 1] = d
                    assert cseq == cseq_1

                elif p == y < z:
                    assert p + 1 == z <= q
                    cseq = [cseq_0[0][:], cseq_0[1][:]]
                    cseq[0][p - 1], cseq[0][p] = cseq[0][p], cseq[0][p - 1]
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


def test_random_bump_differential(bound=30):
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


def help_test_disjoint_cap(a, u, v):
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


def test_disjoint_cap_simple():
    a = (8, 10, 9, 5, 4, 13, 14, 6, 2, 10, 7, 6, 12, 11, 3, 4, 5, 8, 12)
    u = 4
    v = 6
    help_test_disjoint_cap(a, u, v)

    a = (5, 8, 7, 14, 11, 10, 3, 9, 4, 1, 12, 13, 16, 8, 6, 14, 2, 11, 15, 5, 3, 7, 6, 16, 10, 8, 12, 5, 14, 9, 2, 13, 8, 7) 
    u = 6
    v = 8
    help_test_disjoint_cap(a, u, v)

    a = (1, 14, 12, 11, 13, 8, 10, 4, 7, 2, 14, 11, 6, 3, 9, 12, 7, 8, 5)
    u = 2
    v = 4
    help_test_disjoint_cap(a, u, v)

    a = (5, 1, 4, 10, 3, 12, 7, 11, 9, 6, 10, 2)
    u = 5
    v = 7
    help_test_disjoint_cap(a, u, v)

    a = (3, 9, 7)
    u = 4
    v = 8
    help_test_disjoint_cap(a, u, v)

    a = (3, 5, 6, 1, 4, 5, 6, 2, 3)
    u = 1
    v = 4
    help_test_disjoint_cap(a, u, v)


def test_random_disjoint_cap(bound=30):
    count = None
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        for i in range(len(w) - 1):
            a, u, v = w[:i], w[i], w[i + 1]
            if u + 1 < v:
                incr = help_test_disjoint_cap(a, u, v)
                count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                print(count)


def test_disjoint_cap(bound=7):
    wseen = set()
    seen = set()
    count = None
    for n in range(bound):
        pi = Permutation.longest_element(n)
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
                    incr = help_test_disjoint_cap(a, u, v)
                    count = incr if count is None else [incr[s] + t for (s, t) in enumerate(count)] if incr else count
                    print(count)


def help_test_disjoint(a, u, v):
    t1, w1, cw1, s1, cs1 = partial_insert(a, u)
    t2, w2, cw2, s2, cs2 = partial_insert(a, v)

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
        e1 = {_.number for _ in q1.entries()}
        e2 = {_.number for _ in q2.entries()}
        diag1 = {q1[box].number for box in q1 if box[0] == box[1]}
        diag2 = {q2[box].number for box in q2 if box[0] == box[1]}

        (x, y) = sorted(strict)[0]

        try:
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
            print('(1) cseq:', cseq(t1, (v,)), '==', cseq(t2, (u,)))
            assert False

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


def test_random_disjoint(bound=30):
    c1, c2, c3 = 0, 0, 0
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        for i in range(len(w) - 1):
            a = w[:i]
            u = w[i]
            v = w[i + 1]
            if u + 1 < v:
                count = help_test_disjoint(a, u, v)
                c1 += count[0]
                c2 += count[1]
                c3 += count[2]
                print(c1, c2, c3)


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
                if u + 1 >= v:
                    continue
                if (a, u, v) in wseen:
                    continue
                else:
                    wseen.add((a, u, v))
                t = Tableau.from_rows(insert(a)[0])
                if u + 1 < v and (t, u, v) not in seen:
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
    rows = insert(word[:j])[0]
    b = word[j:]
    gamma = gamma_map(rows, b)
    return [
        [gamma[(i, i)] for i in range(1, 1 + len(rows))],
        [rows[i - 1][0] for i in range(1, 1 + len(rows))],
    ]


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


def test_random_gamma(bound=30):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_gamma(w)


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


def test_random_bac(bound=30):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_bac(w)


def test_bac(bound=7):
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            help_test_bac(w)


def help_test_acb(w):
    for i in range(len(w) - 2):
        a, c, b = w[i: i + 3]
        if a < b < c:
            u = w[:i + 3]
            v = w[:i] + (c, a, b)
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


def test_random_acb(bound=30):
    for n in range(bound):
        w = Permutation.random_involution_word(n)
        help_test_acb(w)


def test_acb(bound=7):
    for n in range(bound):
        pi = Permutation.longest_element(n)
        for w in pi.get_involution_words():
            help_test_acb(w)
