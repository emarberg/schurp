from permutations import Permutation


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


def extract_cycle(a, rows):
    p = {a, a + 1}
    for row in reversed(rows):
        for i in row:
            s = Permutation.s_i(i)
            p = {s(x) for x in p}
    return tuple(sorted(p))


def insert(iword):
    rows = []
    diagonal_cycles = set()
    permuted_cycles = {}
    for n, a in enumerate(iword):
        transposed = False
        a = abs(a)
        i = 0
        while i < len(rows):
            j = [t for t in range(len(rows[i])) if a <= rows[i][t]]
            if j:
                j = j[0]
                ####
                if j == 0 and not transposed and a + 1 < rows[i][j]:
                    p = extract_cycle(rows[i][j], rows[:i + 1])
                    q = extract_cycle(a, rows[:i])
                if j == 0 and not transposed and a == rows[i][j]:
                    p = extract_cycle(a, rows[:i] + [rows[i] + [a]])
                    q = extract_cycle(a + 2, rows[:i] + [rows[i] + [a]] + [rows[i + 1]])
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
        ####
        if a is not None and not transposed:
            diagonal_cycles.add(extract_cycle(a, rows))
        ####
    return rows, diagonal_cycles, permuted_cycles
