from stable.tableaux import Tableau


def row_bump(tab, i, x):
    j = 1
    while tab.get(i, j, default=x + 1) <= x:
        j += 1
    if (i, j) not in tab:
        return tab.set(i, j, x), i
    y, tab = tab.get(i, j), tab.set(i, j, x)
    return row_bump(tab, i + 1, y)


def reverse_row_bump(tab, i_start, i_stop, x):
    assert 1 <= i_stop <= i_start
    if i_start == i_stop:
        return tab, x
    i = i_start
    j = 0
    while tab.get(i, j + 1, default=x) < x:
        j += 1
    assert j > 0
    y, tab = tab.get(i, j), tab.set(i, j, x)
    return reverse_row_bump(tab, i - 1, i_stop, y)


def row_extract_and_bump(tab, start=None):
    if start is None:
        start = max([(i, j) for i, j, v in tab if len(v) > 1])
    i_start, j_start = start
    g = tab.get(i_start, j_start, unpack=False)
    assert len(g) > 1
    x = max(g)
    tab = tab.remove(i_start, j_start, x)
    i = i_start + 1
    return row_bump(tab, i, x)


def reverse_row_extract_and_bump(tab, i_start, i_stop):
    assert 1 <= i_stop < i_start
    j_start = tab.max_column(i_start)
    g = tab.get(i_start, j_start, unpack=False)
    assert len(g) == 1

    x = max(g)
    tab = tab.remove(i_start, j_start, x)
    i = i_start - 1
    tab, x = reverse_row_bump(tab, i, i_stop, x)

    i = i_stop
    j = 0
    while max(tab.get(i, j + 1, unpack=False, default=x)) < x:
        j += 1
    assert j > 0
    return tab.add(i, j, x), j


def slide(tab, i, j):
    assert not tab.contains_box(i, j)
    a = tab.get(i + 1, j)
    b = tab.get(i, j + 1)

    if a is None and b is None:
        return (tab, i, j)
    if a is None:
        tab = tab.add(i, j, b).remove(i, j + 1)
        return slide(tab, i, j + 1)
    if b is None:
        tab = tab.add(i, j, a).remove(i + 1, j)
        return slide(tab, i + 1, j)

    assert type(a) == type(b) == int
    if a > b:
        tab = tab.add(i, j, b).remove(i, j + 1)
        return slide(tab, i, j + 1)
    else:
        tab = tab.add(i, j, a).remove(i + 1, j)
        return slide(tab, i + 1, j)


def reverse_slide(tab, i_start, i_stop, j=None):
    assert 1 <= i_stop <= i_start
    if j is None:
        j = tab.max_column(i_start) + 1
    assert not tab.contains_box(i_start, j)
    if i_stop == i_start:
        return tab, j
    i = i_start

    a = tab.get(i - 1, j)
    b = tab.get(i, j - 1)
    assert a is not None or b is not None

    if a is None:
        tab = tab.add(i, j, b).remove(i, j - 1)
        return reverse_slide(tab, i, i_stop, j - 1)
    if b is None:
        tab = tab.add(i, j, a).remove(i - 1, j)
        return reverse_slide(tab, i - 1, i_stop, j)
    if a < b:
        tab = tab.add(i, j, b).remove(i, j - 1)
        return reverse_slide(tab, i, i_stop, j - 1)
    else:
        tab = tab.add(i, j, a).remove(i - 1, j)
        return reverse_slide(tab, i - 1, i_stop, j)


def extract_and_slide(tab, start=None):
    if start is None:
        start = max([(i, j) for i, j, _ in tab if tab.get(i, j, unpack=True) is tab.get(i + 1, j, unpack=True)])
    i, j = start
    return slide(tab.remove(i, j), i, j)


def test_row_bump():
    p = Tableau("""
        1   1,2 2   3,4 4 6
        2,3 3   3,4 5   6 8
        3   4   5   6   7
        4   6   7   8
        6                   """)  # noqa
    q = Tableau("""
        1   1,2 2   3,4 4 6
        2,3 3   3 5   6 8
        3   4   4   6   7
        4   5   7   8
        6   6               """)  # noqa
    assert row_extract_and_bump(p) == (q, 5)
    assert reverse_row_extract_and_bump(q, 5, 2) == (p, 3)


def test_slide():
    p = Tableau("""
        1 1 2 4 5 5
        1 2 3 4 4 6
        2 2 6 6 7
        3 4 7 9
        5 """)
    x = p.remove(2, 2)
    q = Tableau("""
        1 1 2 4 5 5
        1 2 3 4 4 6
        2 4 6 6 7
        3 7 9
        5 """)
    assert slide(x, 2, 2) == (q, 4, 4)
    assert extract_and_slide(p) == (q, 4, 4)
    assert reverse_slide(q, 4, 2) == (x, 2)
