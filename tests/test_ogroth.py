from ogroth import grothendieck_transitions
import itertools
import time


def test_grothendieck_transitions():
    w = (1, 3, 4, 5, 2)
    j = 3
    assert set(grothendieck_transitions(w, j)) == {
        ((1, 3, 4, 5, 2), 1),
        ((1, 3, 5, 4, 2), 1),
        ((1, 4, 3, 5, 2), -1),
        ((1, 4, 5, 3, 2), -1),
        ((3, 4, 1, 5, 2), 1),
        ((3, 4, 5, 1, 2), 1),
        ((3, 4, 2, 5, 1), 1),
        ((3, 4, 5, 2, 1), 1)
    }


def test_all(n):
    total = 1
    for i in range(1, n + 1):
        total *= i
    t = total // 100
    ###
    t0 = t1 = time.time()
    w0 = tuple(i + 1 for i in range(n))
    for i, w in enumerate(itertools.permutations(w0)):
        for j in range(1, n + 2):
            grothendieck_transitions(w, j)
        ###
        if n > 7 and (i + 1) % t == 0:
            print('  ', (i + 1) // t, '%', time.time() - t1)
            t1 = time.time()
        ###
    print('n =', n, 'time =', time.time() - t0)