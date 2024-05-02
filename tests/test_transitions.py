from transitions import grothendieck_transitions
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
    t0 = time.time()
    w0 = tuple(i + 1 for i in range(n))
    for w in itertools.permutations(w0):
        for j in range(1, n + 2):
            grothendieck_transitions(w, j)
    t1 = time.time()
    print('n =', n, 'time =', t1 - t0)