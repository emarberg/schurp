from collections import deque

GROTHENDIECK_TRANSITIONS_CACHE = {}
O_DOMINANT_CACHE = {(): [((), 1)]}


def bruhat_cover(y, i, j):
    assert i != j
    if i > j:
        i, j = j, i
    a = y[i - 1] if i - 1 < len(y) else i
    b = y[j - 1] if j - 1 < len(y) else j
    if a > b:
        return False
    for t in range(i, j - 1):
        e = y[t] if t < len(y) else t + 1
        if a < e < b:
            return False
    return True


def multiply(y, i, j):
    assert i != j
    if i > j:
        i, j = j, i
    y += tuple(a + 1 for a in range(len(y), j))
    z = y[:i - 1] + (y[j - 1],) + y[i:j - 1] + (y[i - 1],) + y[j:]
    n = len(z)
    while n > 0 and z[n - 1] == n:
        n -= 1
    return z[:n]


def grothendieck_transitions(w, j):
    """
    Yields pairs (z, sign) where

        z = w(a_1,j)(a_2,j)···(a_p,j)(j,b_1)(j,b_2)···(j,b_q)
        sign = (-1)**p

    and a_p < ... < a_2 < a_1 < j < b_q < ... < b_2 < b_1, and the length
    increases by exactly one upon multiplication by each transposition.
    """
    if (w, j) not in GROTHENDIECK_TRANSITIONS_CACHE:
        n = max(j, len(w))
        ans = [(w, 1)]
        queue = deque([(w, j - 1, 0)])
        while queue:
            (y, i, q) = queue.popleft()
            if i <= 0:
                queue.append((y, n + 1, q))
                continue
            if i == j:
                continue
            if bruhat_cover(y, i, j):
                z = multiply(y, i, j)
                b = q + 1 if i < j else q
                ans.append((z, (-1)**b))
                queue.append((z, i - 1, b))
            queue.append((y, i - 1, q))
        GROTHENDIECK_TRANSITIONS_CACHE[w, j] = ans
    return GROTHENDIECK_TRANSITIONS_CACHE[w, j]


def grothendieck_double_transitions(w, i, j):
    (i, j) = (j, i) if i > j else (i, j)
    if (w, i, j) not in GROTHENDIECK_TRANSITIONS_CACHE:
        ans = {w: -1}
        for y, c in grothendieck_transitions(w, i):
            for z, d in grothendieck_transitions(y, j):
                ans[z] = ans.get(z, 0) + c * d
        ans = [(y, ans[y]) for y in ans if ans[y]]
        GROTHENDIECK_TRANSITIONS_CACHE[w, i, j] = ans
    return GROTHENDIECK_TRANSITIONS_CACHE[w, i, j]


def ogroth_expansion(mu):
    if mu not in O_DOMINANT_CACHE:
        j = len(mu)
        i = mu[j - 1] + j - 1

        nu = list(mu)
        nu[-1] -= 1
        if nu[-1] == 0:
            nu = nu[:-1]
        nu = tuple(nu)

        ans = {}
        for y, c in ogroth_expansion(nu):
            for z, d in grothendieck_double_transitions(y, i, j):
                ans[z] = ans.get(z, 0) + c * d
        ans = [(y, ans[y]) for y in ans if ans[y]]
        O_DOMINANT_CACHE[mu] = ans
    return O_DOMINANT_CACHE[mu]
