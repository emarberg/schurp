from permutations import Permutation


def fpfdes(w):
    return len([i for i in w.right_descent_set if w(i) != i + 1])


def ides(w):
    return len(w.right_descent_set)
    # ans = set()
    # for i in w.right_descent_set:
    #     s = Permutation.s_i(i)
    #     ws = w * s
    #     if s * ws != w:
    #         ws = s * ws
    #     ans.add(ws)
    # return len(ans)


def fpf_weak_order_intervals(n):
    assert n % 2 == 0
    g = [Permutation.s_i(i) for i in range(1, n)]
    e = Permutation()
    for i in range(1, n, 2):
        e *= g[i - 1]
    lower = {e: {e}}
    level = {e}
    while level:
        newlevel = set() 
        for w in level:
            for i in range(1, n):
                if w(i) < w(i + 1):
                    s = g[i - 1]
                    sws = s * w * s
                    if sws not in lower:
                        lower[sws] = {sws}
                    lower[sws] |= lower[w]
                    newlevel.add(sws)
        level = newlevel
    return lower


def inv_weak_order_intervals(n):
    g = [Permutation.s_i(i) for i in range(1, n)]
    e = Permutation()
    lower = {e: {e}}
    level = {e}
    while level:
        newlevel = set() 
        for w in level:
            for i in range(1, n):
                if w(i) < w(i + 1):
                    s = g[i - 1]
                    ws = w * s
                    if s * ws != w:
                        ws = s * ws
                    if ws not in lower:
                        lower[ws] = {ws}
                    lower[ws] |= lower[w]
                    newlevel.add(ws)
        level = newlevel
    return lower


def fpfexpectx(n):
    ans = {}
    lower = fpf_weak_order_intervals(n)
    for w in lower:
        e = 0
        for u in lower[w]:
            e += fpfdes(u)
        ans[w] = (e, len(lower[w]))
    return ans


def fpfexpecty(n):
    ans = {}
    for w in Permutation.fpf_involutions(n):
        q = (w.fpf_involution_length() + 1) * len(w.get_fpf_involution_words())
        p = 0
        for word in w.get_fpf_involution_words():
            for i in range(len(word) + 1):
                u = Permutation.from_fpf_involution_word(*word[:i])
                p += fpfdes(u)
        ans[w] = (p, q)
    return ans


def invexpectx(n):
    ans = {}
    lower = inv_weak_order_intervals(n)
    for w in lower:
        e = 0
        for u in lower[w]:
            e += ides(u)
        ans[w] = (e, len(lower[w]))
    return ans


def invexpecty(n):
    ans = {}
    for w in Permutation.involutions(n):
        q = (w.involution_length() + 1) * len(w.get_involution_words())
        p = 0
        for word in w.get_involution_words():
            for i in range(len(word) + 1):
                u = Permutation.from_involution_word(*word[:i])
                p += ides(u)
        ans[w] = (p, q)
    return ans


def right_weak_order_intervals(n):
    g = [Permutation.s_i(i) for i in range(1, n)]
    e = Permutation()
    lower = {e: {e}}
    level = {e}
    while level:
        newlevel = set() 
        for w in level:
            for i in range(1, n):
                if w(i) < w(i + 1):
                    ws = w * g[i - 1]
                    if ws not in lower:
                        lower[ws] = {ws}
                    lower[ws] |= lower[w]
                    newlevel.add(ws)
        level = newlevel
    return lower


def expectx(n):
    ans = {}
    lower = right_weak_order_intervals(n)
    for w in lower:
        e = 0
        for u in lower[w]:
            e += len(u.right_descent_set)
        ans[w] = (e, len(lower[w]))
    return ans


def expecty(n):
    ans = {}
    for w in Permutation.all(n):
        q = (w.length() + 1) * len(w.get_reduced_words())
        p = 0
        for word in w.get_reduced_words():
            for i in range(len(word) + 1):
                u = Permutation.from_word(*word[:i])
                p += len(u.right_descent_set)
        ans[w] = (p, q)
    return ans


def test_weak_order(n):
    def bprint(*args):
        args = [a if type(a) != bool else str(a) + (' ' if a else '') for a in args]
        print(*args)

    x = expectx(n)
    y = expecty(n)

    def cde(w):
        p1, q1 = x[w]
        p2, q2 = y[w]
        return p1 * q2 == p2 * q1

    for w in sorted(x, key=lambda z: (cde(z), z.shape())):
        if w.is_vexillary():
            bprint(cde(w), w.oneline_repr(n), w.is_vexillary(), w.is_grassmannian(), w.inverse().is_grassmannian(), w.is_dominant(), w.shape(), w.shape().is_balanced())
            assert cde(w) == w.shape().is_balanced()
            

def test_inv_weak_order(n):
    def bprint(*args):
        args = [a if type(a) != bool else str(a) + (' ' if a else '') for a in args]
        print(*args)

    x = invexpectx(n)
    y = invexpecty(n)

    def cde(w):
        p1, q1 = x[w]
        p2, q2 = y[w]
        return p1 * q2 == p2 * q1

    for w in sorted(x, key=lambda z: (cde(z), z.involution_shape())):
        if w.is_vexillary():
            bprint(cde(w), w.oneline_repr(n), w.is_vexillary(), w.is_inv_grassmannian(), w.is_dominant(), w.involution_shape(), w.involution_shape().is_shifted_balanced() or w.involution_shape().is_trapezoidal())


def test_fpf_weak_order(n):
    def bprint(*args):
        args = [a if type(a) != bool else str(a) + (' ' if a else '') for a in args]
        print(*args)

    x = fpfexpectx(n)
    y = fpfexpecty(n)

    def cde(w):
        p1, q1 = x[w]
        p2, q2 = y[w]
        return p1 * q2 == p2 * q1

    for w in sorted(x, key=lambda z: (cde(z), z.fpf_involution_shape())):
        bprint(cde(w), w, w.is_fpf_dominant(), w.fpf_involution_shape(), w.fpf_involution_shape().is_shifted_balanced() or w.fpf_involution_shape().is_trapezoidal())
