from schubert import *
from collections import defaultdict


def _test_topdegree(n):
    for w in Permutation.fpf_involutions(n):
        gp = FPFGrothendieck.get(w)
        print(w, gp)
        print()
    return gp


def up(w, n):
    for i in range(len(w) - 2):
        c, a, b = w[i: i + 3]
        if a < b < c and c != n + 1:
            for v in [
                w[:i] + (b, c, a) + w[i + 3:],
            ]:
                yield v


def chinese_class(w, n):
    def span(w, seen):
        if len(w) == n:
            for i in range(n // 2, len(w)):
                if max((0,) + w[i:]) <= n // 2:
                    v = w[:i] + (n + 1,) + w[i:]
                    if v not in seen:
                        seen.add(v)
                        yield v

        if len(w) == n + 1:
            v = tuple(i for i in w if i <= n)
            if v not in seen:
                seen.add(v)
                yield v

        for i in range(len(w) - 2):
            c, a, b = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (b, c, a) + w[i + 3:],
                    w[:i] + (c, b, a) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

            b, c, a = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (c, a, b) + w[i + 3:],
                    w[:i] + (c, b, a) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

            c, b, a = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (b, c, a) + w[i + 3:],
                    w[:i] + (c, a, b) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

    seen = {w}
    add = {w}
    while add:
        nextadd = set()
        for w in add:
            yield w
            nextadd |= set(span(w, seen))
        add = nextadd


def read_grothendieck(n):
    with open('gdata/g%s.py' % n, 'r') as f:
        return eval(f.read())


def _nested(w):
    y = 0
    for j in range(2, len(w) + 1):
        for i in range(len(w) - j + 1):
            if all(w[i + k] > w[i + k + 1] for k in range(j - 1)):
                if i + j >= len(w) or w[i + j - 1] < w[i + j]:
                    if i == 0 or w[i - 1] < w[i]:
                        for p in _nested(w[:i] + w[i + j:]):
                            a = ((w[i:i + j]),) + p
                            yield a
                            y += 1
    if y == 0:
        yield ()


def nested(w):
    def num(s):
        ans = 0
        while s:
            ans = 10 * ans + s[0]
            s = s[1:]
        return ans

    pairs = {
        tuple([num(t) for t in seq])
        for seq in set(_nested(w))
    }
    # b = True
    # while b:
    #     b = False
    #     remove = set()
    #     for x in pairs:
    #         for y in pairs:
    #             if x != y and set(x).issubset(set(y)):
    #                 b = True
    #                 remove.add(x)
    #     pairs -= remove
    return pairs


def longest_grothendieck(n):
    s = InvGrothendieck.top(n)
    m = Permutation.longest_element(n).involution_length()
    d = Grothendieck.decompose(s)
    f = {
        tuple(w.inverse().oneline): c * (-X(0))**(w.length() - m) for w, c in d.dictionary.items()
    }
    a = set(f)
    classes = []
    for w in a:
        if any(w in cl for cl in classes):
            continue
        c = set(chinese_class(w, n))
        classes.append(c)
        t = c.issubset(a)
        # print(t, ':', c, ' is subset')
        assert t
    assert len(classes) == 1
    print()
    d = defaultdict(int)
    a = sorted(a, key=lambda x: (f[x].substitute(0, 1), f[x].degree(), f[x], x))
    for w in a:
        if len(w) > n:
            continue
        d[f[w]] += 1
        print('  ', w, ':', f[w].set(0, 1))  # , '  :  ', Permutation(*w).get_reduced_words())
        # for v in sorted(a):
        #     if v != w and w == tuple(i for i in v if i <= n):
        #         print('  ', v, ':', f[v])
        # print()
    print(d)
    return a
