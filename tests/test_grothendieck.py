from schubert import *

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
        print(t, ':', c, ' is subset')
        assert t
    assert len(classes) == 1
    print()
    a = sorted(a, key=lambda x:(Permutation(*x).rank,Permutation(*x).length()))
    return a
    # for w in a:
    #     for v in a:
    #         if v < w:
    #             if len(v) == len(w) and len([i for i in range(len(v)) if v[i] != w[i]]) == 3:
    #                 print(v, '<->', w)

