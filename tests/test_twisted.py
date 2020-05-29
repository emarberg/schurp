from permutations import Permutation
import subprocess


def test_twisted_involutions(n=5):
    w0 = Permutation.longest_element(n)
    assert set(Permutation.twisted_involutions(n)) == {
        w0 * w for w in Permutation.involutions(n)
    }
    w0.get_twisted_atoms(n)


def test_twisted(n=5):
    for w in Permutation.twisted_involutions(n):
        for a in w.get_twisted_atoms(n):
            v = a.star(n).inverse() % a
            print(w, a, a.get_reduced_word(), v)
            assert v == w
            assert a.length() == w.twisted_involution_length(n)


def _test_shape(w, n):
    y = w.star(n).inverse() % w
    assert y.twisted_involution_length(n) == w.length()

    w0 = Permutation.longest_element(n)
    y = w0 * y
    aword = reversed(w.get_reduced_word())

    # yfixed = {i for i in range(1, n + 1) if y(i) == i}
    v = Permutation()
    sh = set()
    for a in aword:
        if a > 0 and y(a) == a and y(a + 1) == a + 1:
            e, f = tuple(sorted([v(a), v(a + 1)]))
            sh |= {(e, f)}
        s = Permutation.s_i(a)
        v *= s
        z = s % y % s
        assert z.involution_length() == y.involution_length() + 1
        y = z
    f = {i for p in sh for i in p}
    return sh  # | {(i, i) for i in yfixed - f}


def test_twisted_shape(n=5):
    w0 = Permutation.longest_element(n)
    for w in Permutation.twisted_involutions(n):
        shapes = {}
        for a in w.get_twisted_atoms(n):
            sh = tuple(sorted(a.twisted_shape(n)))
            shapes[sh] = shapes.get(sh, []) + [a]
        print(w, ' :: ', (w0 * w).cycle_repr(), ' :: ', (w0 * w).fixed(n))
        print()
        for sh, atoms in shapes.items():
            minima = [a.inverse() for a in atoms if is_min_atom(a, n)]
            maxima = [a.inverse() for a in atoms if is_max_atom(a, n)]
            u = w.get_max_twisted_atom(n, sh).inverse()
            v = w.get_min_twisted_atom(n, sh).inverse()
            print(' ', set(sh), '->', v, '=', minima, '<', [a.inverse() for a in atoms], '<', maxima, '=', u)
            print()

            test = {v.inverse()} | {Permutation(*q).inverse() for (_, q, _) in span(v, n)}
            assert len(minima) == 1
            assert len(maxima) == 1
            assert minima == [v]
            assert maxima == [u]
            assert sorted(test) == sorted(atoms)
        print()
        print()
        for sh, atoms in shapes.items():
            sh = set(sh)
            for a in atoms:
                if _test_shape(a, n) != sh:
                    print(' *', a, ':', sh, '!=', _test_shape(a, n))
                    raise Exception


def is_min_atom(w, n):
    v = w.inverse().oneline
    for i in range(len(v) + 1, n + 1):
        v += (i,)
    for i in range(n):
        if i + 1 < n - 2 - i:
            d, a, b, c = v[i], v[i + 1], v[n - 2 - i], v[n - 1 - i]
            if a < d and b < c:
                return False
    return True


def is_max_atom(w, n):
    v = w.inverse().oneline
    for i in range(len(v) + 1, n + 1):
        v += (i,)
    for i in range(n):
        if i + 1 < n - 2 - i:
            a, d, c, b = v[i], v[i + 1], v[n - 2 - i], v[n - 1 - i]
            if a < d and b < c:
                return False
    return True


def span(v, n, strong=False):
    v = tuple(v.oneline)
    for i in range(len(v) + 1, n + 1):
        v += (i,)
    level = {(None, v, None)}
    while level:
        nextlevel = set()
        for u, v, label in level:
            if u is not None:
                yield (u, v, label)
            for i in range(n):
                if i + 1 < n - 2 - i:
                    a, d, c, b = v[i], v[i + 1], v[n - 2 - i], v[n - 1 - i]
                    if a < d and b < c:
                        u = list(v)
                        u[i], u[i + 1], u[n - 2 - i], u[n - 1 - i] = d, a, b, c
                        nextlevel.add((v, tuple(u), False))
            if strong:
                if n % 2 == 0 and n >= 4:
                    # 3412 -> 2431
                    k = n // 2
                    for i in range(k):
                        for j in range(i + 1, k):
                            b, d, c, a = v[i], v[j], v[n - 1 - j], v[n - 1 - i]
                            if a < b < c < d and all(x > c for x in v[i + 1:j]) and all(x > d for x in v[j + 1:k]) and all(x > c for x in v[k:n - 1 - j]) and all(x > b for x in v[n - j:n - i - 1]):
                                u = list(v)
                                u[i], u[j], u[n - 1 - j], u[n - 1 - i] = c, d, a, b
                                nextlevel.add((v, tuple(u), True))
                if n % 2 != 0 and n >= 3:
                    # 312 -> 231
                    k = n // 2
                    for j in range(k):
                        b, c, a = v[j], v[k], v[n - 1 - j]
                        if a < b < c and all(c < x for x in v[j + 1:k]) and all(b < x for x in v[k + 1:n - j - 1]):
                            u = list(v)
                            u[j], u[k], u[n - 1 - j] = c, a, b
                            nextlevel.add((v, tuple(u), True))
                    for i in range(k):
                        for j in range(i + 1, k):
                            b, d, c, a = v[i], v[j], v[n - 1 - j], v[n - 1 - i]
                            if a < b < c < d < v[k] and all(x > c for x in v[i + 1:j]) and all(x > d for x in v[j + 1:k]) and all(x > c for x in v[k + 1:n - 1 - j]) and all(x > b for x in v[n - j:n - i - 1]):
                                u = list(v)
                                u[i], u[j], u[n - 1 - j], u[n - 1 - i] = c, d, a, b
                                nextlevel.add((v, tuple(u), True))
        level = nextlevel


def print_atoms_span(n=3):
    for w in Permutation.twisted_involutions(n):
        v = w.get_min_twisted_atom(n).inverse()
        edges = list(span(v, n, True))
        if len(edges) == 0:
            continue
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [fontname="courier"];']
        for x in set(w.get_twisted_atoms(n)):
            s += ['    "%s";' % str(x.inverse())]
        s += ['    "%s" -> "%s" [style="%s"];' % (str(Permutation(*x)), str(Permutation(*y)), 'dotted' if b else 'bold') for (x, y, b) in edges]
        s += ['}']
        s = '\n'.join(s)
        name = ''.join([str(v(i)) for i in range(1, n + 1)])
        file = '/Users/emarberg/Desktop/examples/atoms/'
        dotfile = file + 'dot/AIII/' + name + '.dot'
        pngfile = file + 'png/AIII/' + name + '.png'
        with open(dotfile, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])


def test_atoms_span(n=5):
    cls = Permutation
    for w in cls.twisted_involutions(n):
        v = w.get_min_twisted_atom(n).inverse()
        test = sorted({v} | {cls(*u) for (_, u, _) in span(v, n, True)})
        sest = sorted([u.inverse() for u in w.get_twisted_atoms(n)])
        if test != sest:
            print(w)
            print('  ', v, '->', test)
            x = []
            for u in sest:
                if is_min_atom(u.inverse(), n) and u not in test:
                    x += ['*']
                x += [u]
            print('  ', x)
            assert False
