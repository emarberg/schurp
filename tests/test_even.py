from signed import SignedPermutation
from even import EvenSignedPermutation
import subprocess


def test_is_perfect(m=5):
    def is_perfect(w, n):
        for t in EvenSignedPermutation.reflections(n):
            wtw = w * t.star() * w.star()
            if t * wtw != wtw * t:
                return False
        return True

    for n in range(3, m + 1):
        seen = {}
        for w in EvenSignedPermutation.twisted_involutions(n):
            if is_perfect(w, n):
                seen[abs(w)] = seen.get(abs(w), []) + [w]

        for w, enum in seen.items():
            print(w.cycle_repr(), enum)
        assert len(seen) == 1 and list(seen)[0].is_identity()


def test_get_minimal_fpf_involution():
    cls = EvenSignedPermutation
    assert cls.get_minimal_fpf_involution(1) == cls(1)
    assert cls.get_minimal_fpf_involution(2) == cls(2, 1)
    assert cls.get_minimal_fpf_involution(3) == cls(1, 3, 2)
    assert cls.get_minimal_fpf_involution(4) == cls(2, 1, 4, 3)


def test_get_involution_word():
    s = EvenSignedPermutation.s_i(0, 3)
    t = EvenSignedPermutation.s_i(2, 3)
    assert (t * s * t).get_involution_word() in {(0, 2), (2, 0)}
    assert (s * t * s).get_involution_word() in {(0, 2), (2, 0)}
    assert set((t * s * t).get_involution_words()) == {(0, 2), (2, 0)}
    assert set((s * t * s).get_involution_words()) == {(0, 2), (2, 0)}


def test_fpf_involution_words():
    for n in [2, 3]:
        for w in EvenSignedPermutation.fpf_involutions(n):
            words = set(w.get_fpf_involution_words())
            if words:
                print(w, '->', len(words))


def test_get_atoms():
    for n in [3, 4]:
        i = list(EvenSignedPermutation.involutions(n))
        for w in i:
            words = set()
            for a in w.get_atoms():
                assert len(a) == w.involution_length()
                assert a.inverse() % a == w
                words |= set(a.get_reduced_words())
            assert words == set(w.get_involution_words())


def test_get_twisted_atoms():
    for n in [3, 4]:
        i = list(EvenSignedPermutation.twisted_involutions(n))
        for w in i:
            words = set()
            for a in w.get_twisted_atoms():
                assert len(a) == w.twisted_involution_length()
                assert a.inverse().star() % a == w
                words |= set(a.get_reduced_words())
            assert words == set(w.get_twisted_involution_words())


def test_twisted_involutions():
    # n must be odd
    for n in [3, 5]:
        w0 = EvenSignedPermutation.longest_element(n)
        assert set(EvenSignedPermutation.twisted_involutions(n)) == {
            w0 * w for w in EvenSignedPermutation.involutions(n)
        }


def test_twisted():
    for n in [4, 5]:
        for w in EvenSignedPermutation.twisted_involutions(n):
            for a in w.get_twisted_atoms():
                v = a.star().inverse() % a
                print(w, a, a.get_reduced_word(), v)
                assert v == w
                assert a.length() == w.twisted_involution_length()


def _nneg(w):
    des = 0
    for i in range(len(w) - 1):
        a, b = w[i:i + 2]
        if a > b:
            des += 1
            v = w[:i] + w[i + 2:]
            for nn in _nneg(v):
                yield {(abs(a), b)} | nn
    if des == 0:
        yield {w}


def test_nneg(nn=3):
    for n in [nn, nn + 1]:
        for w in EvenSignedPermutation.involutions(n):
            for atom in w.get_atoms():
                o = atom.inverse().oneline
                nn = {tuple(sorted(x)) for x in _nneg(o)}

                ndes, fix = atom.ndes(), atom.nfix()
                neg = []
                for a, b in ndes:
                    if a < -b:
                        neg += [a, -b]
                neg = tuple(sorted(neg))

                sh = {(a, -b) for a, b in ndes if 0 < a < -b}
                sh |= {(b, -a) for a, b in ndes if 0 < a < -b}
                sh = tuple(sorted(sh))

                pair = sorted((b, a) for a, b in ndes if not (0 < a < -b))

                cfix = tuple(i for i in w.fixed_points() if i > 0)
                cneg = tuple(i for i in w.negated_points() if i > 0)
                csh = tuple(sorted(atom.shape()))
                cpair = w.pair()

                print(w, ':', list(o))
                print('         shape:', csh, '==', sh)
                print('         pairs:', cpair, '==', pair)
                print('  fixed points:', cfix, '==', fix)
                print('negated points:', cneg, '==', neg)
                print()
                assert fix == cfix
                assert neg == cneg
                assert len(nn) == 1
                assert csh == sh
                assert cpair == pair


def test_twisted_nneg(nn=3):
    for n in [nn, nn + 1]:
        for w in EvenSignedPermutation.twisted_involutions(n):
            for atom in w.get_twisted_atoms():
                o = atom.inverse().oneline[1:]
                nn = {tuple(sorted(x)) for x in _nneg(o)}

                ndes, fix = atom.twisted_ndes(), atom.twisted_nfix()
                neg = []
                for a, b in ndes:
                    if a < -b:
                        neg += [a, -b]
                    elif a == -b:
                        neg += [a]
                neg = tuple(sorted(neg))

                sh = {(-a, a) for a, b in ndes if a + b == 0}
                sh |= {(a, -b) for a, b in ndes if 0 < a < -b}
                sh |= {(b, -a) for a, b in ndes if 0 < a < -b}
                sh = tuple(sorted(sh))

                pair = sorted((b, a) for a, b in ndes if not (0 < a <= -b))

                cfix = tuple(i for i in w.twisted_fixed_points() if i > 0)
                cneg = tuple(i for i in w.twisted_negated_points() if i > 0)
                csh = tuple(sorted(atom.twisted_shape()))
                cpair = w.pair()

                print(w, ':', list(o))
                print('         shape:', csh, '==', sh)
                print('         pairs:', cpair, '==', pair)
                print('  fixed points:', cfix, '==', fix)
                print('negated points:', cneg, '==', neg)
                print()
                assert fix == cfix
                assert neg == cneg
                assert len(nn) == 1
                assert csh == sh
                assert cpair == pair


def is_max_atom(w):
    v = w.inverse().oneline
    for i in range(len(v) - 2):
        b, c, a = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 3:
        b, c, a = v[:3]
        if a < -b < c:
            return False
    return True


def is_min_atom(w):
    v = w.inverse().oneline
    for i in range(len(v) - 2):
        c, a, b = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 3:
        c, a, b = v[:3]
        if a < b < -c:
            return False
    return True


def span(v, strong=False):
    v = v.oneline
    level = {(None, v, None)}
    while level:
        nextlevel = set()
        for u, v, label in level:
            if u is not None:
                yield (u, v, label)
            for i in range(len(v) - 2):
                b, c, a = v[i: i + 3]
                if a < b < c:
                    w = v[:i] + (c, a, b) + v[i + 3:]
                    nextlevel.add((v, w, False))
            if len(v) >= 3:
                b, c, a = v[:3]
                if a < -b < c:
                    w = (-c, a, -b) + v[3:]
                    nextlevel.add((v, w, False))
            if strong:
                for j in range(1, len(v)):
                    for k in range(j + 1, len(v) - 1):
                        b, c, d = v[j], v[k], v[k + 1]
                        if 0 < -b < c < -d and all(abs(x) <= abs(b) for x in v[:k]):
                            w = v[:j] + (d,) + v[j + 1:k] + (-b, -c) + v[k + 2:]
                            nextlevel.add((v, w, True))
        level = nextlevel


def print_atoms_span(n=3):
    def printer(oneline):
        w = SignedPermutation(*oneline.oneline) if type(oneline) == EvenSignedPermutation else SignedPermutation(*oneline)
        sh = w.inverse().dshape()
        return str(w) + '\n' + str(sh)
    
    cls = EvenSignedPermutation
    for w in cls.involutions(n):
        v = w.get_min_atom().inverse()
        edges = list(span(v, True))
        atoms = {a.inverse() for a in w.get_atoms()}
        
        if len(edges) == 0:
            continue

        sources = {x for x in atoms if not any(tuple(x) == tuple(v) for u, v, b in edges)}
        assert len(sources) == 1
        assert next(iter(sources)) == v
        expected = {x for x in atoms if not any(not b and tuple(x) == tuple(v) for u, v, b in edges)}
        minima = {x for x in atoms if x.inverse() == w.get_min_atom(SignedPermutation(*x).inverse().dshape(offset=0))}
        assert expected == minima
        
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        for x in atoms:
            if x in minima:
                s += ['    "%s";' % printer(x)]
            else:
                s += ['    "%s" [fillcolor=white];' % printer(x)]
        s += ['    "%s" -> "%s" [style="%s"];' % (printer(x), printer(y), 'dotted' if b else 'bold') for (x, y, b) in edges]
        s += ['}']
        s = '\n'.join(s)

        name = ''.join([str(w(i)) for i in range(1, n + 1)])
        name = 'n' + str(n) + '_' + str(len(atoms)) + '_' + name

        file = '/Users/emarberg/examples/atoms/'
        dotfile = file + 'dot/DI/' + name + '.dot'
        pngfile = file + 'png/DI/' + name + '.png'
        with open(dotfile, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])
        #subprocess.run(["open", pngfile])
        #input('')


def test_atoms_span(nn=4):
    def shape(u):
        return SignedPermutation(*u).inverse().dshape()

    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.involutions(n):
            v = w.get_min_atom().inverse()
            test = sorted({v} | {cls(*u) for (_, u, _) in span(v, True)})
            sest = sorted([u.inverse() for u in w.get_atoms()])
            if test != sest:
                print(w)
                print('  ', test)
                x = []
                for u in sest:
                    if is_min_atom(u.inverse()) and u not in test:
                        x += ['*']
                    x += [u]
                print('  ', x)
                assert False
            assert all(b == (shape(u) != shape(x)) for (u, x, b) in span(v, True))    


def test_shape(nn=4):
    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.involutions(n):
            shapes = {}
            for a in w.get_atoms():
                sh = tuple(sorted(a.shape()))
                shapes[sh] = shapes.get(sh, []) + [a]

                _a = SignedPermutation(*a.oneline)
                _sh = tuple(sorted(_a.dshape()))
                assert sh == _sh     
            for sh, atoms in shapes.items():
                minima = [a.inverse() for a in atoms if is_min_atom(a)]
                maxima = [a.inverse() for a in atoms if is_max_atom(a)]
                v = w.get_min_atom(sh).inverse()
                test = {v.inverse()} | {cls(*q).inverse() for (_, q, _) in span(v)}
                assert sorted(test) == sorted(atoms)

                assert not any(i == j for i, j in sh)
                assert minima == [v]


def is_max_twisted_atom(w):
    v = w.inverse().oneline
    for i in range(1, len(v) - 2):
        b, c, a = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 4:
        b, c, a = v[1:4]
        if a < -b < c:
            return False
    if len(v) >= 2:
        b, a = v[:2]
        if abs(a) < abs(b) and b > 0:
            return False
    return True


def is_min_twisted_atom(w):
    v = w.inverse().oneline
    for i in range(1, len(v) - 2):
        c, a, b = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 4:
        c, a, b = v[1:4]
        if a < b < -c:
            return False
    if len(v) >= 2:
        b, a = v[:2]
        if abs(a) < abs(b) and b < 0:
            return False
    return True


def twisted_span(v, strong=False):
    v = v.oneline
    level = {(None, v, None)}
    seen = set()
    while level:
        nextlevel = set()
        for u, v, label in level:
            if u is not None:
                yield (u, v, label)
            
            if v in seen:
                continue
            seen.add(v)
            
            for i in range(1, len(v) - 2):
                b, c, a = v[i: i + 3]
                if a < b < c and not (i == 1 and abs(c) < -v[0]):
                    w = v[:i] + (c, a, b) + v[i + 3:]
                    nextlevel.add((v, w, False))

            if len(v) >= 2:
                b, a = v[:2]
                if abs(a) < abs(b) == b:
                    w = (-b, -a) + v[2:]
                    nextlevel.add((v, w, False))
            if len(v) >= 4:
                x, b, c, a = v[:4]
                if a < b < c and abs(c) < -x:
                    w = (-x, -c, a, b) + v[4:]
                    nextlevel.add((v, w, False))
            if strong:
                for j in [0]:
                    for k in range(j + 1, len(v) - 1):
                        a, b, c = v[j], v[k], v[k + 1]
                        if 0 < abs(a) < b < -c and all(x <= a for x in v[:k]):
                            w = v[:j] + (-c,) + v[j + 1:k] + (a, -b) + v[k + 2:]
                            nextlevel.add((v, w, True))

                for j in range(2, len(v)):
                    for k in range(j + 1, len(v) - 1):
                        b, c, d = v[j], v[k], v[k + 1]
                        if 0 < -b < c < -d and all(abs(x) <= abs(b) for x in v[:k]):
                            w = v[:j] + (d,) + v[j + 1:k] + (-b, -c) + v[k + 2:]
                            nextlevel.add((v, w, True))

        level = nextlevel


def test_twisted_shape(nn=4):
    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.twisted_involutions(n):
            shapes = {}
            for a in w.get_twisted_atoms():
                sh = tuple(sorted(a.twisted_shape()))
                shapes[sh] = shapes.get(sh, []) + [a]

                _a = SignedPermutation(*a.oneline)
                _sh = tuple(sorted(_a.dshape(1)))
                assert sh == _sh
            for sh, atoms in shapes.items():
                minima = [a.inverse() for a in atoms if is_min_twisted_atom(a)]
                maxima = [a.inverse() for a in atoms if is_max_twisted_atom(a)]
                v = w.get_min_twisted_atom(sh).inverse()
                test = {v.inverse()} | {cls(*u).inverse() for (_, u, _) in twisted_span(v)}
                assert sorted(test) == sorted(atoms)
                assert len([i for i, j in sh if i == -j]) == 1
                assert len(minima) == 1
                assert minima == [v]


def print_twisted_atoms_span(n):
    def printer(oneline):
        w = SignedPermutation(*oneline.oneline) if type(oneline) == EvenSignedPermutation else SignedPermutation(*oneline)
        sh = w.inverse().dshape(offset=1)
        return str(w) + '\n' + str(sh)

    cls = EvenSignedPermutation
    for w in cls.twisted_involutions(n):
        v = w.get_min_twisted_atom().inverse()
        
        edges = list(twisted_span(v, True))
        atoms = {a.inverse() for a in w.get_twisted_atoms()}

        sources = {x for x in atoms if not any(tuple(x) == tuple(v) for u, v, b in edges)}
        assert len(sources) == 1
        assert next(iter(sources)) == v
        expected = {x for x in atoms if not any(not b and tuple(x) == tuple(v) for u, v, b in edges)}
        minima = {x for x in atoms if x.inverse() == w.get_min_twisted_atom(SignedPermutation(*x).inverse().dshape(offset=1))}
        assert expected == minima

        for x, y, b in edges:
            if cls(*x) not in atoms or cls(*y) not in atoms:
                print(x, y, b)
                print(atoms)
                return
        
        if len(edges) == 0:
            continue
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        for x in atoms:
            if x in minima:
                s += ['    "%s";' % printer(x)]
            else:
                s += ['    "%s" [fillcolor=white];' % printer(x)]
        s += ['    "%s" -> "%s" [style="%s"];' % (printer(x), printer(y), 'dotted' if b else 'bold') for (x, y, b) in edges]
        s += ['}']
        s = '\n'.join(s)
        name = ''.join([str(v(i)) for i in range(1, n + 1)])

        name = 'n' + str(n) + '_' + str(len(atoms)) + '_' + name
        
        file = '/Users/emarberg/examples/atoms/'
        dotfile = file + 'dot/DII/' + name + '.dot'
        pngfile = file + 'png/DII/' + name + '.png'
        with open(dotfile, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])


def test_twisted_atoms_span(nn=4):
    def shape(u):
        return SignedPermutation(*u).inverse().dshape(1)

    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.twisted_involutions(n):
            v = w.get_min_twisted_atom().inverse()
            test = sorted({v} | {cls(*u) for (_, u, _) in twisted_span(v, True)})
            sest = sorted([u.inverse() for u in w.get_twisted_atoms()])
            if test != sest:
                print(w)
                print('  ', v, '<', test)
                x = []
                for u in sest:
                    if is_min_twisted_atom(u.inverse()):
                        x += ['*']
                    x += [u]
                print('  ', v, '<', x)
                print('  ', [u for u in test if u not in sest])
                assert False
            assert all(b == (shape(u) != shape(x)) for (u, x, b) in twisted_span(v, True))    


def is_max_fpf_atom(w):
    v = w.inverse().oneline
    for i in range(len(v) - 2):
        b, c, a = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 3:
        b, c, a = v[:3]
        if a < -b < c:
            return False
    return True


def is_min_fpf_atom(w):
    v = w.inverse().oneline
    for i in range(len(v) - 2):
        c, a, b = v[i: i + 3]
        if a < b < c:
            return False
    if len(v) >= 3:
        c, a, b = v[:3]
        if a < b < -c:
            return False
    return True


def fpf_span(v, strong=False):
    v = v.oneline
    level = {(None, v, None)}
    while level:
        nextlevel = set()
        for u, v, label in level:
            if u is not None:
                yield (u, v, label)
            for i in range(len(v) - 2):
                b, c, a = v[i: i + 3]
                if a < b < c:
                    w = v[:i] + (c, a, b) + v[i + 3:]
                    nextlevel.add((v, w, False))
            if len(v) >= 3:
                b, c, a = v[:3]
                if a < -b < c:
                    w = (-c, a, -b) + v[3:]
                    nextlevel.add((v, w, False))
            if strong:
                for j in range(1, len(v)):
                    for k in range(j + 1, len(v) - 1):
                        b, c, d = v[j], v[k], v[k + 1]
                        if 0 < -b < c < -d and all(abs(x) <= abs(b) for x in v[:k]):
                            w = v[:j] + (d,) + v[j + 1:k] + (-b, -c) + v[k + 2:]
                            nextlevel.add((v, w, True))
        level = nextlevel


def print_fpf_atoms_span(n=3):
    def printer(oneline):
        w = SignedPermutation(*oneline.oneline) if type(oneline) == EvenSignedPermutation else SignedPermutation(*oneline)
        sh = w.inverse().fpf_dshape()
        return str(w) + '\n' + str(sh)
    
    cls = EvenSignedPermutation
    for w in cls.fpf_involutions(n):
        v = w.get_min_fpf_atom().inverse()
        edges = list(span(v, True))
        atoms = {a.inverse() for a in w.get_fpf_atoms()}
        
        if len(edges) == 0:
            continue

        sources = {x for x in atoms if not any(tuple(x) == tuple(v) for u, v, b in edges)}
        assert len(sources) == 1
        assert next(iter(sources)) == v
        expected = {x for x in atoms if not any(not b and tuple(x) == tuple(v) for u, v, b in edges)}
        minima = {x for x in atoms if x.inverse() == w.get_min_fpf_atom(SignedPermutation(*x).inverse().fpf_dshape())}
        assert expected == minima
        
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        for x in atoms:
            if x in minima:
                s += ['    "%s";' % printer(x)]
            else:
                s += ['    "%s" [fillcolor=white];' % printer(x)]
        s += ['    "%s" -> "%s" [style="%s"];' % (printer(x), printer(y), 'dotted' if b else 'bold') for (x, y, b) in edges]
        s += ['}']
        s = '\n'.join(s)

        name = ''.join([str(w(i)) for i in range(1, n + 1)])
        name = 'n' + str(n) + '_' + str(len(atoms)) + '_' + name

        file = '/Users/emarberg/examples/atoms/'
        dotfile = file + 'dot/DIII/' + name + '.dot'
        pngfile = file + 'png/DIII/' + name + '.png'
        with open(dotfile, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])
        #subprocess.run(["open", pngfile])
        #input('')


def test_fpf_atoms_span(nn=4):
    def shape(u):
        return SignedPermutation(*u).inverse().fpf_dshape()

    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.fpf_involutions(n):
            v = w.get_min_fpf_atom().inverse()
            test = sorted({v} | {cls(*u) for (_, u, _) in span(v, True)})
            sest = sorted([u.inverse() for u in w.get_fpf_atoms()])
            if test != sest:
                print(w)
                print('  ', test)
                x = []
                for u in sest:
                    if is_min_fpf_atom(u.inverse()) and u not in test:
                        x += ['*']
                    x += [u]
                print('  ', x)
                assert False
            assert all(b == (shape(u) != shape(x)) for (u, x, b) in span(v, True))    


def test_fpf_shape(nn=4):
    for n in [nn, nn + 1]:
        cls = EvenSignedPermutation
        for w in cls.fpf_involutions(n):
            shapes = {}
            for a in w.get_fpf_atoms():
                sh = tuple(sorted(a.fpf_shape()))
                shapes[sh] = shapes.get(sh, []) + [a]

                _a = SignedPermutation(*a.oneline)
                _sh = tuple(sorted(_a.fpf_dshape()))
                assert sh == _sh     
            for sh, atoms in shapes.items():
                minima = [a.inverse() for a in atoms if is_min_fpf_atom(a)]
                maxima = [a.inverse() for a in atoms if is_max_fpf_atom(a)]
                v = w.get_min_fpf_atom(sh).inverse()
                test = {v.inverse()} | {cls(*q).inverse() for (_, q, _) in fpf_span(v)}
                assert sorted(test) == sorted(atoms)
                assert not any(i == j for i, j in sh)
                assert minima == [v]