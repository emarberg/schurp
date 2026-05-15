from signed import SignedPermutation
from tests.test_crystals import draw_graph
from tests.test_brion import shape_CII


def vzero(n):
    assert n >= 0
    if n == 0:
        return SignedPermutation()
    elif n % 2 == 0:
        a = vzero(n - 2).inflate(n)
        return a * SignedPermutation.t_ij(n - 1, -n, n)
    else:
        a = vzero(n - 1).inflate(n)
        return a * SignedPermutation.t_ij(n, -n, n)


def vzero_plus(n):
    if n == 1:
        return SignedPermutation(1)
    if n % 4 == 0:
        return vzero(n)
    elif n % 4 == 2:
        return vzero(n) * SignedPermutation.t_ij(1, -1, n) * SignedPermutation.t_ij(2, -2, n)
    else:
        return vzero(n) * SignedPermutation.t_ij(2, -2, n)


colorfn = lambda x: 'white' if x.length() % 2 ==0 else 'gray'


def test(n, levels):
    w = SignedPermutation.longest_element(n)
    atoms = {v.reduce() for v in w.get_atoms()}

    sources = atoms.copy()
    gv = set()
    ge = set()
    while levels > 0:
        add = set()
        for sigma in SignedPermutation.all(n + 1):
            for j in [1]:
                top, bot, v, e = sigma.transition_graph(j, flip=True)
                if set(bot) and set(bot).issubset(sources):
                    print('* ', sigma, j, bot)
                    add |= set(top)
                    gv |= v
                    ge |= e
        levels -= 1
        sources |= add

    draw_graph(gv, ge)

    w0 = SignedPermutation.longest_element

    print()
    print('q  =', w0(n - 1), ':', {v.reduce() for v in w0(n - 1).get_atoms()})
    print()
    print('w =', w, ':', atoms)
    print()
    print('missing:', atoms - sources)


def test_v(n, levels):
    w = vzero(n)
    atoms = {v.reduce() for v in w.get_atoms()}

    sources = atoms.copy()
    gv = set()
    ge = set()
    while levels > 0:
        add = set()
        for sigma in SignedPermutation.all(n + 1):
            for j in [1]:
                top, bot, v, e = sigma.transition_graph(j, flip=True)
                if set(bot) and set(bot).issubset(sources):
                    print('* ', sigma, j, bot)
                    add |= set(top)
                    gv |= v
                    ge |= e
        levels -= 1
        sources |= add

    draw_graph(gv, ge, colors=colorfn)

    w0 = SignedPermutation.longest_element
    v0 = vzero

    print()
    print('p  =', v0(n - 1), ':', {v.reduce() for v in v0(n - 1).get_atoms()})
    print('q  =', w0(n - 1), ':', {v.reduce() for v in w0(n - 1).get_atoms()})
    print()
    print('w =', w, ':', atoms)
    print()
    print('missing:', atoms - sources)


def test_fpf(n, levels):
    wp = SignedPermutation.longest_element(n - 1)
    w = SignedPermutation.longest_element(n)
    atoms = {v.reduce() for v in w.get_fpf_atoms()}

    sources = atoms.copy()
    sinks = set()
    gv = set()
    ge = set()
    while levels > 0:
        add = set()
        for sigma in SignedPermutation.all(n + 1):
            for j in [1,2]:
                top, bot, v, e = sigma.transition_graph(j, flip=False)
                if set(bot) and set(bot).issubset(sources):
                    print('* ', sigma, j, bot)
                    add |= set(top)
                    gv |= v
                    ge |= e
                    sinks |= set(bot)
        levels -= 1
        sources |= add

    p = n if n % 2 == 0 else n - 1
    q = 2 * n - p
    shape = shape_CII(p, q)
    printer = lambda v: str(v) + '\n' + str(shape(v.inflate(n).inverse()))
    draw_graph(gv, ge, colors=colorfn, printer=printer)

    print()
    print('p  =', wp, ':', {v.reduce() for v in wp.get_fpf_atoms()})
    print()
    print('w =', w, ':', atoms)
    print()
    print('missing:', atoms - sinks)


def test_fpf_v(n, levels):
    w = vzero(n)
    atoms = {v.reduce() for v in w.get_fpf_atoms()}

    sources = atoms.copy()
    gv = set()
    ge = set()
    while levels > 0:
        add = set()
        for sigma in SignedPermutation.all(n + 1):
            for j in [1,2]:
                top, bot, v, e = sigma.transition_graph(j, flip=True)
                if set(bot) and set(bot).issubset(sources):
                    print('* ', sigma, j, bot)
                    add |= set(top)
                    gv |= v
                    ge |= e
        levels -= 1
        sources |= add

    draw_graph(gv, ge)

    print()
    print('p  =', vzero(n - 2), ':', {v.reduce() for v in vzero(n - 2).get_fpf_atoms()})
    print('q  =', vzero(n - 1), ':', {v.reduce() for v in vzero(n - 1).get_fpf_atoms()})
    print()
    print('w =', w, ':', atoms)


def test_d(n, levels):
    twisted = n % 2 != 0
    w = vzero_plus(n)
    atoms = {v.reduce() for v in w.get_atoms_d(twisted=twisted)}

    sources = atoms.copy()
    gv = set()
    ge = set()
    i = 0
    while levels > 0:
        add = set()
        for sigma in SignedPermutation.all(n + 1, dtype=True):
            for j in [1 + i]:
                top, bot, v, e = sigma.dtransition_graph(j, flip=True)
                if set(bot) and set(bot).issubset(sources):
                    print('* ', sigma, j, bot)
                    add |= set(top)
                    gv |= v
                    ge |= e
        levels -= 1
        sources |= add
        i += 1

    draw_graph(gv, ge)

    v0 = vzero_plus

    print()
    print('p  =', v0(n - 1), ':', {v.reduce() for v in v0(n - 1).get_atoms_d(twisted=not twisted)})
    print()
    print('w =', w, ':', atoms)
