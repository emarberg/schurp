from schubert import FPFGrothendieck, InvGrothendieck, Grothendieck
from permutations import Permutation
import subprocess
import time

O_GROTHENDIECK_DEGREE = {}
SP_GROTHENDIECK_DEGREE = {}
BASE_DIRECTORY = '/Users/emarberg/Downloads/symplectic-cmr/'


def nice_str(u, n):
    return ' '.join([str(u(i)) for i in range(1, n + 1)])


def get_orthogonal_grothendieck_degree(w):
    if w not in O_GROTHENDIECK_DEGREE:
        t0 = time.time()
        print('... computing inv Grothendieck')
        f = InvGrothendieck.get(w)
        t1 = time.time()
        print('... done', t1 - t0, '\n')

        t0 = time.time()
        print('... computing Grothendieck terms')
        dec = Grothendieck.decompose(f)
        a = {u: Grothendieck.get(u).degree() for u in dec}
        t1 = time.time()
        print('... done', t1 - t0, '\n')        
        
        print('dec =', dec)
        print()

        deg = max(a.values())
        
        found = None
        for u in a:
            if a[u] == deg:
                assert found is None
                found = u
        O_GROTHENDIECK_DEGREE[w] = (found, deg)
    return O_GROTHENDIECK_DEGREE[w]


def test_inv(n=6):
    ans = []
    for case, w in enumerate(Permutation.involutions(n)):
        if w.is_vexillary():
            print('CASE', case + 1, ':', 'z =', w.cycle_repr())
            print()
            a, degw = get_orthogonal_grothendieck_degree(w)
            print('  degree =', degw, 'atom =', nice_str(a, n + 1))
            print()
            w.print_involution_rothe_diagram(sep='.')
            print()
            ans += [(degw, '$' + w.cycle_repr().replace(',', '\\,') +'$ & $' + nice_str(a, n) + '$ & $' + str(degw) + '$ \\\\ & & \\\\ ')]
        # v = w.standardize({i for i in range(1, n + 1)} - {1, w(1)})
        # b, degv = get_symplectic_grothendieck_degree(v)
        # print(w.cycle_repr(), degw, '-->', v.cycle_repr(), degv)
        # print('  ', nice_str(a, n), '  ', nice_str(b, n - 2), '==', nice_str(a.standardize({i for i in range(1, n + 1)} - {1, w(1)}), n - 2), b == a.standardize({i for i in range(1, n + 1)} - {1, w(1)}))
        # print()
        # ell = len(w.get_symplectic_hecke_atoms())
        # if ell > 3:
        #    draw(w)
        #    input(str(ell))
    ans = [b for (a, b) in sorted(ans)]
    print('\n'.join(ans))


def get_symplectic_grothendieck_degree(w):
    if w not in SP_GROTHENDIECK_DEGREE:
        t0 = time.time()
        print('... computing fpf Hecke atoms')
        w.get_symplectic_hecke_atoms()
        t1 = time.time()
        print('... done', t1 - t0, '\n')

        t0 = time.time()
        print('... computing Grothendieck terms')
        a = {u: Grothendieck.get(u).degree() for u in w.get_symplectic_hecke_atoms()}
        t1 = time.time()
        print('... done', t1 - t0, '\n')        
        
        deg = max(a.values())
        
        # t0 = time.time()
        # print('... computing fpf Grothendieck of w =', w.cycle_repr())
        # deg = FPFGrothendieck.get(w).degree()
        # t1 = time.time()
        # print('... done', t1 - t0, '\n')  
        #
        # f = 0
        # for u in a:
        #     f += Grothendieck.get(u)
        # assert f == FPFGrothendieck.get(w)

        found = None
        for u in a:
            if a[u] == deg:
                assert found is None
                found = u
        SP_GROTHENDIECK_DEGREE[w] = (found, deg)
    return SP_GROTHENDIECK_DEGREE[w]


def test_fpf(n=6):
    ans = []
    for case, w in enumerate(Permutation.fpf_involutions(n)):
        print('CASE', case + 1, ':', 'z =', w.cycle_repr())
        a, degw = get_symplectic_grothendieck_degree(w)
        print('  degree =', degw, 'atom =', nice_str(a, n))
        print()
        w.print_fpf_rothe_diagram(sep='.')
        ans += [(degw, '$' + w.cycle_repr().replace(',', '\\,') +'$ & $' + nice_str(a, n) + '$ & $' + str(degw) + '$ \\\\ & & \\\\ ')]
        # v = w.standardize({i for i in range(1, n + 1)} - {1, w(1)})
        # b, degv = get_symplectic_grothendieck_degree(v)
        # print(w.cycle_repr(), degw, '-->', v.cycle_repr(), degv)
        # print('  ', nice_str(a, n), '  ', nice_str(b, n - 2), '==', nice_str(a.standardize({i for i in range(1, n + 1)} - {1, w(1)}), n - 2), b == a.standardize({i for i in range(1, n + 1)} - {1, w(1)}))
        # print()
        # ell = len(w.get_symplectic_hecke_atoms())
        # if ell > 3:
        #    draw(w)
        #    input(str(ell))
    ans = [b for (a, b) in sorted(ans)]
    print('\n'.join(ans))


def draw(w):
    n = len(w.oneline)

    def printer(e):
        g = Grothendieck.get(Permutation(*e).inverse())
        return ' '.join([str(i) for i in e] + [str(i) for i in range(len(e) + 1, n + 1)]) + '  |  ' + str(g.degree())

    seen = set()
    edges = set()
    for e in w.get_symplectic_hecke_atoms():
        e = e.inverse().oneline
        while len(e) < n:
            e.append(len(e) + 1)
        e = tuple(e)
        seen.add(e)
        for i in range(0, n - 3, 2):
            a, d, b, c = e[i:i + 4]
            if a < b < c < d:
                f = list(e[:])
                f[i:i + 4] = [b, c, a, d]
                f = tuple(f)
                edges.add((e, f))

                f = list(f)
                f[i:i + 4] = [b, d, a, c]
                f = tuple(f)
                edges.add((e, f))
            b, c, a, d = e[i:i + 4]
            if a < b < c < d:
                f = list(e[:])
                f[i:i + 4] = [b, d, a, c]
                f = tuple(f)
                edges.add((e, f))

    s = ['digraph G {']
    s += ['    overlap=false;']
    s += ['    splines=true;']
    for x in seen:
        s += ['    "%s";' % printer(x)]
    for (e, f) in edges:
        s += ['    "%s" -> "%s";' % (printer(e), printer(f))]
    s += ['}']
    s = '\n'.join(s)

    filename = '_'.join([str(i) for i in w.oneline])
    dot_filename = BASE_DIRECTORY + 'graphs/' + '%s.dot' % filename
    png_filename = BASE_DIRECTORY + 'graphs/' + '%s.png' % filename
    with open(dot_filename, 'w') as f:
        f.write(s)
    subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
    subprocess.run(["open", png_filename])
