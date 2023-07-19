from schubert import FPFGrothendieck, Grothendieck
from permutations import Permutation
import subprocess

SP_GROTHENDIECK_DEGREE = {}
BASE_DIRECTORY = '/Users/emarberg/Downloads/symplectic-cmr/'


def nice_str(u, n):
    return ' '.join([str(u.inverse()(i)) for i in range(1, n + 1)])


def get_symplectic_grothendieck_degree(w):
    if w not in SP_GROTHENDIECK_DEGREE:
        a = {u: Grothendieck.get(u).degree() for u in w.get_symplectic_hecke_atoms()}
        deg = FPFGrothendieck.get(w).degree()
        found = None
        for u in a:
            if a[u] == deg:
                assert found is None
                found = u
        SP_GROTHENDIECK_DEGREE[w] = (found, deg)
    return SP_GROTHENDIECK_DEGREE[w]


def test_fpf(n=6):
    for w in Permutation.fpf_involutions(n):    
        a, degw = get_symplectic_grothendieck_degree(w)
        print(w.cycle_repr(), degw)
        print('  ', nice_str(a, n))
        
        # v = w.standardize({i for i in range(1, n + 1)} - {1, w(1)})
        # b, degv = get_symplectic_grothendieck_degree(v)
        # print(w.cycle_repr(), degw, '-->', v.cycle_repr(), degv)
        # print('  ', nice_str(a, n), '  ', nice_str(b, n - 2), '==', nice_str(a.standardize({i for i in range(1, n + 1)} - {1, w(1)}), n - 2), b == a.standardize({i for i in range(1, n + 1)} - {1, w(1)}))
        print()
        #ell = len(w.get_symplectic_hecke_atoms())
        #if ell > 3:
        #    draw(w)
        #    input(str(ell))


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
