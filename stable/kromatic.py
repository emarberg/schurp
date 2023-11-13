from .symmetric import SymmetricPolynomial
from .partitions import Partition
from .polynomials import beta, X, Y, Polynomial
from .vectors import Vector
from .tableaux import Tableau
import itertools
import math


Y_INDEX = Polynomial.MIN_INT


def is_connected_graph(v, e):
    if len(v) == 0:
        return True
    component = {v[0]}
    size = len(component)
    while True:
        component |= {b for (a, b) in e if a in component} | {a for (a, b) in e if b in component}
        if len(component) == size:
            break
        size = len(component)
    return len(component) == len(v)


def is_cluster_graph(e):
    e = {tuple(sorted(edge)) for edge in e}
    for (i, j) in e:
        for (k, l) in e:
            v = {i, j, k, l}
            if len(v) == 3:
                a, b, c = sorted(v)
                if (a, b) not in e or (a, c) not in e or (b, c) not in e:
                    return False
    return True


SMALL_MULTIPERMUTATIONS_CACHE = {}


def is_small_multipermutation(n):
    def ans(w):
        if sorted(set(w)) != list(range(1, n + 1)):
            return False
        for i in range(len(w) - 1):
            if w[i] == w[i + 1]:
                return False
        return True
    return ans


def words(n, m):
    for v in range(n**m):
        ans = []
        for _ in range(m):
            ans.append((v % n) + 1)
            v = v // n
        yield tuple(ans)


def small_multipermutations(n, m):
    assert n >= 0 and m >= n
    if (n, m) not in SMALL_MULTIPERMUTATIONS_CACHE:
        ans = list(filter(is_small_multipermutation(n), words(n, m)))
        SMALL_MULTIPERMUTATIONS_CACHE[n, m] = ans
    for w in SMALL_MULTIPERMUTATIONS_CACHE[n, m]:
        yield w


def posets(n):
    edges = sorted([(i, j) for i in range(1, n + 1) for j in range(1, n + 1) if i != j])
    for k in range(1 + (n * (n - 1)) // 2):
        for subset in itertools.combinations(edges, k):
            if all((i, l) in subset for (i, j) in subset for (k, l) in subset if j == k):
                yield subset


def incomparability_graph(n, poset):
    v = list(range(1, n + 1))
    e = [(i, j) for i in v for j in v if i < j and (i, j) not in poset and (j, i) not in poset]
    return (v, e)


q_multinomial_cache = {}


def q_multinomial(n, *args):
    if len(args) == 0:
        return 0 if n != 0 else 1

    args = tuple(args) + ((n - sum(args),) if sum(args) < n else ())
    
    if n < 0 or min(args) < 0:
        return 0
    if n == 0:
        return 1
    
    key = (n,) + args
    cache = q_multinomial_cache
    if key not in cache:
        ans = 0
        m = len(args)
        for k in range(1, m + 1):
            for J in itertools.combinations(list(range(m)), k):
                newargs = list(args)
                for j in J:
                    newargs[j] -= 1
                term = (-1)**(k - 1) * q_multinomial(n - k, *newargs)
                for i in range(n - k + 1, n):
                    term *= (1 - Y()**i)
                ans += term
        cache[key] = ans
    return cache[key]


def weak_compositions(r, n, strict=False):
    assert r >= 0
    if r > 0 >= n:
        return
    elif r == n == 0:
        yield ()
    else:
        for i in range(1 if strict else 0, r + 1):
            for alpha in weak_compositions(r - i, n - 1, strict):
                yield (i,) + alpha


def strict_compositions(r, n):
    for alpha in weak_compositions(r, n, True):
        yield alpha


galois_number_cache = {}


def galois_number(r, n):
    if (r, n, False) not in galois_number_cache:
        ans = 0
        if r >= 0:
            for alpha in weak_compositions(r, n):
                ans += q_multinomial(r, *alpha)
        galois_number_cache[r, n, False] = ans
    return galois_number_cache[r, n, False]


def strict_galois_number(r, n):
    if (r, n, True) not in galois_number_cache:
        ans = 0
        if r >= 0:
            for alpha in strict_compositions(r, n):
                ans += q_multinomial(r, *alpha)
        galois_number_cache[r, n, True] = ans
    return galois_number_cache[r, n, True]


def y_terms(f):
    if type(f) == int:
        f = f * Polynomial.one()
    if type(f) == SymmetricPolynomial:
        max_y_degree = max([0] + [d.get(Y_INDEX, 0) for _, c in f.items() for d in ([] if type(c) == int else c.coeffs)])
    else:
        max_y_degree = max([0] + [d.get(Y_INDEX, 0) for d in f.coeffs])
    terms = []
    for _ in range(max_y_degree + 1):
        terms.append(f.set_variable(Y_INDEX, 0))
        f = (f - terms[-1]) * Y()**-1
    return terms


def oriented_expander(expander, f):
    terms = y_terms(f)
    ans = 0
    for i, t in enumerate(terms):
        ans += expander(t) * Y()**i
    return ans


def from_kmonomial_expansion(n, a):
    return sum(kmonomial(n, mu) * coeff for (mu, coeff) in a.items())


def from_p_expansion(n, a):
    fac = math.lcm(*[t.denominator for t in a.dictionary.values()])
    b = a * fac
    assert all(coeff.denominator == 1 for (mu, coeff) in b.items())
    ans = sum(p(n, mu) * coeff.numerator for (mu, coeff) in b.items())
    assert all(coeff % fac == 0 for (mu, coeff) in ans.items())
    return (ans // fac).polynomial()


def tex_partition(mu):
    if sum(mu) == 0:
        return '\\emptyset'
    ans = []
    for i, m in enumerate(mu):
        if i == 0 or m != mu[i - 1]:
            ans.append([m, 1])
        else:
            ans[-1][1] += 1
    ans = [str(a) if b == 1 else str(a) + '^{' + str(b) + '}' for (a, b) in ans]
    return '(' + ','.join(ans) + ')'


def tex_expansion(vec, basis_symbol, degree_bound=None):
    def get_key(k):
        return (sum(k), k, basis_symbol + '_{' + tex_partition(k) +'}')

    newvec = Vector(
        {get_key(k): v for (k, v) in vec.items() if degree_bound is None or sum(k) <= degree_bound},
        printer=lambda k: k[2],
        sorter=lambda k: k[:2]
    )
    ans = str(newvec).replace('*', ' ').replace('β', '\\beta')
    if len(newvec) < len(vec):
        ans += ' + \\dots'
    return '$ ' + ans + ' $'


def collect_by_numbers_of_parts(vec):
    ans = []
    for mu, c in vec.items():
        while len(mu) > len(ans):
            ans.append(0)
        ans[len(mu) - 1] += c
    return ans


def circle_graph(vertices):
    return {(vertices[i], vertices[i + 1]) for i in range(len(vertices) - 1)} | {(vertices[-1], vertices[0])}


def path_graph(vertices):
    return {(vertices[i], vertices[i + 1]) for i in range(len(vertices) - 1)}


def complete_graph(vertices):
    return {(a, b) for a in vertices for b in vertices if a != b}


def _kromatic_helper(num_variables, coloring, vertices, edges, weights, oriented, chromatic):
    if vertices:
        v = vertices[0]
        vertices = vertices[1:]
        subset = set(range(1, 1 + num_variables))
        for w in edges.get(v, []):
            subset -= coloring.get(w, set())
        for k in range(1, 1 + (1 if chromatic else len(subset))):
        # for k in range(1, 1 + (1 if chromatic else 2)):
            for s in itertools.combinations(subset, k):
                coloring[v] = set(s)
                for ans in _kromatic_helper(num_variables, coloring, vertices, edges, weights, oriented, chromatic): # or k == 2):
                    yield ans
                del coloring[v]
    else:
        ans = 1
        for v, subset in coloring.items():
            for i in subset:
                ans *= (X(i) if chromatic else (beta * X(i)))**weights.get(v, 1)
            if oriented == 1:
                for w in edges.get(v, []):
                    if v < w:
                        ans *= Y()**len([(i, j) for i in coloring[v] for j in coloring[w] if i < j])
            elif oriented == 2:
                for w in edges.get(v, []):
                    if v < w and max(subset) < max(coloring[w]):
                        ans *= Y()
            elif oriented == 3:
                for w in edges.get(v, []):
                    if v < w and min(subset) < min(coloring[w]):
                        ans *= Y()
        yield ans


def chromatic(num_variables, vertices, edges, weights=None, oriented=False):
    return kromatic(num_variables, vertices, edges, weights, False, True)


def oriented_chromatic(num_variables, vertices, edges, weights=None, oriented=False):
    return kromatic(num_variables, vertices, edges, weights, True, True)


def oriented_chromatic_polynomial(num_variables, vertices, edges, weights=None, oriented=False):
    return kromatic(num_variables, vertices, edges, weights, True, True, symmetrize=False)


def oriented_kromatic(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 1)


def oriented_kromatic_polynomial(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 1, symmetrize=False)


def max_oriented_kromatic(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 2)


def max_oriented_kromatic_polynomial(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 2, symmetrize=False)


def min_oriented_kromatic(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 3)


def min_oriented_kromatic_polynomial(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, 3, symmetrize=False)


def kromatic(num_variables, vertices, edges, weights=None, oriented=0, chromatic=False, symmetrize=True):
    weights = weights or {}
    total_weight = 0
    for v in vertices:
        total_weight += weights.get(v, 1)

    e = {v: set() for v in vertices}
    for a, b in edges:
        e[a].add(b)
        e[b].add(a)
    ans = 0

    for a in _kromatic_helper(num_variables, {}, vertices, e, weights, oriented, chromatic):
        ans += a * (1 if chromatic else beta**(-total_weight))
    return SymmetricPolynomial.from_polynomial(ans) if symmetrize else ans


KMONOMIAL_CACHE = {}
KPOWERSUM_CACHE = {}


def kmonomial(num_variables, mu):
    key = (num_variables, mu)
    cache = KMONOMIAL_CACHE
    if key not in cache:
        n = len(mu)
        weights = {i + 1: mu[i] for i in range(n)}
        v = list(range(1, n + 1))
        cache[key] = kromatic(num_variables, v, complete_graph(v), weights)
    return cache[key]


def kpowersum(num_variables, mu):
    key = (num_variables, mu)
    cache = KPOWERSUM_CACHE
    if key not in cache:
        n = len(mu)
        weights = {i + 1: mu[i] for i in range(n)}
        v = list(range(1, n + 1))
        cache[key] = kromatic(num_variables, v, [], weights)
    return cache[key]


def monic_kmonomial(num_variables, mu):
    return kmonomial(num_variables, mu) // Partition.stabilizer_order(mu)


def kmonomial_expansion(f):
    exp = SymmetricPolynomial._expansion(f, monic_kmonomial, SymmetricPolynomial._get_term_from_lowest_degree)
    ans = Vector({mu: val // Partition.stabilizer_order(mu) for mu, val in exp.items()}, multiplier=exp.multiplier, sorter=exp.sorter)
    return ans.set_beta(1)


def kpowersum_expansion(f, degree_bound):
    if not f:
        return Vector(sorter=lambda x: (sum(x), x))
    
    f = f.set_beta(1)

    m = sorted(f, key=lambda x: (x.degree(), x.index()))[0]
    nvars = m.n
    mu = m.mu

    if sum(mu) > degree_bound:
        return Vector(sorter=lambda x: (sum(x), x))

    coeff = f[m]
    denom = Partition.stabilizer_order(mu)
    assert coeff % denom == 0
    coeff = coeff // denom
    
    exp = kpowersum_expansion(f - kpowersum(nvars, mu).set_beta(1) * coeff, degree_bound) 
    return Vector({mu: coeff}, sorter=exp.sorter) + exp
