from .symmetric import SymmetricPolynomial
from .partitions import Partition
from .polynomials import beta, X, Y, Polynomial
from .vectors import Vector
import itertools
import math


Y_INDEX = Polynomial.MIN_INT


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


def _kromatic_helper(num_variables, coloring, vertices, edges, weights, oriented):
    if vertices:
        v = vertices[0]
        vertices = vertices[1:]
        subset = set(range(1, 1 + num_variables))
        for w in edges.get(v, []):
            subset -= coloring.get(w, set())
        for k in range(1, len(subset) + 1):
            for s in itertools.combinations(subset, k):
                coloring[v] = set(s)
                for ans in _kromatic_helper(num_variables, coloring, vertices, edges, weights, oriented):
                    yield ans
                del coloring[v]
    else:
        ans = 1
        for v, subset in coloring.items():
            for i in subset:
                ans *= (beta * X(i))**weights.get(v, 1)
            if oriented:
                for w in edges.get(v, []):
                    if v < w and max(subset) < max(coloring[w]):
                        ans *= Y()
                    # if v < w:
                    #    ans *= Y()**len([i for i in subset if i < max(coloring[w])])
                    # if v < w:
                    #     ans *= Y()**len([(i, j) for i in coloring[v] for j in coloring[w] if i < j])
        yield ans


def oriented_kromatic(num_variables, vertices, edges, weights=None):
    return kromatic(num_variables, vertices, edges, weights, True)


def kromatic(num_variables, vertices, edges, weights=None, oriented=False):
    weights = weights or {}
    total_weight = 0
    for v in vertices:
        total_weight += weights.get(v, 1)

    e = {v: set() for v in vertices}
    for a, b in edges:
        e[a].add(b)
        e[b].add(a)
    ans = 0

    for a in _kromatic_helper(num_variables, {}, vertices, e, weights, oriented):
        ans += a * beta**(-total_weight)
    # if oriented:
    #    ans = ans if type(ans) == int else ans.set_variable(0, 0)
    # if oriented:
    #    return ans
    return SymmetricPolynomial.from_polynomial(ans)


KMONOMIAL_CACHE = {}


def kmonomial(num_variables, mu):
    key = (num_variables, mu)
    cache = KMONOMIAL_CACHE
    if key not in cache:
        n = len(mu)
        weights = {i + 1: mu[i] for i in range(n)}
        v = list(range(1, n + 1))
        cache[key] = kromatic(num_variables, v, complete_graph(v), weights)
    return cache[key]


def monic_kmonomial(num_variables, mu):
    return kmonomial(num_variables, mu) // Partition.stabilizer_order(mu)


def kmonomial_expansion(f):
    exp = SymmetricPolynomial._expansion(f, monic_kmonomial, SymmetricPolynomial._get_term_from_lowest_degree)
    ans = Vector({mu: val // Partition.stabilizer_order(mu) for mu, val in exp.items()}, multiplier=exp.multiplier, sorter=exp.sorter)
    return ans.set_beta(1)
