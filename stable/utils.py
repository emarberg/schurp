from .symmetric import SymmetricPolynomial
from .permutations import Permutation
from .polynomials import beta, X # noqa
import itertools


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


def _kromatic_helper(num_variables, coloring, vertices, edges):
    if vertices:
        v = vertices[0]
        vertices = vertices[1:]
        subset = set(range(1, 1 + num_variables))
        for w in edges.get(v, []):
            subset -= coloring.get(w, set())
        for k in range(1, len(subset) + 1):
            for s in itertools.combinations(subset, k):
                coloring[v] = set(s)
                for ans in _kromatic_helper(num_variables, coloring, vertices, edges):
                    yield ans
                del coloring[v]
    else:
        ans = 1
        for _, subset in coloring.items():
            for i in subset:
                ans *= X(i)
        yield ans


def kromatic(num_variables, vertices, edges):
    e = {v: set() for v in vertices}
    for a, b in edges:
        e[a].add(b)
        e[b].add(a)
    ans = 0
    for a in _kromatic_helper(num_variables, {}, vertices, e):
        ans += a * beta**(a.total_degree() - len(vertices))
    return SymmetricPolynomial.from_polynomial(ans)


def grothendieck(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.stable_grothendieck(num_variables, degree_bound=degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_Q(num_variables, mu, nu=(), degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_q(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_P(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.symplectic_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_p(num_variables, mu, nu, degree_bound)

def grothendieck_S(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.signed_involution_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_s(num_variables, mu, nu, degree_bound)


def shifted_ribbon(alpha):
    nu = []
    mu = []
    for i in range(1, len(alpha)):
        nu.append(sum(alpha[i:]))
        mu.append(nu[-1] + alpha[i - 1])
    mu.append(alpha[-1])
    return tuple(mu), tuple(nu)


GE = SymmetricPolynomial.ktheoretic_e
GE_expansion = SymmetricPolynomial.ktheoretic_e_expansion

m = SymmetricPolynomial.monomial
p = SymmetricPolynomial.powersum
p_expansion = SymmetricPolynomial.powersum_expansion

s = SymmetricPolynomial.schur
P = SymmetricPolynomial.schur_p
Q = SymmetricPolynomial.schur_q
S = SymmetricPolynomial.schur_s

G = grothendieck
GP = grothendieck_P
GQ = SymmetricPolynomial.stable_grothendieck_q
GS = grothendieck_S

mn_G = SymmetricPolynomial.mn_stable_grothendieck
mn_GP = SymmetricPolynomial.mn_stable_grothendieck_p
mn_GQ = SymmetricPolynomial.mn_stable_grothendieck_q

G_doublebar = SymmetricPolynomial.stable_grothendieck_doublebar
GP_doublebar = SymmetricPolynomial.stable_grothendieck_p_doublebar
GQ_doublebar = SymmetricPolynomial.stable_grothendieck_q_doublebar
GS_doublebar = SymmetricPolynomial.stable_grothendieck_s_doublebar

g = SymmetricPolynomial.dual_stable_grothendieck
gp = SymmetricPolynomial.dual_stable_grothendieck_p
gq = SymmetricPolynomial.dual_stable_grothendieck_q
gs = SymmetricPolynomial.dual_stable_grothendieck_s

mp_g = SymmetricPolynomial.mp_dual_stable_grothendieck
mp_gp = SymmetricPolynomial.mp_dual_stable_grothendieck_p
mp_gq = SymmetricPolynomial.mp_dual_stable_grothendieck_q

j = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck
jp = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_p
jq = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_q
js = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_s

decomposition_jp = SymmetricPolynomial.slow_decomposition_jp

schur_expansion = SymmetricPolynomial.schur_expansion

G_expansion = SymmetricPolynomial.grothendieck_expansion

g_expansion = SymmetricPolynomial.dual_grothendieck_expansion
j_expansion = SymmetricPolynomial.j_expansion
gp_expansion = SymmetricPolynomial.gp_expansion
gq_expansion = SymmetricPolynomial.gq_expansion
gs_expansion = SymmetricPolynomial.gs_expansion

mp_g_expansion = SymmetricPolynomial.mp_dual_grothendieck_expansion
mp_gp_expansion = SymmetricPolynomial.mp_gp_expansion
mp_gq_expansion = SymmetricPolynomial.mp_gq_expansion

mn_G_expansion = SymmetricPolynomial.mn_grothendieck_expansion
mn_GP_expansion = SymmetricPolynomial.mn_GP_expansion
mn_GQ_expansion = SymmetricPolynomial.mn_GQ_expansion

gp_free_expansion = SymmetricPolynomial.gp_free_expansion

jp_expansion = SymmetricPolynomial.jp_expansion
jq_expansion = SymmetricPolynomial.jq_expansion
js_expansion = SymmetricPolynomial.js_expansion

GP_expansion = SymmetricPolynomial.GP_expansion
GQ_expansion = SymmetricPolynomial.GQ_expansion
GS_expansion = SymmetricPolynomial.GS_expansion

P_expansion = SymmetricPolynomial.P_expansion
Q_expansion = SymmetricPolynomial.Q_expansion
S_expansion = SymmetricPolynomial.S_expansion
