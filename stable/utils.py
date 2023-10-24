from .symmetric import SymmetricPolynomial
from .permutations import Permutation
from .polynomials import beta, X, Y, Polynomial # noqa
from .partitions import Partition
from .vectors import Vector
from fractions import Fraction


def shortest_expansion(f):
    expanders = [
        e_expansion,
        GE_expansion,
        p_expansion,
        schur_expansion,
        G_expansion,
        g_expansion,
        j_expansion,
        gp_expansion,
        gq_expansion,
        gs_expansion,
        mp_g_expansion,
        mp_gp_expansion,
        mp_gq_expansion,
        mn_G_expansion,
        mn_GP_expansion,
        mn_GQ_expansion,
        gp_free_expansion,
        jp_expansion,
        jq_expansion,
        js_expansion,
        GP_expansion,
        GQ_expansion,
        GS_expansion,
        P_expansion,
        Q_expansion,
        S_expansion,
    ]
    best = None
    ans = []
    for exp in expanders:
        try:
            res = exp(f)
            ell = len(res)
            if best is None or ell < best:
                best = ell
                ans = [(res, exp.__name__)]
            elif best == ell:
                ans += [(res, exp.__name__)]
        except:
            continue
    return ans


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


def pprod(a, b):
    return [(Partition.sort(a + b, trim=True), 1)]


e_to_p_cache = {0: Vector({(): 1}, multiplier=pprod)}


def e_to_p(k):
    assert k >= 0

    def combine(v, i, k):
        return Vector({Partition.sort(key + (i,), trim=True): coeff * (-1)**(i - 1) * Fraction(1, k) for (key, coeff) in v.items()})

    if k not in e_to_p_cache:
        ans = Vector()
        for i in range(1, k + 1):
            ans += combine(e_to_p(k - i), i, k)
        ans.multiplier = pprod
        ans.sorter = lambda tup: (sum(tup), tup)
        e_to_p_cache[k] = ans

    return e_to_p_cache[k]


def p_expansion(f):
    exp = e_expansion(f)
    x = {key: val.set(0, 1) for key, val in exp.items()}
    assert all(val.is_integer() for _, val in x.items())
    exp = Vector({key: val.constant_term() for key, val in x.items()}, multiplier=exp.multiplier, sorter=exp.sorter)

    ans = 0
    for mu, coeff in exp.items():
        term = e_to_p(0) * coeff
        for part in mu:
            term *= e_to_p(part)
        ans += term
    return ans


e = SymmetricPolynomial.e
e_expansion = SymmetricPolynomial.e_expansion

GE = SymmetricPolynomial.ktheoretic_e
GE_expansion = SymmetricPolynomial.ktheoretic_e_expansion

m = SymmetricPolynomial.monomial
p = SymmetricPolynomial.powersum

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
