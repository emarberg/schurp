from stable.partitions import Partition
from stable.polynomials import beta
from stable.vectors import Vector
from stable.utils import g, g_expansion, gp, gp_expansion


g_skew_cache = {'computed': set()}
gp_skew_cache = {'computed': set()}


def print_antipode_test(alpha):
    v=barShRibbon_expansion(antipode_in_ribbon_basis(barShRibbon(alpha)))
    for a in v:
        print(v[a])
        Partition.print_shifted(*Partition.shifted_ribbon(a))
        print()
    print()
    print()
    gamma = list(reversed(alpha))
    gamma[-1] -= 1
    gamma[0] += 1
    gamma = tuple(gamma)
    Partition.print_shifted(*Partition.shifted_ribbon(gamma))


def transpose(alpha):
    n = sum(alpha)
    alpha = reversed(alpha)
    des = []
    for a in alpha:
        des.append(a + (des[-1] if des else 0))
    des = [0] + sorted(set(range(1, n)) - set(des))
    ans = []
    for i in range(1, len(des)):
        ans.append(des[i] - des[i - 1])
    if n > 0:
        ans.append(n - des[-1])
    return tuple(ans)


def generate_g_skew_cache(nvars, nboxes):
    if (nvars, nboxes) in g_skew_cache['computed']:
        return
    for lam in Partition.all(nboxes):
        for mu in Partition.subpartitions(lam):
            f = g_expansion(g(nvars, lam, mu))
            if f not in g_skew_cache:
                g_skew_cache[f] = []
            g_skew_cache[f].append((nvars, lam, mu))
    g_skew_cache['computed'].add((nvars, nboxes))


def find_g_skew(nvars, f, limit=8):
    seen = set()
    for nboxes in range(limit + 1):
        print('nboxes =', nboxes)
        generate_g_skew_cache(nvars, nboxes)
        for _, lam, nu in g_skew_cache.get(f, []):
            if (lam, nu) not in seen:
                yield lam, nu
                seen.add((lam, nu))


def generate_gp_skew_cache(nvars, nboxes):
    if (nvars, nboxes) in gp_skew_cache['computed']:
        return
    for lam in Partition.all(nboxes, strict=True):
        for mu in Partition.subpartitions(lam, strict=True):
            f = gp_expansion(gp(nvars, lam, mu))
            if f not in gp_skew_cache:
                gp_skew_cache[f] = []
            gp_skew_cache[f].append((nvars, lam, mu))
    gp_skew_cache['computed'].add((nvars, nboxes))


def find_gp_skew(nvars, f, limit=8):
    seen = set()
    for nboxes in range(limit + 1):
        print('nboxes =', nboxes)
        generate_gp_skew_cache(nvars, nboxes)
        for _, lam, nu in gp_skew_cache.get(f, []):
            if (lam, nu) not in seen:
                yield lam, nu
                seen.add((lam, nu))


def is_peak_composition(alpha):
    return all(a > 0 for a in alpha) and all(a >= 2 for a in alpha[:-1])


def I(alpha):
    des = [0]
    for a in alpha[:-1]:
        des.append(des[-1] + a)
    des = des[1:]
    return des


def Lambda(alpha):
    if sum(alpha) == 0:
        return ()
    des = I(alpha)
    des = [0] + [i for i in des if 0 < i - 1 and i - 1 not in des] + [sum(alpha)]
    return tuple(des[i] - des[i - 1] for i in range(1, len(des)))


def ribbon_multiplier_homogeneous(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])
        b = tuple(b)

        yield a, 1
        yield b, 1


def ribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])
        b = tuple(b)

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])
        c = tuple(c)

        yield a, 1
        yield b, 1
        yield c, beta


def shribbon_multiplier_homogeneous(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g

        if is_peak_composition(a):
            yield tuple(a), 1
        if is_peak_composition(b):
            yield tuple(b), 1
        if is_peak_composition(d):
            yield tuple(d), 1


def shribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g

        e = list(alpha)
        e += g

        if is_peak_composition(a):
            yield tuple(a), 1
        if is_peak_composition(b):
            yield tuple(b), 1
        if is_peak_composition(c):
            yield tuple(c), beta
        if is_peak_composition(d):
            yield tuple(d), 1
        if is_peak_composition(e):
            yield tuple(e), beta


def bar_shribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g

        e = list(alpha)
        e += g

        f = list(alpha)
        f[-1] += gamma[0] - 2
        f += list(gamma[1:])

        p = int(alpha[-1] > 1)
        q = int(gamma[0] > 2 or gamma == (2,))

        if is_peak_composition(a):
            yield tuple(a), 1
        if is_peak_composition(b):
            yield tuple(b), 2
        if is_peak_composition(c):
            yield tuple(c), (1 + p + q) * beta
        if is_peak_composition(d):
            yield tuple(d), 1
        if is_peak_composition(e):
            yield tuple(e), beta
        if is_peak_composition(f):
            yield tuple(f), p * q * beta**2


def Ribbon(alpha, multiplier=ribbon_multiplier):
    return Vector({alpha: 1}, multiplier=multiplier)


def barShRibbon(alpha):
    ans = 0
    n = len(alpha)
    for v in range(2**n):
        delta = list(alpha)
        for i in range(n):
            delta[i] -= v % 2
            v = v // 2
        if any(d < 2 for d in delta[:-1]) or (n > 0 and delta[-1] < 1):
            continue
        ans += ShRibbon(tuple(delta)) * 2**(n + sum(delta) - sum(alpha)) * beta**(sum(alpha) - sum(delta))
    return ans


def ShRibbon(alpha):
    ans = 0
    n = sum(alpha)
    if any(d < 2 for d in alpha[:-1]) or (n > 0 and alpha[-1] < 1):
        raise Exception()
    for gamma in Partition.compositions(n):
        if Lambda(gamma) == alpha:
            ans += Ribbon(gamma)
    return ans


def ShRibbon_homogeneous(alpha):
    ans = 0
    n = sum(alpha)
    for gamma in Partition.compositions(n):
        if Lambda(gamma) == alpha:
            ans += Ribbon(gamma, ribbon_multiplier_homogeneous)
    return ans


def barShRibbon_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        two = 2**len(alpha)
        coeff = f[alpha]
        assert coeff % two == 0
        ans += Ribbon(alpha) * (coeff // two)
        f -= barShRibbon(alpha) * (coeff // two)
    return ans


def ShRibbon_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        coeff = f[alpha]
        ans += Ribbon(alpha) * coeff
        f -= ShRibbon(alpha) * coeff
    return ans


SIMPLE_ANTIPODE_CACHE = {}


def simple_antipode_in_ribbon_basis(n):
    if n in SIMPLE_ANTIPODE_CACHE:
        return SIMPLE_ANTIPODE_CACHE[n]
    
    ans = 0
    g = Ribbon_free_expansion(Ribbon(n * (1,)))
    for alpha, coeff in g.items():
        if type(coeff) != int and coeff.substitute(0, 0) == 0:
            continue
        term = 1
        for a in alpha:
            term *= Ribbon((a,))
        ans += term * coeff
    ans = ans * (-1)**n

    SIMPLE_ANTIPODE_CACHE[n] = ans
    return ans


def antipode_in_ribbon_basis(f):
    g = Ribbon_free_expansion(f)
    ans = 0
    for alpha, coeff in g.items():
        term = 1
        for a in reversed(alpha):
            term *= simple_antipode_in_ribbon_basis(a)
        ans += term * coeff
    return ans


def test_antipode(n=5):
    for p in range(1, n + 1):
        for alpha in Partition.peak_compositions(p):
            pass


def Ribbon_free_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        term = 1
        for a in alpha:
            term *= Ribbon((a,))
        coeff = f[alpha]
        ans += Ribbon(alpha) * coeff
        f -= term * coeff
    return ans


def map_g(n, vec):
    bns = 0
    for alpha, coeff in vec.items():
        term = 1
        for a in alpha:
            term *= g(n, (a,))
        bns += term * coeff
    return g_expansion(bns)   


def ShRibbon_free_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        coeff = f[alpha]
        term = 1
        
        gamma = []
        for a in alpha:
            if a % 2 == 0:
                gamma += [1, a - 1]
            else:
                gamma += [a]
        alpha = tuple(gamma)

        for a in alpha:
            term *= ShRibbon((a,))

        ans += Ribbon(alpha) * coeff
        f -= term * coeff
    return ans


def map_gp(n, vec):
    bns = 0
    for alpha, coeff in vec.items():
        term = 1
        for a in alpha:
            term *= gp(n, (a,))
        bns += term * coeff
    return gp_expansion(bns)   


def barShRibbon_free_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        term = 1
        for a in alpha:
            term *= barShRibbon((a,))
        
        two = 2**len(alpha)
        coeff = f[alpha]
        assert coeff % two == 0
        coeff = coeff // two
        
        ans += Ribbon(alpha) * coeff
        f -= term * coeff
    return ans


def test_shribbon_product(n=5):
    for p in range(1, n + 1):
        for q in range(1, n + 1):
            for alpha in Partition.peak_compositions(p):
                for gamma in Partition.peak_compositions(q):
                    expected = Ribbon(alpha, shribbon_multiplier) * Ribbon(gamma, shribbon_multiplier)
                    computed = ShRibbon_expansion(ShRibbon(alpha) * ShRibbon(gamma))
                    if expected != computed:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('computed =', computed)
                        input('\n')
                    else:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('.')
                        print('\n')
                    # assert expected == computed


def test_bar_shribbon_product(n=5):
    for p in range(1, n + 1):
        for q in range(1, n + 1):
            for alpha in Partition.peak_compositions(p):
                for gamma in Partition.peak_compositions(q):
                    expected = Ribbon(alpha, bar_shribbon_multiplier) * Ribbon(gamma, bar_shribbon_multiplier)
                    computed = barShRibbon_expansion(barShRibbon(alpha) * barShRibbon(gamma))
                    if expected != computed:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('computed =', computed)
                        input('\n')
                    else:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('.')
                        print('\n')


def omega(n):
    return ShRibbon((n,)) if n != 2 else ShRibbon((2,)) + beta * ShRibbon((1,))


def test_omega(n=4):
    def r(k):
        return ShRibbon((k,))

    for k in range(1, n + 1):
        lhs = r(1) * r(k)
        rhs = 0 if k % 2 == 0 else r(k + 1)
        for i in range(k - 1):
            rhs += (-1) ** i * r(i + 2) * r(k - 1 - i)
        if lhs != rhs:
            print('R_1 * R_' + str(k))
            print()
            print(lhs, '==', rhs)
        assert lhs == rhs

