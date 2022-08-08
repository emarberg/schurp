from permutations import Permutation
from qp_utils import beissinger_rsk, rsk


def config_tikz(config):
    w, i = config
    sigma = Permutation(*w)

    ans = []
    ans += ["\\arcstart{"]

    vertices = []
    letter = 'a'
    for j, a in enumerate(w):
        loc = j + 1
        curr = "*{" + (letter if loc not in [i - 1, i, i + 1] else "\\bullet") + "}"
        letter = chr(ord(letter) + 1) if loc not in [i - 1, i, i + 1] else letter
        if loc < a:
            dist = (a - loc) / 2
            curr += "\\arc{" + str(dist) + "}{" + (a - loc) * "r" + "}"
        vertices += [curr]
    ans += " & ".join(vertices)
    ans += ["""}\\arcstop"""]

    return "".join(ans)

def triple_configuration(w, i):
    base = sorted({i - 1, i, i + 1} | {w(i - 1), w(i), w(i + 1)})
    n = len(base)
    psi = {a: j + 1 for j, a in enumerate(base)}
    return tuple(psi[w(a)] for a in base), psi[i]


def look_for_conjugacy(i, u, v):
    s = Permutation.s_i(i - 1)
    t = Permutation.s_i(i)
    r = s * t * s

    sus = s * u * s
    tut = t * u * t
    tsust = t * s * u * s * t
    stuts = s * t * u * t * s

    print('D_%s' % i, u.cycle_repr(), '-->', v.cycle_repr())
    print()
    print(i, i + 1, i + 2, '-->')
    print(u(i), u(i + 1), u(i + 2))
    print()
    print('  ell(s u s) =', sus.length() - u.length(), '  ell(t u t) =', tut.length() - u.length(), '  ell(s t u t s) =', stuts.length() - tut.length())
    print()
    print('  ell(t u t) =', tut.length() - u.length(), '  ell(s u s) =', sus.length() - u.length(), '  ell(t s u s t) =', tsust.length() - sus.length())
    print()
    if u == v:
        print('  u == v')
    elif r * u * r == v:
        print('  r * u * r == v')
    elif s * u * s == v:
        print('  s * u * s == v')
        assert (tut.length() <= u.length() < sus.length() < tsust.length()) or (tut.length() > u.length() > sus.length() >= tsust.length())
    elif t * u * t == v:
        print('  t * u * t == v')
        assert (sus.length() <= u.length() < tut.length() < stuts.length()) or (sus.length() > u.length() > tut.length() >= stuts.length())

    else:
        print('  ? * u * ? == v')
        input('')
    print()


def test_beissinger_ops_alt(n=7):
    test_beissinger_ops(n, True)


def test_beissinger_ops(n=7, sgn=False):
    seen = {}
    for m in range(n + 1):
        for w in Permutation.involutions(n):
            w = [w(i + 1) for i in range(n)]
            btab = beissinger_rsk(w, sgn=sgn)
            seen[btab] = tuple(w)

    configs = {}
    for tab in seen:
        m = len(tab)
        for i in range(2, m):
            alt = tab.dual_equivalence_operator(i - 1)
            
            u = Permutation(*seen[tab])
            v = Permutation(*seen[alt])
            base = set(range(1, m + 1)) - {i - 1, i, i + 1, u(i - 1), u(i), u(i + 1)}
            assert all(u(j) == v(j) for j in base)

            tabconfig = triple_configuration(u, i)
            altconfig = triple_configuration(v, i)
            if tabconfig not in configs:
                configs[tabconfig] = altconfig
            assert configs[tabconfig] == altconfig

    for tabconfig, altconfig in configs.items():
        if tabconfig == altconfig:
            continue
        w, i = tabconfig
        v, j = altconfig
        if Permutation.s_i(i) * Permutation(*w) * Permutation.s_i(i) == Permutation(*v):
            continue
        print("\\[")
        print(config_tikz(tabconfig))
        print("\\qquad\\longrightarrow\\qquad")
        print(config_tikz(altconfig))
        print("\\]")
        a, b, c = w[i - 2], w[i - 1], w[i]
        # print(a, b, c)
        print()
        # assert not(a < b < c) and not(a > b > c)
    print()
    print("\\newpage\\newpage")
    print()
    for tabconfig, altconfig in configs.items():
        if tabconfig == altconfig:
            continue
        w, i = tabconfig
        v, j = altconfig
        if Permutation.s_i(i - 1) * Permutation(*w) * Permutation.s_i(i - 1) == Permutation(*v):
            continue
        print("\\[")
        print(config_tikz(tabconfig))
        print("\\qquad\\longrightarrow\\qquad")
        print(config_tikz(altconfig))
        print("\\]")
        a, b, c = w[i - 2], w[i - 1], w[i]
        # print(a, b, c)
        print()
        # assert not(a < b < c) and not(a > b > c)
    print()
    print("\\newpage\\newpage")
    print()
    for tabconfig, altconfig in configs.items():
        if tabconfig != altconfig:
            continue
        w, i = tabconfig
        print("\\[")
        print(config_tikz(tabconfig))
        print("\\qquad\\longrightarrow\\qquad")
        print(config_tikz(altconfig))
        print("\\]")
        a, b, c = w[i - 2], w[i - 1], w[i]
        # print(a, b, c)
        print()
        # assert a < b < c or a > b > c
    print(len(configs))
            

def test_all_row(n=7):
    seen = {}
    for w in Permutation.involutions(n):
        w = [w(i + 1) for i in range(n)]
        btab = beissinger_rsk(w, sgn=False)
        p, q = rsk(w)
        assert p == q
        rsktab = p
        # print(w)
        # print(btab)
        # print(rsktab)
        assert btab == rsktab 
        seen[btab] = tuple(w)

    for tab in seen:
        print(seen[tab])
        print(tab)
        for i in range(2, n):
            alt = tab.dual_equivalence_operator(i - 1)

            u = Permutation(*seen[tab])
            v = Permutation(*seen[alt])
            # look_for_conjugacy(i, u, v)
            # if alt != tab:
            #    print(alt)
            #    print()

            a, b, c = u(i - 1), u(i), u(i + 1)

            if u == v:
                assert a < b < c or a > b > c
            else:
                assert not (a < b < c) and not (a > b > c)
        print()

def test_all_col(n=7):
    seen = {}
    for w in Permutation.involutions(n):
        w = [w(i + 1) for i in range(n)]
        btab = beissinger_rsk(w, sgn=True)
        rsktab = rsk(w)[0]
        # print(w, '=', Permutation(*w).cycle_repr())
        # print(btab)
        # print(rsktab)
        # print()
        assert btab not in seen
        seen[btab] = tuple(w)

    for tab in seen:
        print(seen[tab])
        print(tab)
        for i in range(2, n):
            alt = tab.dual_equivalence_operator(i - 1)

            u = Permutation(*seen[tab])
            v = Permutation(*seen[alt])
            # look_for_conjugacy(i, u, v)
            # if alt != tab:
            #    print(alt)
            #    print()

            # a, b, c = u(i - 1), u(i), u(i + 1)

            # if u == v:
            #     assert a < b < c or a > b > c
            # else:
            #     assert not (a < b < c) and not (a > b > c)
        print()