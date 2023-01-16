from partitions import Partition
from permutations import Permutation
from qp_utils import beissinger_rsk, rsk
from collections import defaultdict


def config_tikz(config):
    w, i = config
    sigma = Permutation(*w)

    ans = []
    ans += ["\\xy\\xymatrix@R=.1cm@C=.5cm{"]

    vertices = []
    letter = 'a'
    for j, a in enumerate(w):
        loc = j + 1
        curr = "*{" + (letter if loc not in [i - 1, i, i + 1] else "\\bullet") + "}"
        letter = chr(ord(letter) + 1) if loc not in [i - 1, i, i + 1] else letter
        if loc < a:
            dist = (a - loc) * 0.3
            curr += "\\arc{" + str(dist) + "}{" + (a - loc) * "r" + "}"
        vertices += [curr]
    ans += " & ".join(vertices)
    ans += ["""}\\endxy"""]

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

    round_a, round_b, round_c, round_d = {}, {}, {}, {}
    for tabconfig, altconfig in configs.items():
        if tabconfig == altconfig:
            round_a[tabconfig] = altconfig
            continue

        w, i = tabconfig
        v, _ = altconfig

        w = Permutation(*w)
        v = Permutation(*v)
        s = Permutation.s_i(i - 1)
        t = Permutation.s_i(i)
        r = s * t * s

        if v == s * w * s:
            round_b[tabconfig] = altconfig
            assert sgn or min(v(i - 1), v(i)) < v(i + 1) < max(v(i), v(i - 1))
        elif v == t * w * t:
            round_c[tabconfig] = altconfig
            assert sgn or min(v(i), v(i + 1)) < v(i - 1) < max(v(i), v(i + 1))
        elif v == r * w * r:
            round_d[tabconfig] = altconfig
            assert {v(i - 1), v(i), v(i + 1)} == {i - 1, i, i + 1}
        else:
            raise Exception

    for rnd, tag in [(round_a, "y"), (round_b, "s_{i-1}ys_{i-1}"), (round_c, "s_iys_i"), (round_d, "(i-1,i+1)y(i-1,i+1)")]:
        seen = set()
        print("\\[y \\qquad\\sim\\qquad z=%s \\]" % tag)
        for tabconfig, altconfig in rnd.items():
            if tabconfig not in seen:
                seen.add(tabconfig)
                seen.add(altconfig)
            else:
                continue
            print()
            print("\\smallskip")
            print()
            print("\\[")
            print(config_tikz(tabconfig))
            print("\\qquad\\sim\\qquad")
            print(config_tikz(altconfig))
            print("\\]")

            w, i = tabconfig
            w = Permutation(*w)
            v, _ = altconfig
            v = Permutation(*v)

            def is_des(u, a):
                x = -a if u(a) == a else a if u(a) in [i - 1, i, i + 1] else u(a)
                y = -a - 1 if u(a + 1) == a + 1 else a + 1 if u(a + 1) in [i - 1, i, i + 1] else u(a + 1)
                return x > y

            print("\\[", "\\text{" + ("des" if is_des(w, i - 1) else "asc") + "}", ",", "\\text{" + ("des" if is_des(w, i) else "asc") + "}")
            print("\\qquad\\sim\\qquad")
            print("\\text{" + ("des" if is_des(v, i - 1) else "asc") + "}", ",", "\\text{" + ("des" if is_des(v, i) else "asc") + "}", "\\]")

        print()
        print("\\newpage")
        print()
    print(list(map(len, [round_a, round_b, round_c, round_d])))


def test_all_row(n=7):
    seen = {}
    for w in Permutation.involutions(n):
        f = len([i for i in range(n) if i + 1 == w(i + 1)])
        
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

        oddcols = sum([i % 2 for i in btab.partition().transpose()])
        assert f == oddcols

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


SYT = [1, 1, 2, 4, 10, 26, 76, 232, 764, 2620, 9496, 35696, 140152, 568504, 2390480, 10349536, 46206736, 211799312, 997313824, 4809701440, 23758664096, 119952692896, 618884638912, 3257843882624, 17492190577600, 95680443760576, 532985208200576, 3020676745975552]


def test_all_col(n=7):
    seen = {}
    for count, w in enumerate(Permutation.involutions(n)):
        if count % (SYT[n] // 100 + 1) == 0:
            print('left:', int((SYT[n] - count) / SYT[n] * 100), '%')
        f = len([i for i in range(n) if i + 1 == w(i + 1)])

        w = [w(i + 1) for i in range(n)]
        btab = beissinger_rsk(w, sgn=True)
        rsktab = rsk(w)[0]
        # print(w, '=', Permutation(*w).cycle_repr())
        # print(btab)
        # print(rsktab)
        # print()
        assert btab not in seen
        seen[btab] = tuple(w)

        oddrows = sum([i % 2 for i in btab.partition()])
        assert f == oddrows

    print()
    mapping = {}
    for tab, w in seen.items():
        v = seen[beissinger_rsk(w, sgn=False).transpose()]
        if v == w:
            print('S_%s' % n, ' case, Psi fixed point:', w)
        mapping[w] = v

    print()
    print('done first part')
    print()

    orders = defaultdict(int)
    for tab, w in seen.items():
        v = w
        o = 1
        while True:
            v = seen[beissinger_rsk(v, sgn=False).transpose()]
            if v == w:
                break
            o += 1
        orders[o] += 1
    cycletype = ""
    mu = []
    for o, v in sorted(orders.items()):
        assert v % o == 0
        for _ in range(v // o):
            mu += [o]
        cycletype += str(o) + "^{" + str(v // o) + "} "
    print(n, ':', cycletype)
    print(mu)
    mu = Partition(*reversed(mu))
    print(mu.shape)

    for tab in seen:
        # print(seen[tab], '=', Permutation(*seen[tab]).cycle_repr())
        # print(tab)
        # print()
        for i in range(2, n):
            alt = tab.dual_equivalence_operator(i - 1)

            u = Permutation(*seen[tab])
            v = Permutation(*seen[alt])

            # print("D_%s" % i)
            # print()
            # print(seen[alt], '=', v.cycle_repr())
            # print(alt)
            # print()
            
            s = Permutation.s_i(i - 1)
            t = Permutation.s_i(i)

            def tilde(u, j):
                if u(j) == j:
                    return -j
                if u(j) in [i - 1, i, i + 1]:
                    return j
                return u(j)

            a, b, c = tilde(u, i - 1), tilde(u, i), tilde(u, i + 1)

            if u == v:
                assert a < b < c or a > b > c
            elif s * u * s == v:
                assert a < c < b or b < c < a
            elif t * u * t == v:
                assert b < a < c or c < a < b
            else:
                assert False

    return mapping
