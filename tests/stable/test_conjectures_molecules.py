from stable.permutations import Permutation
from stable.insertion import InsertionAlgorithm
from stable.tableaux import Tableau
from stable.partitions import Partition
from collections import defaultdict


MAPPING_CACHE = {}


def irsk(pi, n=None):
    if n is None:
        n = pi.rank
    cycles = sorted([(pi(i), i) for i in range(1, n + 1) if i <= pi(i)])
    tab = Tableau()
    for b, a in cycles:
        if a == b:
            tab = tab.add(1, tab.max_column() + 1, a)
        else:
            p, q = InsertionAlgorithm.hecke(tab.row_reading_word() + (a,))
            i, j = q.find(len(q))[0]
            while j > 1 and (i + 1, j - 1) not in p:
                j = j - 1
            tab = p.add(i + 1, j, b)
    return tab


def irsk_inverse(tab):
    return Permutation(*InsertionAlgorithm.inverse_hecke(tab, tab)[0])


def dual_irsk(pi, n=None):
    if n is None:
        n = pi.rank
    cycles = sorted([(pi(i), i) for i in range(1, n + 1) if i <= pi(i)])
    tab = Tableau()
    for b, a in cycles:
        if a == b:
            tab = tab.add(tab.max_row() + 1, 1, a)
        else:
            p, q = InsertionAlgorithm.hecke(tab.row_reading_word() + (a,))
            i, j = q.find(len(q))[0]
            while i > 1 and (i - 1, j + 1) not in p:
                i = i - 1
            tab = p.add(i, j + 1, b)
    return tab


def _dual_irsk_inverse_helper(tab):
    if len(tab) == 0:
        return []

    n = max(tab.values())
    i, j = tab.find(n)[0]

    if j == 1:
        return _dual_irsk_inverse_helper(tab.remove(i, j)) + [(n, n)]
    else:
        tab = tab.remove(i, j)
        j = j - 1

        i = 1
        while (i + 1, j) in tab:
            i += 1

        a = tab.get(i, j)
        tab = tab.remove(i, j)
        while i > 1:
            i = i - 1
            j = 1
            while (i, j + 1) in tab and tab.get(i, j + 1) < a:
                j += 1
            a, tab = tab.get(i, j), tab.set(i, j, a)

        return _dual_irsk_inverse_helper(tab) + [(a, n)]


def dual_irsk_inverse(tab):
    ans = Permutation()
    for a, b in _dual_irsk_inverse_helper(tab):
        if a != b:
            ans *= Permutation.transposition(a, b)
    return ans


def test_dual_irsk_inverse(n=6):
    for w in Permutation.involutions(n):
        tab = dual_irsk(w, n)
        v = dual_irsk_inverse(tab)
        assert v == w


def print_involution_operation(n=6):
    for w in Permutation.involutions(n):
        tab = irsk(w).transpose()
        v = dual_irsk_inverse(tab)
        print(w, '->', v)
        print(tab)
        print()
        print()
        print()
        print()


def print_tableau_operation(n=6):
    mapping = {}
    for mu in Partition.generate(n):
        for tab in Tableau.standard(mu):
            w = irsk_inverse(tab)
            ntab = dual_irsk(w, sum(mu)).transpose()
            mapping[tab] = ntab
            print(tab)
            print('w =', w)
            print(ntab)
            #lines = str(tab).split('\n')
            #nlines = str(ntab).split('\n')
            #assert len(lines) == len(nlines)
            #print('\n'.join([lines[i] + ' -> ' + nlines[i] for i in range(len(lines))]))
            print()
            print()
            print()
            print()
    orders = {}
    for tab in mapping:
        o = 1
        x = mapping[tab]
        while x != tab:
            o += 1
            x = mapping[x]
        orders[tab] = o
    return orders, mapping


def rsk(pi, n=None):
    if n is None:
        n = pi.rank
    oneline = list(pi.oneline)
    while len(oneline) < n:
        oneline += [len(oneline) + 1]
    return InsertionAlgorithm.hecke(oneline)


def print_mn_operators(n=8):
    def wstr(w):
        return ' $\\barr{c}' + str(w) + ' \\\\ ' + w.oneline_repr() + ' \\earr$ '

    ans = []
    for n in range(2, n + 1, 2):
        for mu in Partition.generate(n):
            s = []
            i = 0
            for t in Tableau.standard(mu):
                u = dual_irsk(irsk_inverse(t), n).transpose()
                s += ['\n&\n'.join([t.tex(), u.tex()])]
                i += 1
                if i * t.max_row() >= 24:
                    s += ['\n\\end{tabular} \\newpage \\begin{tabular}{ccccc}  Tableau &  Tableau \\\\ \\hline ']
                    i = 0
            s = '\n \\\\ \\\\ \n'.join(s)
            ans += ['\\begin{tabular}{ccccc} Tableau &  Tableau \\\\ \\hline \\\\ \\\\ \n' + s + '\n\\end{tabular}']

    ans = '\n\n\\newpage\n\n'.join(ans + [''])
    with open('/Users/emarberg/Dropbox/projects/affine-transitions/notes/eric_notes/examples.tex', 'w') as f:
        f.write(ans)
        f.close()


def test_irsk(n=6):
    for w in Permutation.involutions(n):
        p, q = rsk(w)
        t = irsk(w)
        print(w)
        print(p)
        print(q)
        print(t)
        assert t == p == q


def test_dual_irsk(n=7):
    seen = {}
    for w in Permutation.involutions(n):
        t = dual_irsk(w, n)
        assert t.is_standard()
        shape = t.shape()
        print(w)
        print(t)
        print(shape)
        print()
        print()
        assert len([x for x in shape if x % 2 != 0]) == len([i for i in range(1, n + 1) if w(i) == i])
        assert t not in seen
        seen[t] = w


def des_m(pi):
    return pi.right_descent_set


def des_n(pi):
    ans = set()
    for i in range(1, pi.rank):
        s = Permutation.s_i(i)
        if len(s * pi * s) >= len(pi):
            ans.add(i)
    return ans


def dual_equivalence(tab, i):
    word = tab.row_reading_word()

    subword = tuple(a for a in word if a in [i - 1, i, i + 1])
    if subword in [(i, i + 1, i - 1), (i - 1, i + 1, i)]:
        word = tuple(a + (1 if a == i - 1 else -1 if a == i else 0) for a in word)

    elif subword in [(i, i - 1, i + 1), (i + 1, i - 1, i)]:
        word = tuple(a + (1 if a == i else -1 if a == i + 1 else 0) for a in word)

    else:
        assert subword in [(i - 1, i, i + 1), (i + 1, i, i - 1)]

    return Tableau.from_row_reading_word(word)


def representative_m(mu):
    a = 0
    w = Permutation()
    for b in Partition.transpose(mu):
        for i in range(b // 2):
            w *= Permutation.transposition(a + i + 1, a + b - i)
        a += b
    return w


def anti_representative_m(mu):
    t = Tableau({box: i + 1 for i, box in enumerate(sorted(Partition.shape(mu)))})
    oneline, _ = InsertionAlgorithm.inverse_hecke(t, t)
    return Permutation(*oneline)


def anti_representative_n(mu):
    mapping = construct_molecular_correspondence(mu)
    return mapping[anti_representative_m(mu)]


def tilde_representative_n(mu):
    shape = {
        (i, j): (i, min(j, mu[i - 1] + 1 - j)) if j != (mu[i - 1] + 1 - j) else (i, j)
        for (i, j) in Partition.shape(mu)
    }
    word = [shape[key] for key in sorted(shape, key=lambda ij: (ij[1], -ij[0]))]
    pairs = defaultdict(list)
    for i, key in enumerate(word):
        pairs[key].append(i + 1)
    w = Permutation()
    for pair in pairs.values():
        if len(pair) == 2:
            i, j = tuple(pair)
            w *= Permutation.transposition(i, j)
    return w


def representative_n(mu):
    shape = {
        (i, j): (i, min(j, mu[i - 1] + 1 - j)) if j != (mu[i - 1] + 1 - j) else (2 * ((i + 1) // 2), j)
        for (i, j) in Partition.shape(mu)
    }
    word = [shape[key] for key in sorted(shape, key=lambda ij: (ij[1], -ij[0]))]
    pairs = defaultdict(list)
    for i, key in enumerate(word):
        pairs[key].append(i + 1)
    w = Permutation()
    for pair in pairs.values():
        i, j = tuple(pair)
        w *= Permutation.transposition(i, j)
    return w

    # if mu == (2, 2, 2, 2):
    #     return Permutation(5, 6, 7, 8, 1, 2, 3, 4)
    # if mu == (3, 3, 2, 2):
    #     return Permutation(5, 6, 9, 10, 1, 2, 8, 7, 3, 4)
    # if mu == (2, 2, 2, 2, 1, 1):
    #     return Permutation(2, 1, 7, 8, 9, 10, 3, 4, 5, 6)
    # if mu == (4, 4, 2, 2):
    #     return Permutation(5, 6, 11, 12, 1, 2, 9, 10, 7, 8, 3, 4)
    # if mu == (3, 3, 3, 3):
    #     return Permutation(9, 10, 11, 12, 6, 5, 8, 7, 1, 2, 3, 4)
    # if mu == (3, 3, 2, 2, 1, 1):
    #     return Permutation(2, 1, 7, 8, 11, 12, 3, 4, 10, 9, 5, 6)
    # if mu == (2, 2, 2, 2, 2, 2):
    #     return Permutation(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)
    # if mu == (2, 2, 2, 2, 1, 1, 1, 1):
    #     return Permutation(2, 1, 4, 3, 9, 10, 11, 12, 5, 6, 7, 8)
    # if mu == (5, 5, 2, 2):
    #     return Permutation(5, 6, 13, 14, 1, 2, 11, 12, 10, 9, 7, 8, 3, 4)
    # if mu == (3, 3, 3, 3, 3, 3):
    #     return Permutation(13, 14, 15, 16, 17, 18, 8, 7, 10, 9, 12, 11, 1, 2, 3, 4, 5, 6)
    # if mu == (4, 4, 3, 3):
    #     return Permutation(9, 10, 13, 14, 6, 5, 11, 12, 1, 2, 7, 8, 3, 4)
    # if mu == (5, 5, 3, 3):
    #     return Permutation(9, 10, 15, 16, 6, 5, 13, 14, 1, 2, 12, 11, 7, 8, 3, 4)
    # if mu == (4, 4, 3, 3, 1, 1):
    #     return Permutation(2, 1, 11, 12, 15, 16, 8, 7, 13, 14, 3, 4, 9, 10, 5, 6)
    # if mu == (3, 3, 2, 2, 2, 2):
    #     q = sorted([w for w in Permutation.fpf_involutions(14) if des_n(w) == des_m(representative_m(mu))])
    #     x = q[0]
    #     print(q)
    #     return x

    # a = 0
    # w = Permutation()
    # assert len(mu) % 2 == 0
    # mu = [mu[i] for i in range(0, len(mu), 2)]
    # for b in mu:
    #     for i in range(b):
    #         w *= Permutation.transposition(a + i + 1, a + 2 * b - i)
    #         if i % 2 != 0:
    #             s = Permutation.s_i(a + i)
    #             w = s * w * s
    #     a += 2 * b
    # return w.star()


def bidirected_edges_m(w):
    n = w.rank

    for i in range(2, n):
        s = Permutation.s_i(i - 1)
        t = Permutation.s_i(i)

        y = s * w * s
        if len(t * w * t) <= len(w) < len(y) < len(t * y * t):
            yield y, i

        if len(t * w * t) > len(w) > len(y) >= len(t * y * t):
            yield y, i

        y = t * w * t
        if len(s * w * s) <= len(w) < len(y) < len(s * y * s):
            yield y, i

        if len(s * w * s) > len(w) > len(y) >= len(s * y * s):
            yield y, i


def molecule_m(mu):
    w = representative_m(mu)
    ans = set()
    level = {w}
    while level:
        ans |= level
        level = {z for y in level for z, _ in bidirected_edges_m(y) if z not in ans}
    return ans


def get_molecules_m(n, verbose=False):
    ans = {}
    nn = double_fac(n - 1)
    for i, w in enumerate(Permutation.fpf_involutions(n)):
        mu = rsk(w)[0].shape()
        if mu not in ans:
            ans[mu] = set()
        ans[mu].add(w)

        if verbose:
            a = nn - i
            if a % 100 == 0:
                print(a)
    return ans


def bidirected_edges_n(w):
    assert w.is_fpf_involution()
    n = w.rank

    for i in range(2, n):
        s = Permutation.s_i(i - 1)
        t = Permutation.s_i(i)

        y = s * w * s
        if len(t * w * t) < len(w) < len(y) <= len(t * y * t):
            yield y, i

        if len(t * w * t) >= len(w) > len(y) > len(t * y * t):
            yield y, i

        y = t * w * t
        if len(s * w * s) < len(w) < len(y) <= len(s * y * s):
            yield y, i

        if len(s * w * s) >= len(w) > len(y) > len(s * y * s):
            yield y, i


def molecule_n(mu):
    w = representative_n(mu)
    ans = set()
    level = {w}
    while level:
        ans |= level
        level = {z for y in level for z, _ in bidirected_edges_n(y) if z not in ans}
    return ans


def get_molecules_n(n, verbose=False):
    ans = {}
    for mu in Partition.generate(n, even_parts=True):
        mu = Partition.transpose(mu)
        molecule = molecule_n(mu)
        ans[mu] = molecule
    return ans


def double_fac(n):
    if n <= 2:
        return 1
    return double_fac(n - 2) * n


def construct_molecular_correspondence(mu):
    if mu in MAPPING_CACHE:
        return MAPPING_CACHE[mu]
    mapping = {}
    w = representative_m(mu)
    mapping[w] = representative_n(mu)
    level = {w}
    a = 0
    while level:
        nextlevel = set()
        for w in level:
            edges_m = {i: y for y, i in bidirected_edges_m(w)}
            edges_n = {i: y for y, i in bidirected_edges_n(mapping[w])}
            if set(edges_m) != set(edges_n):
                print()
                print('shape:', mu)
                print('level:', a)
                print()
                print(w, set(edges_m), '?=', set(edges_n), mapping[w])
                print()
                assert set(edges_m) == set(edges_n)
            for i in set(edges_m) & set(edges_n):
                y = edges_m[i]
                z = edges_n[i]
                if y not in mapping:
                    nextlevel.add(y)
                    mapping[y] = z
                else:
                    assert mapping[y] == z
        level = nextlevel
        a += 1
    MAPPING_CACHE[mu] = mapping
    return mapping


def descent_dicts(n):
    ndes = {}
    mdes = {}
    for w in Permutation.involutions(n):
        d = tuple(sorted(des_m(w)))
        if d not in mdes:
            mdes[d] = {w}
        else:
            mdes[d].add(w)

        d = tuple(sorted(des_n(w)))
        if d not in mdes:
            ndes[d] = {w}
        else:
            ndes[d].add(w)
    return mdes, ndes


def test_get_molecules_m(n=8):
    ans = {}
    inv = set(Permutation.fpf_involutions(n))
    while inv:
        w = inv.pop()
        mu = rsk(w)[0].shape()
        molecule = molecule_m(mu)
        inv -= molecule
        assert mu not in ans
        ans[mu] = molecule
    bns = get_molecules_m(n)
    assert all(ans[mu] == bns[mu] for mu in bns)
    assert all(representative_m(mu) in bns[mu] for mu in bns)


def test_bidirected_edges_m(n=8):
    for w in Permutation.fpf_involutions(n):
        p, _ = rsk(w)
        for y, i in set(bidirected_edges_m(w)):
            print(w, '<--', i, '-->', y)
            q, _ = rsk(y)
            print(p)
            print(q)
            print()
            assert q == dual_equivalence(p, i)


def print_molecular_correspondence(n=8):
    def wstr(w):
        return ' $\\barr{c}' + str(w) + ' \\\\ ' + w.oneline_repr() + ' \\earr$ '

    ans = []
    for n in range(2, n + 1, 2):
        for mu in Partition.generate(n, even_parts=True):
            s = []
            mu = Partition.transpose(mu)
            mapping = construct_molecular_correspondence(mu)
            i = 0
            for v in mapping:
                w = mapping[v]
                # w = w * Permutation.longest_element(n)
                t, _ = rsk(v)
                p, q = rsk(w, n)
                s += ['\n&\n'.join([wstr(v), t.tex(), wstr(w), dual_irsk(w).tex()])]
                i += 1
                if i * t.max_row() >= 24:
                    s += ['\n\\end{tabular} \\newpage \\begin{tabular}{ccccc}\nM & Tableau &  N & & \\\\ \\hline ']
                    i = 0
            s = '\n \\\\ \\\\ \n'.join(s)
            ans += ['\\begin{tabular}{ccccc}\nM & Tableau &   N & & \\\\ \\hline \\\\ \\\\ \n' + s + '\n\\end{tabular}']

    ans = '\n\n\\newpage\n\n'.join(ans + [''])
    with open('/Users/emarberg/Dropbox/projects/affine-transitions/notes/eric_notes/examples.tex', 'w') as f:
        f.write(ans)
        f.close()


def test_molecular_correspondence(n=8):
    for mu in Partition.generate(n, even_parts=True):
        mu = Partition.transpose(mu)
        mapping = construct_molecular_correspondence(mu)
        for v in mapping:
            w = mapping[v]
            assert irsk(v).transpose() == dual_irsk(w)


def test_get_molecules_n(n=8):
    ans = {}
    fpf = set(Permutation.fpf_involutions(n))
    bns = get_molecules_n(n)
    for mu in bns:
        molecule = molecule_n(mu)
        fpf -= molecule
        ans[mu] = molecule
    assert len(fpf) == 0
    assert all(ans[mu] == bns[mu] for mu in bns)


def test_tilde_representatives(n=8):
    for nu in Partition.generate(n):
        w = tilde_representative_n(nu)
        v = anti_representative_m(Partition.transpose(nu))
        assert w == v


def test_star(n=8):
    for nu in Partition.generate(n, even_parts=True):
        mu = Partition.transpose(nu)
        lu = Partition.transpose(tuple(2 * mu[i] for i in range(0, len(mu), 2)))
        assert mu == Partition.transpose(tuple(2 * lu[i] for i in range(0, len(lu), 2)))


def test_anti_representatives(n=8):
    for nu in Partition.generate(n, even_parts=True):
        mu = Partition.transpose(nu)
        lu = Partition.transpose(tuple(2 * mu[i] for i in range(0, len(mu), 2)))
        print(mu)
        print(lu)
        print()

        assert representative_m(lu) == anti_representative_n(mu)

        # a = representative_m(mu)
        # b = representative_m(lu)
        # c = representative_n(mu)
        # d = representative_n(lu)

        # print('a =', a, a.oneline_repr())
        # print('b =', b, b.oneline_repr())
        # print('c =', c, c.oneline_repr())
        # print('d =', d, d.oneline_repr())
        # print()

        # a = anti_representative_m(mu)
        # b = anti_representative_m(lu)
        # c = anti_representative_n(mu)
        # d = anti_representative_n(lu)

        # print('a =', a, a.oneline_repr())
        # print('b =', b, b.oneline_repr())
        # print('c =', c, c.oneline_repr())
        # print('d =', d, d.oneline_repr())
        # print()
