from crystals import(
    AbstractGLCrystal,
    AbstractQCrystal,
    AbstractPrimedQCrystal,
)
from keys import (
    sorting_permutation,
    skew_symmetrize_strict_partition,
    symmetrize_strict_partition,
    strict_half_partition,
    weak_half_partition,
    p_key,
    q_key,
    key,
    skew_symmetric_double,
    symmetric_double,
)
from tests.test_keys import try_to_decompose_p, try_to_decompose_q, decompose_p, decompose_q
from partitions import Partition
from symmetric import (
    SchurQ,
    SchurP,
    Schur
)
from tableaux import Tableau
from permutations import Permutation
from schubert import X
from words import Word, weak_eg_insert, eg_insert, fpf_insert, involution_insert
from keys import decompose_into_keys
from tests.test_keys import try_to_decompose_q, try_to_decompose_p
import random
import time


PRINT_DIR = "/Users/emarberg/Downloads/"
FPF_DEMAZURE_TABLEAU_CACHE = {}
INV_DEMAZURE_TABLEAU_CACHE = {}


def print_fpf_demazure_tableau_table(n, thresh=40, columns=3, numtabs=1):
    # last setting: print_fpf_demazure_tableau_table(8, 68, 4, numtabs=1)
    def sorter(r):
        p, q, a, h = r
        mu = p.partition().tuple()
        return (len(mu), mu, tuple(-x for x in a))
    rows = test_fpf_demazure_tableau(n, include_dominant=True)
    rows.sort(key=sorter)
    caption = 'Some $\\Sp$-reduced tableaux with the weak compositions predicted in Conjecture~\\ref{sp-demazure-conj2}.'
    printers = numtabs * (print_tableau,) + (print_composition,)
    if numtabs == 1:
        rows = [(p, a) for (p, q, a, h) in rows]
    elif numtabs == 2:
        rows = [(p, q, a) for (p, q, a, h) in rows]
    elif numtabs == 3:
        rows = [(p, q, eg_insert(*h)[0].transpose(), a) for (p, q, a, h) in rows]
    else:
        raise Exception
    print_table(rows, printers, thresh, columns, caption, 'Sp')


def print_inv_demazure_tableau_table(n, thresh=40, columns=3, numtabs=1):
    # last setting: print_inv_demazure_tableau_table(7, 74, 4, numtabs=1)
    def sorter(r):
        p, q, a, h = r
        mu = p.partition().tuple()
        return (len(mu), mu, tuple(-x for x in a))
    rows = test_inv_demazure_tableau(n, include_dominant=True, exclude_one_row=True)
    rows.sort(key=sorter)
    caption = 'Some $\\O$-reduced tableaux with the weak compositions predicted in Conjecture~\\ref{o-demazure-conj2}.'
    printers = numtabs * (print_tableau,) + (print_composition,)
    if numtabs == 1:
        rows = [(p, a) for (p, q, a, h) in rows]
    elif numtabs == 2:
        rows = [(p, q, a) for (p, q, a, h) in rows]
    elif numtabs == 3:
        rows = [(p, q, eg_insert(*h)[0].transpose(), a) for (p, q, a, h) in rows]
    else:
        raise Exception
    print_table(rows, printers, thresh, columns, caption, 'O')


def print_table(rows, printers, thresh, columns, caption, tag):
    def convert(tables, ans):
        if ans:
            pre = ['\\ytableausetup{smalltableaux}']
            pre += ['\\begin{tabular}[t]{%s}' % (len(printers) * 'l')]
            gaps = ' '.join((len(printers) - 1) * ['&'])
            pre += ['$T$ %s $\\alpha^{\\%s}(T)$ \\\\ \\hline \\\\[-6pt]' % (gaps, tag)]
            ans = ['\n\\\\[-6pt] \\\\\n'.join(ans)]
            ans = pre + ans + ['\\end{tabular}'] 
            ans = '\n'.join(ans)   
            tables.append(ans) 
    tables = []
    space = 2
    count = -space
    ans = []
    for row in rows:
        ans += [[]]
        for a, f in zip(row, printers):
            ans[-1] += [f(a)]
        delta = max([len(x.split('\n')) for x in ans[-1]]) + space
        ans[-1] = ' & '.join(ans[-1])
        if count + delta > thresh:
            convert(tables, ans[:-1])
            ans = [ans[-1]]
            count = -space
        count += delta
    convert(tables, ans)
    ans = []
    while tables:
        curr_columns = min(len(tables), columns)
        ans += ['\\begin{figure}[h]']
        ans += ['\\begin{center}']
        ans += ['\\begin{tabular}{%s}' % '|'.join(curr_columns * ['l'])]
        ans += ['\n&\n'.join(tables[:curr_columns])]
        ans += ['\\end{tabular}']
        ans += ['\\end{center}']
        ans += ['\\caption{%s}' % caption]
        ans += ['\\end{figure}']
        ans += []
        tables = tables[curr_columns:]
    ans = '\n'.join(ans)
    with open(PRINT_DIR + 'tables_' + tag + '.tex', 'w') as f:
        f.write(ans)


def print_tableau(t):
    rows = []
    for i in range(1, t.max_row + 1):
        row = []
        for j in range(1, t.max_column + 1):
            v = t.entry(i, j)
            row += [str(v) if v is not None else '\\none']
        rows += [' & '.join(row)]
    return '$\\begin{ytableau}' + ' \\\\\n'.join(rows) + '\\end{ytableau}$'


def print_composition(alpha):
    while alpha and alpha[-1] == 0:
        alpha = alpha[:-1]
    assert alpha
    return '{\\footnotesize$' + ''.join(map(str, alpha)) + '$}'
    # return '$\\emptyset$'


def test_inv_odd_almost_highest(n=3):
    for permutation_size in range(7):
        for w in Permutation.involutions(permutation_size):
            crystal = AbstractQCrystal.from_involution(w, n, increasing=False)
            for f in crystal:
                i = 0
                while i + 1 < n and crystal.e_operator(i + 1, f) is None:
                    i += 1
                if i == 0:
                    continue
                j = 1
                while j <= i and crystal.e_operator(-j, f) is None:
                    j += 1
                if i < j:
                    continue

                g = crystal.e_operator(-j, f)
                print('i =', i, 'j =', j)
                print(f, '---', -j, '--->', g)
                print()
                print(eg_insert(*f)[0])
                print()
                print(eg_insert(*g)[0])
                # input('')
                assert all(crystal.e_operator(k, g) is None for k in range(1, j + 1))


            
def test_fpf_odd_almost_highest(n=3):
    pass


def test_highest_lowest(n=3, k=6):
    for w in Permutation.all(k):
        crystal = AbstractGLCrystal.from_permutation(w, n, increasing=False)
        highest = [f for f in crystal if crystal.is_highest_weight(f)]
        lowest = [f for f in crystal if crystal.is_lowest_weight(f)]
        
        for f in highest:
            tab = eg_insert(*f)[0]
            expected = tuple(tuple(reversed(c)) for c in tab.get_columns())
            expected += (n - len(expected)) * ((),)
            assert f == expected

        for f in lowest:
            fstar = tuple(tuple(k - i for i in part) for part in f)
            tab = eg_insert(*fstar)[0]
            tab = tab.__class__({b: k - v.number for b, v in tab.mapping.items()})
            expected = tuple(tuple(r) for r in reversed(tab.get_rows()))
            expected = (n - len(expected)) * ((),) + expected
            assert f == expected


def test_highest_lowest_fpf(n=3, k=6):
    for w in Permutation.fpf_involutions(k):
        crystal = AbstractQCrystal.from_fpf_involution(w, n, increasing=False)
        highest = [f for f in crystal if crystal.is_highest_weight(f)]
        lowest = [f for f in crystal if crystal.is_lowest_weight(f)]
        
        for f in highest:
            tab = fpf_insert(*f)[0]
            # expected = tuple(tuple(reversed(c)) for c in tab.get_columns())
            # expected += (n - len(expected)) * ((),)
            print(tab)
            print(f)
            # input('\n?\n')

        for f in lowest:
            fstar = tuple(tuple(k - i for i in part) for part in f)
            tab = fpf_insert(*fstar)[0]
            t = tab
            tab = tab.__class__({b: k - v.number for b, v in tab.mapping.items()})
            expected = tuple(tuple(r) for r in reversed(tab.get_rows()))
            expected = (n - len(expected)) * ((),) + expected


def test_highest_lowest_inv(n=3, k=6):
    for w in Permutation.involutions(k):
        crystal = AbstractQCrystal.from_involution(w, n, increasing=False)
        highest = [f for f in crystal if crystal.is_highest_weight(f)]
        lowest = [f for f in crystal if crystal.is_lowest_weight(f)]
        
        for f in highest:
            tab = involution_insert(*f)[0]
            # expected = tuple(tuple(reversed(c)) for c in tab.get_columns())
            # expected += (n - len(expected)) * ((),)
            print(tab)
            print(f)
            # input('\n?\n')

        for f in lowest:
            fstar = tuple(tuple(k - i for i in part) for part in f)
            tab = involution_insert(*fstar)[0]
            t = tab
            tab = tab.__class__({b: k - v.number for b, v in tab.mapping.items()})
            expected = tuple(tuple(r) for r in reversed(tab.get_rows()))
            expected = (n - len(expected)) * ((),) + expected
            

def test_e0_preserves_strings(n=3, k=6):
    for w in Permutation.involutions(k):
        crystal = AbstractPrimedQCrystal.from_involution(w, n, increasing=False)
        for f in crystal:
            for j in range(0, n * n, n):
                g = crystal.e_operator(j, f)
                assert g is None or all(crystal.e_string(-i, f) == crystal.e_string(-i, g) for i in range(1, n))
                assert g is None or all(crystal.f_string(-i, f) == crystal.f_string(-i, g) for i in range(1, n))
                        

def factorization_character(subset):
    ans = 0
    for v in subset:
        t = X(0)**0
        for i, a in enumerate(v):
            t *= X(i + 1)**len(a)
        ans += t
    return ans


def flags(n):
    seed = tuple(range(1, n + 1))
    level = {seed}
    while level:
        nextlevel = set()
        for f in level:
            yield f
            for i in range(n - 1):
                if f[i] + 1 <= min(f[i + 1], n):
                    g = f[:i] + (f[i] + 1,) + f[i + 1:]
                    nextlevel.add(g)
        level = nextlevel


def partitions(max_num_parts, max_level=None):
    level = {max_num_parts * (0,)}
    while max_level is None or sum(next(iter(level))) <= max_level:
        nextlevel = set()
        for mu in level:
            yield mu
            mu = list(mu)
            for i in range(max_num_parts):
                if i == 0 or mu[i] + 1 <= mu[i - 1]:
                    mu[i] += 1
                    nextlevel.add(tuple(mu))
                    mu[i] -= 1
        level = nextlevel


def strict_partitions(max_num_parts, max_level=None):
    level = {max_num_parts * (0,)}
    while max_level is None or sum(next(iter(level))) <= max_level:
        nextlevel = set()
        for mu in level:
            yield mu
            mu = list(mu)
            for i in range(max_num_parts):
                if i == 0 or mu[i] + 1 < mu[i - 1]:
                    mu[i] += 1
                    nextlevel.add(tuple(mu))
                    mu[i] -= 1
        level = nextlevel


def partition_permutations(mu, n):
    seen = {mu}
    level = {(mu, None, mu)}
    while level:
        nextlevel = set()
        for mu, i, nu in level:
            yield (mu, i, nu)
            for j in range(1, n):
                ku = nu[:j - 1] + (nu[j], nu[j - 1]) + nu[j + 1:]
                if ku not in seen:
                    nextlevel.add((nu, j, ku))
                    seen.add(ku)
        level = nextlevel


def emax(crystal, index, vertex):
    while crystal.e_operator(index, vertex) is not None:
        vertex = crystal.e_operator(index, vertex)
    return vertex


def restrict_variables(ans, n):
    while any(i > n for i in ans.variables()):
        ans = ans.set(max(ans.variables()), 0)
    return ans


def inv_negative_one_operator_test(crystal, subset):
    if crystal.rank == 1:
        return True
    for b in crystal:
        c = crystal.e_operator(-1, b)
        if b in subset and c is not None and c not in subset:
            return False
        if c is None or c == crystal.e_operator(1, b):
            continue
        if b not in subset and c in subset:
            return False
    return True


def inv_zero_operator_test(crystal, subset):
    n = crystal.rank
    for b in subset:
        for i in range(0, n * n, n):
            c = crystal.f_operator(i, b)
            if c is not None and c not in subset:
                return False
            c = crystal.e_operator(i, b)
            if c is not None and c not in subset:
                return False
    return True


def inv_is_bounded(f, flag=None):
    f = tuple(tuple(map(abs, a)) for a in f)
    return _is_bounded(f, flag)


def _is_bounded(f, flag=None):
    def phi(a):
        assert a >= 1
        return flag[a - 1] if flag is not None and a <= len(flag) else a

    return all(len(a) == 0 or all(i + 1 <= phi(x) for x in a) for i, a in enumerate(f))


def fpf_is_bounded(f, flag=None):
    return _is_bounded(f, flag)


def draw_demazure(alpha):
    n = len(alpha)
    mu = list(sorted(alpha, reverse=True))
    w_mu = Permutation.from_shape(*mu)
    crystal = AbstractGLCrystal.from_permutation(w_mu, n, increasing=False)
    brf = [f for f in crystal if _is_bounded(f)]
    for i in sorting_permutation(alpha):
        brf = [f for f in crystal if emax(crystal, i, f) in brf]
    crystal.draw(highlighted_nodes=brf, extended=True)


def draw_inv_demazure(alpha, n=None, exclude=False):
    if type(alpha) == Permutation:
        w_mu = alpha
        s = []
        assert n is not None
    else:
        n = len(alpha) if n is None else n
        mu = list(sorted(alpha, reverse=True))
        for i in range(len(mu)):
            mu[i] = mu[i] - i if mu[i] > i else 0
        w_mu = Permutation.from_involution_shape(*mu)
        print(w_mu.cycle_repr())
        s = sorting_permutation(alpha)
    crystal = AbstractPrimedQCrystal.from_involution(w_mu, n, increasing=False)
    brf = [f for f in crystal if inv_is_bounded(f)]
    for i in s:
        brf = [f for f in crystal if emax(crystal, i, f) in brf]
    # crystal.draw(highlighted_nodes=brf, tex=True)
    crystal.draw(highlighted_nodes=brf, extended=True, exclude=exclude)


def draw_fpf_demazure(alpha, n=None, exclude=False):
    if type(alpha) == Permutation:
        w_mu = alpha
        s = []
        assert n is not None
    else:
        n = len(alpha) if n is None else n
        mu = list(sorted(alpha, reverse=True))
        for i in range(len(mu)):
            mu[i] = mu[i] - i - 1 if mu[i] > i + 1 else 0
        w_mu = Permutation.from_fpf_involution_shape(*mu)
        print(w_mu.cycle_repr())
        s = sorting_permutation(alpha)
    crystal = AbstractQCrystal.from_fpf_involution(w_mu, n, increasing=False)
    brf = [f for f in crystal if fpf_is_bounded(f)]
    for i in s:
        brf = [f for f in crystal if emax(crystal, i, f) in brf]
    # crystal.draw(highlighted_nodes=brf, tex=True)
    crystal.draw(highlighted_nodes=brf, extended=True, exclude=exclude)


def generate_demazure(mu, dictionary):
    n = len(mu)
    alpha = mu
    if alpha in dictionary:
        return

    w_mu = Permutation.from_shape(*mu)
    crystal = AbstractGLCrystal.from_permutation(w_mu, n, increasing=False)
    brf = [f for f in crystal if _is_bounded(f)]

    subsets = {alpha: brf}
    for ku, i, nu in partition_permutations(alpha, n):
        if i is not None:
            subsets[nu] = [f for f in crystal if emax(crystal, i, f) in subsets[ku]]
        if nu not in dictionary:
            dictionary[nu] = crystal.truncate(subsets[nu])


def generate_inv_demazure(mu, dictionary):
    n = len(mu)
    alpha = symmetrize_strict_partition(mu, n)
    if alpha in dictionary:
        return

    w_mu = Permutation.from_involution_shape(*mu)
    crystal = AbstractPrimedQCrystal.from_involution(w_mu, n, increasing=False)
    brf = [f for f in crystal if inv_is_bounded(f)]

    subsets = {alpha: brf}
    for ku, i, nu in partition_permutations(alpha, n):
        if i is not None:
            subsets[nu] = [f for f in crystal if emax(crystal, i, f) in subsets[ku]]
        if nu not in dictionary:
            dictionary[nu] = crystal.truncate(subsets[nu])


def quick_generate_inv_demazure(alpha, dictionary):
    if alpha in dictionary:
        return

    sorter = sorting_permutation(alpha)
    if len(sorter) == 0:
        n = len(alpha)
        mu = weak_half_partition(alpha)
        w_mu = Permutation.from_involution_shape(*mu)
        h = w_mu.get_involution_word()
        t = inv_decreasing_tableau(*h)
        a = tuple(tuple(_) for _ in reversed(t.get_rows()))
        a += (n - len(a)) * ((),)

        t0 = time.time()
        crystal = AbstractPrimedQCrystal.from_inv_factorization(a, n, increasing=False)
        brf = crystal.truncate([f for f in crystal if inv_is_bounded(f)])
        t1 = time.time()

        print('   ', 'created %s:' % str(alpha), t1 - t0)
        dictionary[alpha] = (brf, crystal)
    else:
        i = sorter[0]
        gamma = alpha[:i - 1] + (alpha[i], alpha[i - 1]) + alpha[i + 1:]
        quick_generate_inv_demazure(gamma, dictionary)
        brf, crystal = dictionary[gamma]
        crf = [f for f in crystal if emax(crystal, i, f) in brf]
        dictionary[alpha] = (crystal.truncate(crf), crystal)


def generate_fpf_demazure(mu, dictionary):
    n = len(mu)
    alpha = skew_symmetrize_strict_partition(mu, n)
    if alpha in dictionary:
        return

    w_mu = Permutation.from_fpf_involution_shape(*mu)
    crystal = AbstractQCrystal.from_fpf_involution(w_mu, n, increasing=False)
    brf = [f for f in crystal if fpf_is_bounded(f)]

    subsets = {alpha: brf}
    for ku, i, nu in partition_permutations(alpha, n):
        if i is not None:
            subsets[nu] = [f for f in crystal if emax(crystal, i, f) in subsets[ku]]
        if nu not in dictionary:
            # print('  ', nu)
            dictionary[nu] = crystal.truncate(subsets[nu])


def quick_generate_fpf_demazure(alpha, dictionary):
    if alpha in dictionary:
        return

    sorter = sorting_permutation(alpha)
    if len(sorter) == 0:
        n = len(alpha)
        mu = strict_half_partition(alpha)
        w_mu = Permutation.from_fpf_involution_shape(*mu)
        h = w_mu.get_fpf_involution_word()
        t = fpf_decreasing_tableau(*h)
        a = tuple(tuple(_) for _ in reversed(t.get_rows()))
        a += (n - len(a)) * ((),)

        t0 = time.time()
        crystal = AbstractQCrystal.from_fpf_factorization(a, n, increasing=False)
        brf = crystal.truncate([f for f in crystal if fpf_is_bounded(f)])
        t1 = time.time()

        print('   ', 'created %s:' % str(alpha), t1 - t0)
        dictionary[alpha] = (brf, crystal)
    else:
        i = sorter[0]
        gamma = alpha[:i - 1] + (alpha[i], alpha[i - 1]) + alpha[i + 1:]
        quick_generate_fpf_demazure(gamma, dictionary)
        brf, crystal = dictionary[gamma]
        crf = [f for f in crystal if emax(crystal, i, f) in brf]
        dictionary[alpha] = (crystal.truncate(crf), crystal)


def find_isomorphism(target, highest, dictionary):
    n = target.rank
    ans = {}
    pairs = []
    for h in highest:
        found = False
        for alpha in dictionary:
            if target.isomorphic_highest_weight_crystals(target, h, dictionary[alpha]):
                ans[alpha] = ans.get(alpha, 0) + 1
                pairs.append((h, alpha))
                found = True
                break
        if not found:
            return None
    return ans, pairs


def get_expected_ch(decomposition, fn):
    expected_ch = 0
    if decomposition is not None:
        for alpha, coeff in decomposition.items():
            expected_ch += coeff *  fn(alpha)
    return expected_ch


def verify_fpf_string_lengths(crystal, demazure):
    n = crystal.rank
    for b in demazure:
        for i in range(1, n):
            wt = demazure.weight(b)
            if demazure.f_operator(i, b) is not None or demazure.e_operator(i, b) is not None:
                if wt[i - 1] - wt[i] != demazure.f_string(i, b) - demazure.e_string(i, b):
                    print('')
                    print('failure at element b =', b, 'index i =', i, 'of', demazure.extended_indices)
                    print('')
                    raise Exception
        for i in range(1, n):
            if all(crystal.e_operator(j, b) is None for j in range(-i + 1, i + 1) if j != 0) and crystal.e_operator(-i, b) is not None:
                assert crystal.e_operator(-i, b) in demazure


def verify_inv_string_lengths(crystal, demazure):
    n = crystal.rank
    for b in demazure:
        for i in demazure.extended_indices:
            wt = demazure.weight(b)
            if i >= 0 and i % demazure.rank == 0:
                j = i // demazure.rank
                if wt[j] != 0 and demazure.f_string(i, b) + demazure.e_string(i, b) != 1:
                    print('')
                    print('failure at element b =', b, 'index i =', i, 'of', demazure.extended_indices)
                    print('')
                    raise Exception
            elif i > 0:
                if demazure.f_operator(i, b) is not None or demazure.e_operator(i, b) is not None:
                    if wt[i - 1] - wt[i] != demazure.f_string(i, b) - demazure.e_string(i, b):
                        print('')
                        print('failure at element b =', b, 'index i =', i, 'of', demazure.extended_indices)
                        print('')
                        raise Exception
        for i in range(1, n):
            if all(crystal.e_operator(j, b) is None for j in range(-i + 1, i + 1)) and crystal.e_operator(-i, b) is not None:
                assert crystal.e_operator(-i, b) in demazure


def verify_string_lengths(crystal):
    for b in crystal:
        for i in crystal.indices:
            wt = crystal.weight(b)
            if crystal.f_operator(i, b) is not None or crystal.e_operator(i, b) is not None:
                if wt[i - 1] - wt[i] != crystal.f_string(i, b) - crystal.e_string(i, b):
                    print('')
                    print('failure at element b =', b, 'index i =', i, 'of', crystal.indices)
                    print('')
                    raise Exception


def _test_results(results, flag, decomposition, brf, crystal):
    def s(a, i):
        return a if a[i] <= a[i + 1] else a[:i] + (a[i + 1], a[i]) + a[i + 2:]

    def swap(a, i, j):
        for x in range(i, j):
            for y in range(i, j):
                a = s(a, y)
        return a

    results[flag] = (decomposition, brf)
    standard = tuple(i + 1 for i in range(len(flag)))
    expected, _ = results[standard]
    word = []
    newflag = flag[:]
    while True:
        i = [i for i in range(len(newflag) - 1) if newflag[i] == newflag[i + 1]]
        if not i:
            break
        i = i[0]
        word.append(i + 1)
        newflag = newflag[:i] + (newflag[i] - 1,) + newflag[i + 1:]

        new_expected = {}
        for a in expected:
            b = s(a, i)
            new_expected[b] = new_expected.get(b, 0) + expected[a]
        expected = new_expected

    if expected != decomposition:
        print('  predecessor:', results[standard][0], decomposition, flag, expected)
        input('\n?\n')
    print()
    print('  word =', tuple(reversed(word)))
    print()

    for i in range(len(flag)):
        if flag[i] > i + 1:
            newflag = list(flag)
            newflag[i] -= 1
            newflag = tuple(newflag)
            j = flag[i] - 1
            print(i + 1, j, flag, newflag)
            print()
            newdecomp = {}
            for a, v in results[newflag][0].items():
                b = s(a, j - 1)
                newdecomp[b] = newdecomp.get(b, 0) + v
            assert newdecomp == decomposition

            _, b = results[newflag]
            a = crystal.truncate([f for f in crystal if emax(crystal, j, f) in b])
            if sorted(a) != sorted(brf):
                print()
                print(list(a))
                print(list(brf))
                print()
                assert False
            break


def test_demazure_generic(n=2, permutation_size=5):
    demazure = {}
    is_bounded = _is_bounded
    for w in Permutation.all(permutation_size):
        print(5 * '\n')
        crystal = AbstractGLCrystal.from_permutation(w, n, increasing=False)
        results = {}
        for flag in flags(n):
            brf = crystal.truncate([f for f in crystal if is_bounded(f, flag)])
            highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
            for f in highest:
                generate_demazure(crystal.weight(f), demazure)
            decomposition, pairs = find_isomorphism(brf, highest, demazure)

            ch = factorization_character(brf)
            expected_ch = get_expected_ch(decomposition, key)
            print(w.oneline_repr(permutation_size), 'flag =', flag, ch == expected_ch, decomposition)
            for h, alpha in pairs:
                tab, _ = weak_eg_insert(*[a for f in reversed(h) for a in reversed(f)])
                #print(h)
                if all(flag[i] == i + 1 for i in range(len(flag))):
                    print(tab)
                #print(alpha)
                #print()
            try:
                assert ch == expected_ch and decomposition is not None
            except:
                crystal.draw(highlighted_nodes=brf)
                print()
                print('ch =', ch)
                print('ex =', expected_ch)
                input('\n?\n')
            _test_results(results, flag, decomposition, brf, crystal)


def test_demazure(n=2, limit=8):
    is_bounded = _is_bounded
    for mu in partitions(n, limit):
        w_mu = Permutation.from_shape(*mu)
        crystal = AbstractGLCrystal.from_permutation(w_mu, n, increasing=False)
        brf = crystal.truncate([f for f in crystal if is_bounded(f)])

        highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
        assert len(highest) == 1

        demazure = {mu: brf}
        for ku, i, nu in partition_permutations(mu, n):
            if i is not None:
                demazure[nu] = crystal.truncate([f for f in crystal if emax(crystal, i, f) in demazure[ku]])
            ch = factorization_character(demazure[nu])
            expected_ch = key(nu) #restrict_variables(q_key(nu), n)
            
            print(nu, w_mu)
            assert ch == expected_ch
            verify_string_lengths(demazure[nu])
            # crystal.draw(highlighted_nodes=demazure[nu])
            # input('')


def inv_decreasing_highest(decreasing_tab):
    m = max([0] + [abs(i) for i in decreasing_tab.row_reading_word()]) + 1
    def invert(i):
        return m - i if i > 0 else -(m + i)
    t = Tableau({b: invert(v.number) for (b, v) in decreasing_tab.mapping.items()})   
    ans = Tableau.inverse_involution_insertion(t, Tableau({(i, j): i for (i, j) in t.mapping}))
    ans = tuple(tuple(invert(x) for x in _) for _ in ans)
    return ans


def inv_decreasing_tableau(*word):
    word = [(w,) if type(w) == int else w for w in word]
    m = max([0] + [abs(i) for _ in word for i in _]) + 1
    def invert(i):
        return (m - i if i > 0 else -(m + i)) if type(i) == int else tuple(invert(_) for _ in i)
    t, _ = involution_insert(*(invert(a) for a in word))
    ans = Tableau({b: invert(v.number) for (b, v) in t.mapping.items()})
    return ans


def inv_increasing_tableau(*word):
    return involution_insert(*word)[0]


def test_inv_demazure_tableau(permutation_size=4, include_dominant=False, exclude_one_row=False):
    key = (permutation_size, include_dominant, exclude_one_row)
    if key in INV_DEMAZURE_TABLEAU_CACHE:
        return INV_DEMAZURE_TABLEAU_CACHE[key]

    ans = []
    invdemazure = {}
    is_bounded = inv_is_bounded
    count = 0 
    dominant = 0
    total = len(list(Permutation.involutions(permutation_size)))
    for z in Permutation.involutions(permutation_size):
        total -= 1
        print()
        print()
        print()
        print('#', total, ':', 'z =', z.cycle_repr(), 'dominant?', z.is_dominant())
        if z.is_dominant():
            count += 1
            dominant += 1
            h = z.get_involution_word()
            t, u = inv_decreasing_tableau(*h), inv_increasing_tableau(*h)
            highest = inv_decreasing_highest(t)
            mu = t.partition().tuple()
            alpha = symmetric_double(mu)
            if len(z) > 0 and include_dominant and (not exclude_one_row or t.num_rows() > 1):
                ans.append((t, u, alpha, highest))
            print(t)
            print()
            print('rank = N/A', 'alpha =', alpha, 'highest =', highest)
            assert inv_decreasing_tableau(*highest) == t
        else:
            tabs = {(inv_decreasing_tableau(*h), inv_increasing_tableau(*h)) for h in z.get_primed_involution_words()}
            count += len(tabs)
            for t, u in tabs:
                if not exclude_one_row or t.num_rows() > 1:
                    print(t)
                    a = tuple(tuple(_) for _ in reversed(t.get_rows()))
                    r = len(t.get_rows())
                    s = max([1] + list(map(abs, t.row_reading_word())))
                    for rank in [s]:
                        b = (rank - len(a)) * ((),) + a
                        if rank not in invdemazure:
                            invdemazure[rank] = {}

                        t0 = time.time()
                        crystal = AbstractPrimedQCrystal.from_inv_factorization(b, rank, increasing=False)
                        t1 = time.time()
                        print('.  ', 'elapsed', t1 - t0)

                        t0 = time.time()
                        brf = crystal.truncate([f for f in crystal if is_bounded(f)])
                        highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
                        assert len(highest) == 1
                        highest = highest[0]
                        ch = factorization_character(brf)
                        alpha = decompose_q(ch)
                        alpha += (rank - len(alpha)) * (0,)
                        t1 = time.time()
                        print('.. ', 'elapsed', t1 - t0)

                        t0 = time.time()
                        quick_generate_inv_demazure(alpha, invdemazure[rank])
                        t1 = time.time()
                        print('...', 'elapsed', t1 - t0, 'cache size:', len(invdemazure[rank]))

                        print()
                        print('rank =', rank, 'alpha =', alpha, '# bdd elems =', len(brf), 'highest =', highest)

                        assert AbstractPrimedQCrystal.isomorphic_highest_weight_crystals(brf, highest, invdemazure[rank][alpha][0])
                        ans.append((t, u, alpha, highest))
        print()
        print('I_%s' % permutation_size, ': tableaux seen:', count, 'dominant:', dominant)

    INV_DEMAZURE_TABLEAU_CACHE[key] = ans
    return ans


def do_inv_test(n, w, invdemazure=None):
    is_bounded = inv_is_bounded
    invdemazure = invdemazure or {}
    crystal = AbstractPrimedQCrystal.from_involution(w, n, increasing=False)
    for flag in flags(n):
        brf = crystal.truncate([f for f in crystal if is_bounded(f, flag)])
        highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
        for f in highest:
            generate_inv_demazure(crystal.weight(f), invdemazure)
        decomposition, pairs = find_isomorphism(brf, highest, invdemazure)

        ch = factorization_character(brf)
        expected_ch = get_expected_ch(decomposition, lambda alpha: restrict_variables(q_key(alpha), n))  
        print(w, 'flag =', flag, ch == expected_ch, decomposition)
        assert ch == expected_ch and decomposition is not None


def test_inv_demazure_generic(n=2, permutation_size=5):
    invdemazure = {}
    for w in Permutation.involutions(permutation_size):
        do_inv_test(n, w, invdemazure)


def test_inv_demazure(n=2, limit=8):
    is_bounded = inv_is_bounded
    for mu in strict_partitions(n, limit):
        w_mu = Permutation.from_involution_shape(*mu)
        crystal = AbstractPrimedQCrystal.from_involution(w_mu, n, increasing=False)
        brf = crystal.truncate([f for f in crystal if is_bounded(f)])

        highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
        assert len(highest) == 1

        alpha = symmetrize_strict_partition(mu, n)
        demazure = {alpha: brf}
        for ku, i, nu in partition_permutations(alpha, n):
            if i is not None:
                demazure[nu] = crystal.truncate([f for f in crystal if emax(crystal, i, f) in demazure[ku]])
            ch = factorization_character(demazure[nu])
            expected_ch = restrict_variables(q_key(nu), n)
            
            print(nu, w_mu)

            assert inv_negative_one_operator_test(crystal, demazure[nu])
            assert inv_zero_operator_test(crystal, demazure[nu])
            verify_inv_string_lengths(crystal, demazure[nu])

            try:
                if n < max(nu):
                    assert ch != q_key(nu)
                else:
                    assert ch == q_key(nu)
                assert ch == expected_ch
            except:
                # crystal.draw(highlighted_nodes=demazure[nu], extended=crystal.extended_indices)
                input('\n?\n')


def fpf_negative_one_operator_test(crystal, subset):
    if crystal.rank == 1:
        return True
    for i in [1]:
        for b in subset:
            c = crystal.f_operator(-i, b)
            if not (c is None or c in subset):
                return False
            c = crystal.e_operator(-i, b)
            if not (c is None or c in subset):
                return False
    return True
      

def fpf_decreasing_highest(decreasing_tab):
    m = max([0] + [abs(i) for i in decreasing_tab.row_reading_word()]) + 1
    t = Tableau({b: m - v.number for (b, v) in decreasing_tab.mapping.items()})   
    ans = Tableau.inverse_fpf_insertion(t, Tableau({(i, j): i for (i, j) in t.mapping}))
    ans = tuple(tuple(m - x for x in _) for _ in ans)
    return ans


def fpf_increasing_tableau(*word):
    return fpf_insert(*word)[0]


def fpf_decreasing_tableau(*word):
    word = [(w,) if type(w) == int else w for w in word]
    m = max([0] + [i for _ in word for i in _]) + 1
    m += int(m % 2 != 0)
    t, _ = fpf_insert(*(tuple(m - i for i in a) for a in word))
    return Tableau({b: m - v.number for (b, v) in t.mapping.items()})


def test_fpf_demazure_tableau(permutation_size=4, include_dominant=False):
    key = (permutation_size, include_dominant)
    if key in FPF_DEMAZURE_TABLEAU_CACHE:
        return FPF_DEMAZURE_TABLEAU_CACHE[key]

    ans = []
    fpfdemazure = {}
    is_bounded = fpf_is_bounded
    count = 0 
    dominant = 0
    total = len(list(Permutation.fpf_involutions(permutation_size)))
    for z in Permutation.fpf_involutions(permutation_size):
        total -= 1
        print()
        print()
        print()
        print('#', total, ':', 'z =', z.cycle_repr(), 'dominant?', z.is_fpf_dominant())
        if z.is_fpf_dominant():
            count += 1
            dominant += 1
            h = z.get_fpf_involution_word()
            t, u = fpf_decreasing_tableau(*h), fpf_increasing_tableau(*h)
            highest = fpf_decreasing_highest(t)
            mu = t.partition().tuple()
            alpha = skew_symmetric_double(mu)
            if z.fpf_involution_length() > 0 and include_dominant:
                ans.append((t, u, alpha, highest))
            print(t)
            print()
            print('rank = N/A', 'alpha =', alpha, 'highest =', highest)
            assert fpf_decreasing_tableau(*highest) == t
        else:
            tabs = {(fpf_decreasing_tableau(*h), fpf_increasing_tableau(*h)) for h in z.get_fpf_involution_words()}
            count += len(tabs)
            for t, u in tabs:
                print(t)
                a = tuple(tuple(_) for _ in reversed(t.get_rows()))
                r = len(t.get_rows())
                s = max([0] + list(t.row_reading_word())) + 1
                for rank in range(s, s + 1):
                    b = (rank - len(a)) * ((),) + a
                    if rank not in fpfdemazure:
                        fpfdemazure[rank] = {}

                    t0 = time.time()
                    crystal = AbstractQCrystal.from_fpf_factorization(b, rank, increasing=False)
                    t1 = time.time()
                    print('.  ', 'elapsed', t1 - t0)

                    t0 = time.time()
                    brf = crystal.truncate([f for f in crystal if is_bounded(f)])
                    highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
                    assert len(highest) == 1
                    highest = highest[0]
                    ch = factorization_character(brf)
                    alpha = decompose_p(ch)
                    alpha += (rank - len(alpha)) * (0,)
                    t1 = time.time()
                    print('.. ', 'elapsed', t1 - t0)

                    t0 = time.time()
                    quick_generate_fpf_demazure(alpha, fpfdemazure[rank])
                    t1 = time.time()
                    print('...', 'elapsed', t1 - t0, 'cache size:', len(fpfdemazure[rank]))

                    print()
                    print('rank =', rank, 'alpha =', alpha, '# bdd elems =', len(brf), 'highest =', highest)

                    assert AbstractQCrystal.isomorphic_highest_weight_crystals(brf, highest, fpfdemazure[rank][alpha][0])
                    ans.append((t, u, alpha, highest))
        print()
        print('Ifpf_%s' % permutation_size, ': tableaux seen:', count, 'dominant:', dominant)
    
    FPF_DEMAZURE_TABLEAU_CACHE[key] = ans
    return ans


def test_fpf_demazure_generic(n=2, permutation_size=4):
    fpfdemazure = {}
    is_bounded = fpf_is_bounded

    for w in Permutation.fpf_involutions(permutation_size):
        crystal = AbstractQCrystal.from_fpf_involution(w, n, increasing=False)
        print()
        print()
        for flag in [None]: #flags(n):
            brf = crystal.truncate([f for f in crystal if is_bounded(f, flag)])
            highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
            
            for f in highest:
                generate_fpf_demazure(crystal.weight(f), fpfdemazure)
            decomposition, pairs = find_isomorphism(brf, highest, fpfdemazure)

            ch = factorization_character(brf)
            expected_ch = get_expected_ch(decomposition, lambda alpha: restrict_variables(p_key(alpha), n))     
            print(w.oneline_repr(permutation_size), 'flag =', flag, ch == expected_ch, decomposition)
            for h, alpha in pairs:
                tab, _ = fpf_insert(*h)
                print(tab)
                print(alpha)
                print()
            try:
                assert ch == expected_ch and decomposition is not None
            except:
                # crystal.draw(highlighted_nodes=demazure[nu], extended=crystal.extended_indices)
                input('\n?\n')
        if len(decomposition) > 1:
            # crystal.draw(highlighted_nodes=brf, extended=crystal.extended_indices)
            print('* multiple')
            # input('')


def test_fpf_demazure(n=2, limit=8):
    is_bounded = fpf_is_bounded
    for mu in strict_partitions(n, limit):
        w_mu = Permutation.from_fpf_involution_shape(*mu)
        crystal = AbstractQCrystal.from_fpf_involution(w_mu, n, increasing=False)
        brf = crystal.truncate([f for f in crystal if is_bounded(f)])

        highest = [f for f in brf if all(crystal.e_operator(i, f) not in brf for i in crystal.extended_indices)]
        assert len(highest) == 1

        alpha = skew_symmetrize_strict_partition(mu, n)
        demazure = {alpha: brf}

        for ku, i, nu in partition_permutations(alpha, n):
            if i is not None:
                demazure[nu] = crystal.truncate([f for f in crystal if emax(crystal, i, f) in demazure[ku]])
            ch = factorization_character(demazure[nu])
            expected_ch = restrict_variables(p_key(nu), n)
            
            print(nu, w_mu)

            assert fpf_negative_one_operator_test(crystal, demazure[nu])
            verify_fpf_string_lengths(crystal, demazure[nu])

            try:
                if n < max(nu):
                    assert ch != p_key(nu)
                else:
                    assert ch == p_key(nu)
                assert ch == expected_ch
            except:
                # crystal.draw(highlighted_nodes=demazure[nu], extended=crystal.extended_indices)
                input('\n?\n')


def test_brf(n=5):
    # testing whether Demazure crystals are closed under raising operators e_i
    def bounded_ch(c, z):
        ans = 0
        for v in c.vertices:
            if is_bounded(v):
                t = X(0)**0
                for i, a in enumerate(v):
                    t *= X(i + 1)**len(a)
                ans += t
        return ans

    def is_bounded(f, none_bounded=True):
        if f is None:
            return none_bounded
        return all(len(a) == 0 or i < min(a) for i, a in enumerate(f))

    for z in Permutation.all(n):
        # if not z.is_dominant():
        #    continue
        print(z, z.is_dominant())
        for k in range(1, min(6, z.length()) + 2):
            c = AbstractGLCrystal.from_permutation(z, k, False)
            print('  rank =', k, c.extended_indices)
            for v in c.vertices:
                if is_bounded(v) and not c.is_highest_weight(v):
                    if not any(is_bounded(c.e_operator(i, v), False) for i in c.extended_indices):
                        for i in c.extended_indices:
                            print(v, '--', i, '-->', c.e_operator(i, v), is_bounded(c.e_operator(i, v), False))
                        assert False
                    for i in c.extended_indices:
                        if not is_bounded(c.e_operator(i, v)):
                            print('\n?', v, '--', i, '-->', c.e_operator(i, v))
                            assert False
            ch = bounded_ch(c, z)
            print('ch(BRF) =', decompose_into_keys(ch))
            print()


def test_brf_inv(n=5, kk=6):
    # testing whether inv Demazure crystals are closed under raising operators e_i
    def bounded_ch(c, z):
        ans = 0
        for v in c.vertices:
            if is_bounded(v):
                t = X(0)**0
                for i, a in enumerate(v):
                    t *= X(i + 1)**len(a)
                ans += t
        return ans * 2**z.number_two_cycles()

    def is_bounded(f, none_bounded=True):
        if f is None:
            return none_bounded
        return all(len(a) == 0 or i < min(a) for i, a in enumerate(f))

    def e_operator(c, i, v):
        # print('**', i, v)
        if i >= -1:
            return c.e_operator(i, v)

        x = c.s_operator(-i - 1, v)
        if not is_bounded(x):
            return None
        # x = x if is_bounded(x) else v

        y = c.s_operator(-i, x)
        if not is_bounded(y):
            return None
        x = y  # if is_bounded(y) else x

        x = e_operator(c, i + 1, x)
        if x is None:
            return x

        y = c.s_operator(-i, x)
        if not is_bounded(y):
            return None
        # print('*', i, x, '->', y, is_bounded(y))
        x = y  # if is_bounded(y) else x

        y = c.s_operator(-i - 1, x)
        if not is_bounded(y):
            return None
        # print('*', i, x, '->', y, is_bounded(y))
        x = y  # if is_bounded(y) else x

        return x

    for z in Permutation.involutions(n):
        print(z, z.is_dominant())
        for k in range(1, min(kk, z.involution_length()) + 2):
            c = AbstractQCrystal.from_involution(z, k, False)
            print('  rank =', k, c.extended_indices)
            for v in c.vertices:
                if is_bounded(v) and not c.is_highest_weight(v):
                    if all(e_operator(c, i, v) is None for i in c.extended_indices) or any(not is_bounded(c.e_operator(i, v)) for i in c.extended_indices):
                        for i in c.extended_indices:
                            vv = e_operator(c, i, v)
                            print(v, '--', i, '-->', vv, is_bounded(vv), c.e_operator(i, v), is_bounded(c.e_operator(i, v)))
                        # c.draw(False, [f for f in c if is_bounded(f)])
                        # assert False
                        print()
                    # # if True:
                    # if all(c.e_operator(i, v) is None for i in c.indices): # fails even with this condition
                    #     for i in c.extended_indices:
                    #         if not is_bounded(c.e_operator(i, v)):
                    #             print('\n?', v, '--', i, '-->', c.e_operator(i, v))
                    #             # c.draw(True, [f for f in c if is_bounded(f)])
                    #             assert False
            if k < min(kk, z.involution_length()) + 1:
                continue
            ch = bounded_ch(c, z)
            dec = try_to_decompose_q(ch)
            if dec == []:
                print('ch(BRF) =', ch)
            print('ch(BRF) =', dec)
            print()
            assert ch == 0 or dec != []


def test_brf_fpf(n=6):
    # testing whether fpf Demazure crystals are closed under raising operators e_i
    def bounded_ch(c, z):
        ans = 0
        for v in c.vertices:
            if is_bounded(v):
                t = X(0)**0
                for i, a in enumerate(v):
                    t *= X(i + 1)**len(a)
                ans += t
        return ans

    def is_bounded(f, none_bounded=True):
        if f is None:
            return none_bounded
        return all(len(a) == 0 or i < min(a) for i, a in enumerate(f))

    for z in Permutation.fpf_involutions(n):
        print(z, z.is_fpf_dominant())
        for k in range(1, min(6, z.fpf_involution_length()) + 2):
            c = AbstractQCrystal.from_fpf_involution(z, k, False)
            print('  rank =', k, c.extended_indices)
            for v in c.vertices:
                if is_bounded(v) and not c.is_highest_weight(v):
                    if not any(is_bounded(c.e_operator(i, v), False) for i in c.extended_indices):
                        for i in c.extended_indices:
                            print(v, '--', i, '-->', c.e_operator(i, v), is_bounded(c.e_operator(i, v), False))
                        # c.draw(True)
                        assert False
                    # # if True:
                    # if all(c.e_operator(i, v) is None for i in c.indices): # fails even with this condition
                    #     for i in c.extended_indices:
                    #         if not is_bounded(c.e_operator(i, v)):
                    #             print('\n?', v, '--', i, '-->', c.e_operator(i, v))
                    #             # c.draw(True, [f for f in c if is_bounded(f)])
                    #             assert False
            if k < min(6, z.fpf_involution_length()) + 1:
                continue
            ch = bounded_ch(c, z)
            dec = try_to_decompose_p(ch)
            if dec == []:
                print('ch(BRF) =', ch)
            print('ch(BRF) =', dec)
            print()
            assert ch == 0 or dec != []


def test_starb_operator(rank=3, factors=3):
    # fails
    for n in range(2, rank + 1):
        crystal = AbstractPrimedQCrystal
        b = crystal.standard_object(n)
        u = None
        for f in range(1, 1 + factors):
            print(crystal, 'rank =', n, 'factors =', f)
            u = b if u is None else u.tensor(b)
            high = [v for v, _ in u.get_highest_weights()]
            low = [v for v, _ in u.get_lowest_weights()]
            print(high)
            print()
            print(low)
            print()
            for v in high:
                w = u.star_operator(u.starb_operator(v))
                print('high to low:', v, w)
                assert w == u.star_plus_operator(v)
                assert w in low
            for v in low:
                w = u.star_operator(u.starb_operator(v))
                print('low to high:', v, w)
                assert w == u.star_plus_operator(v)
                assert w in high


def test_star_operator(rank=3, factors=3):
    # fails
    for n in range(1, rank + 1):
        crystal = AbstractGLCrystal
        b = crystal.standard_object(n)
        u = None
        for f in range(1, 1 + factors):
            print(crystal, 'rank =', n, 'factors =', f)
            u = b if u is None else u.tensor(b)
            for i in u.indices:
                for x in u:
                    try:
                        y = u.star_operator(u.e_operator(i, u.star_operator(x)))
                        z = u.f_operator(n - i, x)
                        assert y == z
                    except:
                        print('x =', x, 'i =', i, 'w_0 e_i w_0(x) =', y, 'f_{n-i}(x) =', z, '\n')
                        # u.draw()
                        # assert False
                    try:
                        y = u.star_operator(u.f_operator(i, u.star_operator(x)))
                        z = u.e_operator(n - i, x)
                        assert y == z
                    except:
                        print('x =', x, 'i =', i, 'w_0 f_i w_0(x) =', y, 'e_{n-i}(x) =', z, '\n')
                        # u.draw()
                        # assert False


def _test_components(crystal, rank, factors, x):
    b = crystal.standard_object(rank)
    u = None
    uch = 1
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)

        u = b if u is None else u.tensor(b)

        groups = u.group_highest_weights()
        lowgroups = u.group_lowest_weights()
        try:
            for mu in groups:
                rmu = tuple(reversed(mu))
                assert rmu in lowgroups
                assert len(groups[mu]) == len(lowgroups[rmu])
            for mu in lowgroups:
                rmu = tuple(reversed(mu))
                assert rmu in groups
                assert len(groups[rmu]) == len(lowgroups[mu])
        except:
            print('* highest:', groups)
            print('*  lowest:', lowgroups)
            print()
            assert False

        for g in groups.values():
            for i in range(len(g) - 1):
                assert u.isomorphic_subcrystals(g[i], g[i + 1])
        weights = list(groups)
        for i in range(len(weights) - 1):
            mu, nu = weights[i:i + 2]
            assert not u.isomorphic_subcrystals(groups[mu][0], groups[nu][0])

        uch = uch * x
        expected = {tuple(k.mu): v for k, v in uch.items() if len(k.mu) <= rank}
        actual = {Partition.trim(k): len(v) for k, v in groups.items()}
        try:
            assert expected == actual
        except:
            print('expected =', expected)
            print('  actual =', actual)
            print()
            assert False


def test_gl_components(rank=3, factors=3):
    x = Schur(1)
    for rank in range(1, rank + 1):
        _test_components(AbstractGLCrystal, rank, factors, x)


def test_q_components(rank=3, factors=3):
    x = SchurP(1)
    for rank in range(2, rank + 1):
        _test_components(AbstractQCrystal, rank, factors, x)


def test_primed_q_components(rank=3, factors=3):
    x = SchurQ(1)
    for rank in range(2, rank + 1):
        _test_components(AbstractPrimedQCrystal, rank, factors, x)


def flatten(t):
    if type(t) != tuple:
        return (t,)
    if len(t) == 0:
        return ()
    return (flatten(t[0]) if type(t[0]) == tuple else (t[0],)) + flatten(t[1:])


def test_forgetul_functor(rank=3, factors=3):
    # fails
    crystal = AbstractQCrystal
    b = crystal.standard_object(rank)
    u = None
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)
        u = b if u is None else u.tensor(b)
        for vertex in u:
            flat = flatten(vertex)
            for i in [2]:
                for j in range(1, len(flat)):
                    target = u.f_operator(-i, vertex)
                    w1, w2 = flat[:j], flat[j:]
                    if not any(a in [i, i + 1] for a in w1):
                        w3 = AbstractQCrystal.f_operator_on_words(-i, w2)
                        expected = (w1 + w3) if w3 is not None else None
                    else:
                        w3 = AbstractQCrystal.f_operator_on_words(-i, w1)
                        expected = (w3 + w2) if w3 is not None else None
                    print('f:', 'i =', -i, w1, ':', w2, '-->', expected, 'vs', target)
                    if target:
                        assert flatten(target) == expected
                    else:
                        assert expected is None

                    target = u.e_operator(-i, vertex)
                    if not any(a in [i, i + 1] for a in w1):
                        w3 = AbstractQCrystal.e_operator_on_words(-i, w2)
                        expected = (w1 + w3) if w3 is not None else None
                    else:
                        w3 = AbstractQCrystal.e_operator_on_words(-i, w1)
                        expected = (w3 + w2) if w3 is not None else None
                    print('e:', 'i =', -i, w1, ':', w2, '-->', expected, 'vs', target)
                    if target:
                        assert flatten(target) == expected
                    else:
                        assert expected is None


def _test_operators_on_words(crystal, rank, factors):
    b = crystal.standard_object(rank)
    u = None
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)
        u = b if u is None else u.tensor(b)
        for vertex in u:
            flat = flatten(vertex)
            for i in u.indices:
                target = u.f_operator(i, vertex)
                target_word = crystal.f_operator_on_words(i, flat)
                if target:
                    assert flatten(target) == target_word
                else:
                    assert target_word is None

                target = u.e_operator(i, vertex)
                target_word = crystal.e_operator_on_words(i, flat)
                if target:
                    assert flatten(target) == target_word
                else:
                    assert target_word is None


def test_gl_operators_on_words(rank=3, factors=3):
    u = Word(1, 3, 3, 1, 2, 1, 2)
    v = Word(2, 3, 3, 1, 2, 1, 2)
    w = Word(3, 3, 3, 1, 2, 1, 2)

    assert AbstractGLCrystal.f_operator_on_words(2, u) is None
    assert AbstractGLCrystal.f_operator_on_words(2, v) == w
    assert AbstractGLCrystal.e_operator_on_words(2, w) == v

    for rank in range(1, rank + 1):
        _test_operators_on_words(AbstractGLCrystal, rank, factors)


def test_q_operators_on_words(rank=4, factors=4):
    for rank in range(2, rank + 1):
        _test_operators_on_words(AbstractQCrystal, rank, factors)


def test_primed_q_operators_on_words(rank=4, factors=4):
    for rank in range(2, rank + 1):
        _test_operators_on_words(AbstractPrimedQCrystal, rank, factors)


def all_primed_words(max_letter, length):
    for sgn in range(2**length):
        for v in range(max_letter**length):
            a = []
            for _ in range(length):
                a.append((-1)**(sgn % 2) * ((v % max_letter) + 1))
                sgn = sgn >> 1
                v = v >> 1
            yield Word(a)


def random_primed_words(count):
    def generator(max_letter, length):
        for _ in range(count):
            a = []
            for _ in range(length):
                a.append((-1)**(random.randint(0, 1)) * random.randint(1, max_letter))
            yield Word(a)
    return generator


def _test_operators_on_shifted_tableaux(generator, rank, length, verbose=False):
    for word in generator(rank, length):
        tab = word.mixed_insert()[0]
        for i in range(-1, rank):
            target_word = AbstractPrimedQCrystal.f_operator_on_words(i, word)
            target_tab = tab.shifted_crystal_f(i)
            backtrack_tab = target_tab.shifted_crystal_e(i) if target_tab else None
            try:
                if target_word is None:
                    assert target_tab is None
                else:
                    assert target_tab == target_word.mixed_insert()[0]
                    assert backtrack_tab == tab
                if verbose:
                    print(word, '--', i, '-->', target_word)
                    print(tab)
                    print(target_tab)
                    print(backtrack_tab)
                    print()
            except:
                print(word, '--', i, '-->', target_word)
                print(tab)
                print(target_tab)
                print(backtrack_tab)
                print()
                assert False

            target_word = AbstractPrimedQCrystal.e_operator_on_words(i, word)
            target_tab = tab.shifted_crystal_e(i)
            backtrack_tab = target_tab.shifted_crystal_f(i) if target_tab else None
            try:
                if target_word is None:
                    assert target_tab is None
                else:
                    assert target_tab == target_word.mixed_insert()[0]
                    assert backtrack_tab == tab
                if verbose:
                    print(word, '<--', i, '--', target_word)
                    print(tab)
                    print(target_tab)
                    print(backtrack_tab)
                    print()
            except:
                print(word, '<--', i, '--', target_word)
                print(tab)
                print(target_tab)
                print(backtrack_tab)
                print()
                assert False


def test_primed_q_operators_on_tableaux(rank=4, length=4):
    _test_operators_on_shifted_tableaux(all_primed_words, rank, length)


def test_random_primed_q_operators_on_tableaux(rank=10, length=10):
    _test_operators_on_shifted_tableaux(random_primed_words(100), rank, length, True)
