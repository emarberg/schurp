from symmetric import *
from partitions import *


# cache = {}
# for n in range(13):
#     print(' . . .', n)
#     for mu, nu in StrictPartition.skew_pairs(n):
#         cache[(mu, nu)] = SchurP(mu).skew(nu)


# icache = {}
# for n in range(13):
#     print(' . . .', n)
#     for mu, nu in StrictPartition.skew_pairs(n):
#         icache[(mu, nu)] = SchurQ(mu).skew(nu)

from signed import *

def count_modified_type_c_inv_words(n):
    z = SignedPermutation.longest_element(n)
    a = 0
    for w in z.get_atoms():
        a += 2**w.ell_zero() * len(w.get_reduced_words())
    return a


from permutations import *

def count_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    return Permutation(*oneline).count_pipe_dreams()


def count_involution_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    return Permutation(*oneline).count_involution_pipe_dreams()

def count_fpf_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    w = Permutation(*oneline)
    for i in range(1, shift, 2):
        w *= Permutation.s_i(i)
    return w.count_fpf_involution_pipe_dreams()


def skew_shape(w, transpose=False, inv=False):
    d = w.fpf_rothe_diagram() if not inv else w.involution_rothe_diagram()
    if transpose:
        d = [(j, i) for (i, j) in d]
    k = sorted({i for i, j in d})
    l = len(k)
    m = max(k) if k else 0
    c = [len({j for i, j in d if i == index}) for index in range(1, m + 1)]
    shape = [
        (i, j + (1 + i - k[i - 1] - c[k[i - 1] - 1]))
        for i in range(1, l + 1)
        for j in range(c[k[i - 1] - 1])
    ]
    offset = (1 - min([j for i, j in shape])) if shape else 0
    shape = [(i, j + offset) for i, j in shape]
    return shape


def count_se_property(n):
    c = 0
    for w in Permutation.fpf_involutions(n):
        d = w.fpf_rothe_diagram()
        if all((j, k) in d for (i, k) in d for (j, l) in d if i < j and k > l):
            c += 1
    return c


def test(n):
    for w in Permutation.fpf_involutions(n):
        if w.fpf_involution_length() > 4:
            continue
        f = FPFStanleyExpander(w).expand()
        pairs = [pr for pr, g in cache.items() if f == g]
        if len(pairs) == 0:
            # if has_se_property(w):
            #     print('z =', w)
            #     print('fpf code =', w.fpf_involution_code())
            #     w.print_fpf_rothe_diagram()
            #     input('?')
            continue
        if not has_se_property(w):
            continue
        modified = pairs #[(mu, nu) for (mu, nu) in pairs if mu != nu and len(nu) > 0]
        if len(f) == 1:
            continue
        print('z =', w)
        print('fpf code =', w.fpf_involution_code())
        w.print_fpf_rothe_diagram(sep='-')
        print('\n.\n')
        FPFStanleyExpander(w).dearc().print_involution_rothe_diagram(sep='-')
        print('\n.\n')
        w.get_min_fpf_atom().inverse().print_rothe_diagram(sep='-')
        print('\n.\n')
        print()
        sh = skew_shape(w, True)
        print(w.print_diagram(sh, sep='-'))
        print()
        print(f)
        for mu, nu in modified:
            print(mu.shape - nu.shape)
            print('\n', mu, nu, '\n')
        print()
        print()
        print()


def itest(n):
    for w in Permutation.involutions(n):
        f = SchurP.to_q_basis(InvStanleyExpander(w).expand() * 2**w.number_two_cycles())
        pairs = [pr for pr, g in icache.items() if f == g]
        if len(pairs) == 0:
            continue
        if not has_se_property(w, inv=True):
            continue
        modified = pairs
        if len(modified) > 8:
            continue
        if len(f) == 1:
            continue
        print('z =', w)
        print('inv code =', w.involution_code())
        w.print_involution_rothe_diagram()
        #print('\n.\n')
        #sh = Shape(skew_shape(w, True, inv=True))
        #print(sh)
        print()
        print(f)
        for mu, nu in modified[:]:
            print(mu.shape - nu.shape)
            print('\n', mu, nu, '\n')
        print()
        print()
        print()


def has_se_property(w, inv=False):
    d = w.fpf_rothe_diagram() if not inv else w.involution_rothe_diagram()
    return all((j, k) in d for (i, k) in d for (j, l) in d if i < j and k > l)


def get_se_property(n):
    for w in Permutation.fpf_involutions(n):
        d = w.fpf_rothe_diagram()
        if all((j, k) in d for (i, k) in d for (j, l) in d if i < j and k > l):
            yield w


def Ffpf(w):
    return FPFStanleyExpander(w).expand()


def test_stanley(n):
    for w in Permutation.fpf_involutions(n):
        f = SchurP.to_schur_basis(FPFStanleyExpander(w).expand())
        g = Vector()
        for x in w.get_fpf_atoms():
            g += StanleyExpander(x).expand()
        assert f == g


def P(mu, nu):
    ans = SchurQ(*mu).skew(StrictPartition(*nu))
    ans = SchurQ.to_p_basis(ans) / 2**(len(mu) - len(nu))
    return ans


def print_se_property(n):
    for w in get_se_property(n):
        print(w)
        print(w.fpf_involution_code())
        print(w.print_diagram(skew_shape(w)))
        w.print_fpf_rothe_diagram()
        print()
        print(Ffpf(w))
        print()
        print()
        print()



from schubert import *

ANTIPODES = {}

def antipode(n):
    assert n >= 0
    if n not in ANTIPODES:
        if n == 0:
            ans = MPolynomial.one()
        else:
            a = antipode(n - 1)
            ans = a * X(0) * (n - 1)
            for i in range(a.total_degree() + 1):
                c = a.coefficient(1, i)
                f = (i + 1) * X(1)**(i + 1) + i * X(0) * X(1)**i
                ans += f * c
            ans *= -1
        ANTIPODES[n] = ans
    return ANTIPODES[n]




from schubert import *
import schubert


def test_formal_operators_experiments():
    x = lambda i: MPolynomial.monomial(i) # noqa
    A = lambda i: Operator.create(i) * x(i) * (1 - x(i + 1)) # noqa
    a = lambda i: Operator.create(i) * (1 - x(i + 1)) # noqa
    print(A(2) * A(1) - a(2) * a(1) * x(1)**2)
    print()
    print(a(1) * x(1))


def test_grothendieck_transitions(n):
    def terms(w, j):
        queue = [(w, n + 1)]
        while queue:
            y, k = queue[0]
            queue = queue[1:]

            if k <= j:
                continue

            s = Permutation.transposition(j, k)
            z = y * s
            if z.length() == y.length() + 1:
                yield z
                queue.append((z, k - 1))
            queue.append((y, k - 1))

    g = list(Permutation.all(n))
    for w in g:
        for i in range(1, n + 1):
            var = 1 - schubert.x(i)

            ts = []
            for k in range(1, i):
                t = Permutation.cycle([k, i])
                v = w * t
                if v.length() == w.length() + 1:
                    ts.append(k)

            ttt = [(w, 1)]
            for k in ts:
                t = Permutation.cycle([k, i])
                ttt += [(v * t, -a) for v, a in ttt]

            f = 0
            for v, a in ttt:
                f += Grothendieck.get(v) * a
            f = f * var

            sp = ''.join(['(1 - t_{%s,%s})' % (k, i) for k in ts]) if ts else '1'
            print('G_%s * %s * (%s) = ' % (w, sp, var))
            print()
            try:
                dec = Grothendieck.decompose(f)
                print()
                print('    ', dec)
                print()

                a = Grothendieck.get(w)
                for z in terms(w, i):
                    if (z.length() - w.length()) % 2 == 0:
                        a += Grothendieck.get(z)
                    else:
                        a -= Grothendieck.get(z)
                assert f == a
            except:
                print('     halted computation')
            print()
            print()


def test_fpf_transitions(n):
    def terms(w, j):
        queue = [(w, n + 1)]
        while queue:
            y, k = queue[0]
            queue = queue[1:]

            if k <= j:
                continue

            s = Permutation.transposition(j, k)
            z = s * y * s
            if z.fpf_involution_length() == y.fpf_involution_length() + 1:
                yield z
                queue.append((z, k - 1))
            queue.append((y, k - 1))

    g = list(Permutation.fpf_involutions(n))
    for w in g:
        cyc = [
            (i, j) for i, j in w.cycles
            # if not any(k < i and l < j for k, l in w. cycles)
        ]
        w = w * Permutation.s_i(n + 1)
        for i, j in cyc:
            var = 1 - schubert.x(i) - schubert.x(j) + schubert.x(i) * schubert.x(j)
            ts = []
            for k in range(1, i):
                t = Permutation.cycle([k, i])
                v = t * w * t
                if v.fpf_involution_length() == w.fpf_involution_length() + 1:
                    ts.append(k)
            ttt = [(w, 1)]
            for k in ts:
                t = Permutation.cycle([k, i])
                ttt += [(t * v * t, -a) for v, a in ttt]
            f = 0
            for v, a in ttt:
                f += FPFGrothendieck.get(v) * a
            # ttt = Vector({w.fpf_trim(): a for w, a in ttt})
            f = f * var
            sp = ''.join(['(1 - t_{%s,%s})' % (k, i) for k in ts]) if ts else '1'
            print('G_%s * %s * (%s) = ' % (w, sp, var))
            print()
            try:
                dec = FPFGrothendieck.decompose(f)
                print()
                print('    ', dec)
                print()
                a = FPFGrothendieck.get(w)
                # print('    ', list(terms(w, j)))
                for z in terms(w, j):
                    if (z.fpf_involution_length() - w.fpf_involution_length()) % 2 == 0:
                        a += FPFGrothendieck.get(z)
                    else:
                        a -= FPFGrothendieck.get(z)
                # print('    ', a - f)
                assert f == a
            except:
                print('     halted computation')
            print()
            print()



def ftest(n):
    g = list(Permutation.fpf_involutions(n))
    for v in g[1:]:
        c = max(v.get_fpf_visible_inversions())
        t = Permutation.cycle(c)
        u = t * v * t
        s = set(v.fpf_rothe_diagram()) - set(u.fpf_rothe_diagram())
        #
        #
        assert len(s) == 1
        i, j = s.pop()
        k, l = c
        assert i == k
        assert list(sorted([i, j])) in u.cycles
        #
        ts = []
        for h in range(1, k):
            t = Permutation.cycle((h, k))
            if (t * u * t).fpf_involution_length() == u.fpf_involution_length() + 1:
                ts.append(t)
        #
        terms = [(u, 1)]
        for t in ts:
            terms = terms + [(t * w * t, -a) for w, a in terms]
        terms = Vector({w.fpf_trim(): a for w, a in terms})
        #
        f = FPFGrothendieck.get(v)
        g = FPFGrothendieck.get(u)
        d = g - f
        #
        # print(f)
        # print()
        # print(g)
        # print()
        # print(d)
        # print()
        e = d.divide_linear(i, -1).divide_linear(j, -1)
        # print(e)
        # print()
        decomp = FPFGrothendieck.decompose(e)
        print('u =', u)
        print('v =', v)
        print()
        print('i, j =', (i, j))
        print('k, l =', c)
        v.print_fpf_rothe_diagram()
        u.print_fpf_rothe_diagram()
        print()
        print(ts)
        print()
        print('G_v = G_u - (1 - x_%i) (1 - x_%i) [' % (i, j), decomp, ']')
        print('                                 ', terms)
        assert decomp == terms
        print()
        print()
        print()
        print()
        print()


def dtest(n):
    g = list(Permutation.all(n))
    var = [x(i) for i in range(1, n)]
    for w in g:
        for i, v in enumerate(var):
            print('G_%s * x_%s' % (w, i + 1))
            f = Grothendieck.get(w)
            print(Grothendieck.decompose(f * v))
            print()


def fpftest(n):
    clss = FPFGrothendieck
    g = list(Permutation.fpf_involutions(n))
    for w in g:
        var = [
            x(i) + x(j) - x(i) * x(j) for i, j in w.cycles
            if not any(k < i and l < j for k, l in w. cycles)
        ]
        for v in var:
            print('G_%s * (%s)' % (w, v))
            f = clss.get(w)
            try:
                print(clss.decompose(f * v))
            except RecursionError:
                print('* Recursion error')
            print()




def generate(n):
    # generates FPF "symplectic" Hecke atoms"
    assert n % 2 == 0
    #
    def next(w):
        for i in range(0, len(w) - 3, 2):
            a, d, b, c = w[i: i + 4]
            if a < b < c < d:
                yield w[:i] + (b, c, a, d) + w[i + 4:]
                yield w[:i] + (b, d, a, c) + w[i + 4:]
            b, c, a, d = w[i:i + 4]
            if a < b < c < d:
                yield w[:i] + (a, d, b, c) + w[i + 4:]
                yield w[:i] + (b, d, a, c) + w[i + 4:]
            b, d, a, c = w[i:i + 4]
            if a < b < c < d:
                yield w[:i] + (a, d, b, c) + w[i + 4:]
                yield w[:i] + (b, c, a, d) + w[i + 4:]
    #
    start = []
    for i in range(n // 2):
        start += [i + 1, n - i]
    start = tuple(start)
    #
    seen = set()
    level = {start}
    while level:
        next_level = set()
        for w in level:
            seen.add(w)
            next_level |= set(next(w))
        level = next_level - seen
    return seen


from words import *


def braids(word):
    ans = 0
    for i in range(len(word) - 2):
        if len({word[i], word[i + 1], word[i + 2]}) == 2:
            ans += 1
    return ans


def test(n):
    w = tuple(reversed(range(1, n + 1)))
    words = get_fpf_involution_words(w)
    a = 0
    for u in words:
        if len(u) >= 2 and u[1] % 2 != 0:
            a += 1
        a += 2 * braids(u)
    return (a, len(words))




from symmetric import *;
w=Permutation([4,5,7,1,2,6,3])
i=InvStanleyExpander(w)
i.expand()


from signed import *
w = SignedPermutation(-1,-2,-3,-4)
w.inv_stanley_schur_q_decomposition()


from symmetric import *
mu = StrictPartition(5,3,1)
nu = StrictPartition(2,1)
SchurQ(mu).skew(nu)

from signed import *
SignedPermutation(-1,-2,-3).inv_stanley_schur_q_decomposition()


from symmetric import *
from signed import *

def test(n):
    w = SignedPermutation.longest_element(n)
    mu = StrictPartition(*range(2 * n - 1,0,-2))
    nu = StrictPartition(*range(n - 1,0,-1))
    a = w.inv_stanley_schur_q_decomposition()
    b = SchurQ(mu).skew(nu)
    print('')
    print(mu.shape)
    print(a)
    print('')
    print(nu.shape)
    print(b)
    print('')
    print(a - b)


from collections import defaultdict
from symmetric import *
from signed import *
from partitions import *


def skewdict(n, d=None):
    if d is None:
        d = defaultdict(dict)
    for k in range(n + 1):
        for mu, nu in StrictPartition.skew_pairs(k):
            key = abs(mu) - abs(nu)
            subkey = mu.shape - nu.shape
            if subkey in d[key]:
                continue
            d[key][subkey] = SchurQ(mu).skew(nu)
    return d


def altskewdict(n, d=None):
    if d is None:
        d = defaultdict(dict)
    for p in Partition.all(n):
        k = len(p) - 1
        nu = StrictPartition(*list(range(k, 0, -1)))
        mu = nu + p
        key = abs(mu) - abs(nu)
        subkey = (mu.shape - nu.shape).justify()
        if subkey in d[key]:
            continue
        d[key][subkey] = SchurQ(mu).skew(nu)
    return d


def involdict(n, d=None):
    if d is None:
        d = defaultdict(dict)
    for k in range(n + 1):
        for w in SignedPermutation.involutions(k):
            key = w.involution_length()
            if w in d[key]:
                continue
            v = w.reduce()
            if v in d[key]:
                d[key][w] = d[key][v]
            else:
                d[key][w] = w.inv_stanley_schur_q_decomposition()
    return d


def test(n, invol=None, skew=None):
    invol = involdict(n, invol)
    for length in sorted(invol.keys()):
        print('\nlength %i:\n' % length)
        skew = altskewdict(length, skew)
        for w in invol[length]:
            if w != w.reduce():
                continue
            ans = set()
            for shape in skew[length]:
                if invol[length][w] == skew[length][shape]:
                    ans |= {str(shape.justify())}
            if ans:
                print(('\u011c_%s = \n\n' % w) + '\n\n'.join(sorted(ans)) + '\n')
            else:
                print(('\u011c_%s = ' % w) + 'NO SKEW SCHUR Q-FUNCTION\n')
    return invol, skew




# number of w is B_n with iF_w = single S-function:
#   n 0 1 2 3  4  5  6
#   # 1 2 6 14 30 55 93

# number of w is B_n with iF_w = single S-function:
#   n 0 1 2 3  4  5   6
#   # 1 2 6 18 59 183 547

from signed import *


def s_test(n):
    bns = 0
    g = list(SignedPermutation.involutions(n))
    ans = []
    for w in g:
        v1 = w.inv_stanley_schur_s_decomposition()
        v2 = w.inv_stanley_schur_q_decomposition()
        print('\n\u011c_%s = %s' % (w, v1))
        print((n + 2) * ' ' + ' = %s\n' % (v2))
        if len(v1) == 1:
            ans.append((w, v1))
        if sym_at(w):
            bns += 1
            print('True')
        else:
            print('False')
    #print('\n%i\n' % len(ans))
    return bns


from signed import *
from tableaux import *
# n = 3
# w = SignedPermutation.longest_element(n)
# words = list(w.get_involution_words())


def tableau(word, n):
    parts = []
    for i in range(n):
        parts += [word[:i + 1]]
        word = word[i + 1:]
    parts = reversed(parts)
    string = ';'.join([','.join([str(i) for i in part]) for part in parts])
    return Tableau.from_string(string)




def print_i(n):
    d = defaultdict(list)
    for w in SignedPermutation.involutions(n):
        d[w.involution_length()].append(w)
    for i in d:
        print('------------------------------------')
        for w in d[i]:
            if w(n) == n:
                continue
            print(w)
            print(w.inv_stanley_schur_q_decomposition())
            print('')

def sym_at(w):
    a = w.get_atoms()
    return all(u.inverse() in a for u in a)

def print_s(n):
    for i in range(1, n + 1):
        print('------------------------------------')
        for mu in Partition.all(i):
            print(mu)
            print(SchurQ.s_lambda(mu))
            print('')


from signed import *

def atoms_view(n):
    def appear(v):
        return ' '.join(list(reversed(['%s:%s' % (str(k.mu), str(val)) for k, val in v.items()])))
    #
    for w in SignedPermutation.involutions(n):
        print('------------------------------------')
        least = None
        for a in w.get_atoms():
            v = a.stanley_schur_q_decomposition()
            if least is None:
                least = min(v.keys())
            print('F^B_%s = %s' % (a.inverse(), appear(v)))
            # print((n + 4) * ' ' + ' = %s' % v)
            # print((n + 4) * ' ' + ' = %s\n' % appear)
            print('')
        v = w.inv_stanley_schur_q_decomposition()
        print('G^B_%s = %s' % (w, appear(v)))
        print('')
        if min(v.keys()) != least:
            input('??')



from permutations import *
from words import *



def test_eg_mirror(n):
    special = {}
    for w in Permutation.involutions(n):
        special[w] = 1
        print(w, '=', ''.join([str(i) for i in w.oneline]))
        print('------------------------------------')
        for e in w.get_involution_words():
            f = Permutation.double_involution_word(e)
            g = tuple(reversed(e)) + tuple(e)
            print(e, '==>', f)
            print('')
            p, q = Word(*e).involution_insert(verbose=False)
            a, b = Word(*f).eg_insert()
            x, y = Word(*g).hecke_insert()
            print(p)
            print('')
            print(a)
            print('')
            print(x)
            print('')
            print('')
            # print(q)
            # print('')
            # print(b)
            # print('')
            # print('')
            if x != p.double() and w in special:
                del special[w]
                input('?')
    return list(special.keys())



from permutations import *
from signed import *
from words import *


def correct(model, given, word):
    def find(x, word):
        return [i for i in range(len(word)) if word[i] == x][0]

    n = len(model)
    missing = sorted([i for i in range(1, n + 1) if i not in given.entries()])
    m = len(missing)
    descents = set()
    for x in missing:
        j = find(x, word)
        if j < len(word) - 1 and word[j + 1] < word[j]:
            descents.add(find(x, missing) + 1)


def test_ab(n):
    w = Permutation.longest_element(n + 1)
    z = SignedPermutation.longest_element(n)
    inv = list(z.get_involution_words())
    red = list(w.get_reduced_words())
    invbijection = {}
    bijection = {}
    t = Word(*tuple(n - i for i in inv[0])).eg_insert()[0]
    for e in inv:
        f = tuple(n - i for i in e)
        a, b = Word(*f).eg_insert() if f in red else (None, None)
        x, y = Word(*f).hecke_insert()
        sh = t.shape() - x.shape()
        assert sh == sh.corners()
        if f not in red:
            print(f in red, e, '==>', f)
            print(x)
            print('')
            if a is not None and a != x:
                print(a)
                print('')
            print(y)
            print('')
            print('')
        # if f not in red:
        #     x, y = correct(model=t, given=y, word=f)
        # assert x == t
        # assert y not in dictionary
        # invbijection[y] = f
        # bijection[f] = y
    return bijection, invbijection



from signed import *


def test_transitions(n):
    w = SignedPermutation(*list(range(n + 1, 0, -1)))
    a = w.queue_stanley_decomposition(n)
    a = sorted(list(a), key=str)
    b = sorted(list(SignedPermutation.longest_element(n).get_atoms()), key=str)






from words import *

def test(n, verbose=False):
    for w in Word.permutations(n):
        p, q = w.eg_insert()
        pp, qq = w.mystery_insert(verbose=verbose)
        if pp.is_increasing():
            continue
        print('* w =', w)
        print()
        print(p)
        print()
        print(q)
        print()
        print()
        print(pp)
        print()
        print(qq)
        print()
        print()



from signed import *
from collections import defaultdict

def coincidences_d(n):
    stan = {w: w.stanley_schur_d_decomposition() for w in SignedPermutation.permutations(n)}
    inv = {w: w.inv_stanley_schur_d_decomposition() for w in SignedPermutation.involutions(n) if w.is_even_signed()}
    co = defaultdict(set)
    for w, f in inv.items():
        for x, g in stan.items():
            if f == g:
                co[w.reduce()].add(x.reduce())
    return co

# co = coincidences_d(5)
# for w, l in co.items():
#      print(w,':',l)


def coincidences(n):
    stan = {w: w.stanley_schur_q_decomposition() for w in SignedPermutation.permutations(n)}
    inv = {w: w.inv_stanley_schur_q_decomposition() for w in SignedPermutation.involutions(n)}
    co = defaultdict(set)
    for w, f in inv.items():
        for x, g in stan.items():
            if f == g:
                co[w.reduce()].add(x.reduce())
    return co




from signed import *
# n = 3; w = EvenSignedPermutation.longest_element(n); w.get_atoms()


# len(list(w.get_flattened_involution_words()))


from permutations import *
from words import *

def invword(n):
    for w in Permutation.involutions(n):
        for word in w.get_involution_words():
            double = tuple(reversed(word)) + word
            p = Word(*word).involution_insert()[0]
            q = Word(*double).hecke_insert()[0]
            print(p, '\n')
            print(q, '\n\n')
            assert p.double() == q or (2 * len(q) - q.count_diagonal_cells()) != len(p)
            yield word, double, p.double() == q

def fpfword(n):
    for w in Permutation.fpf_involutions(n):
        for word in get_fpf_involution_words(tuple(w.oneline)):
            double = tuple(reversed(word)) + word
            print(Word(*word).fpf_insert()[0], '\n')
            print(Word(*double).hecke_insert()[0], '\n\n')
            yield word, double


# ans = [(a, b) for (a, b, c) in invword(n) if not c]
# for a, b in ans:
#     p = Word(*a).involution_insert()[0]
#     q = Word(*b).modified_hecke_insert(verbose=False)[0]
#     r = Word(*b).hecke_insert()[0]
#     print(p, '\n')
#     print(r, '\n')
#     print(q, '\n\n')
#     print(a, b)


from permutations import *
from words import *


def inv_subtest(e):
    m = max(e) if e else 0
    p, q = Word(*e).involution_insert()
    w = tuple(Permutation.double_involution_word(e))
    a, b = Word(*reversed(w)).eg_insert()
    r = a.lower_half()
    rr = a.upper_half().transpose()
    s = p.lower_half()
    # t = Tableau()
    # for i, j in s:
    #     t = t.set(i, j, s.entry(i, j) - r.entry(i, j))
    if all(i - 1 not in e and i + 1 not in e for i in e):
        return
    if s != r and s != rr:
        print('word =', ', '.join(map(str, e)))
        print(' ext =', ', '.join(map(str, w)))
        print()
        print(s)
        print()
        print(a)
        print()
        input('?\n')


def inv_test(n):
    g = Permutation.involutions(n)
    for w in g:
        print(w.cycle_repr())
        print()
        for e in get_involution_words(tuple(w.oneline)):
            inv_subtest(e)
        print()


def fpf_test(n):
    assert n % 2 == 0
    g = Permutation.fpf_involutions(n)
    for w in g:
        print(w.cycle_repr())
        print()
        for e in get_fpf_involution_words(tuple(w.oneline)):
            m = max(e) if e else 0
            mid = tuple(i for i in range(m + 2) if i % 2 != 0 and (i in e or i - 1 in e or i + 1 in e))
            p, q = Word(*e).fpf_insert()
            x, y = (Word(*mid) | Word(*e)).involution_insert()
            z = x
            for i in range(1, len(mid) + 1):
                l = len(x.get_column(i))
                _, z = z.pop(l, i)
            if p == z.translate_left():
                continue
            print(p)
            print()
            print(x)
            print()
            input('\n\n')
        print()


def subtest(e):
    m = max(e) if e else 0
    p, q = Word(*e).fpf_insert()
    pre = tuple(i for i in range(m + 2) if i % 2 != 0 and not (i in e or i - 1 in e or i + 1 in e))
    mid = tuple(i for i in range(m + 2) if i % 2 != 0 and (i in e or i - 1 in e or i + 1 in e))
    doubled = Word(*reversed(e)) | Word(*reversed(mid)) | Word(*e)
    a, b = doubled.eg_insert()
    x, y = (Word(*mid) | Word(*e)).involution_insert()
    r = a.strict_lower_half().translate_left()
    assert x == a.transpose().lower_half()
    t = Tableau()
    for i, j in p:
        t = t.set(i, j, p.entry(i, j) - r.entry(i, j))
    #
    # correct r
    seen = []
    corrections = 0
    for i, j in sorted(r, key=lambda ij: (-ij[0], ij[1])):
        v = r.entry(i, j).number
        if v % 2 != 0 and not (v - 1 in seen or v in seen or v + 1 in seen):
            r = r.set(i, j, MarkedNumber(v + 1))
            corrections += 1
        seen += [r.entry(i, j).number]
    if corrections == 0 and a != a.transpose():
        print(tuple(reversed(e)), tuple(reversed(mid)), e)
        print()
        print(a)
        print()
        print(p)
        print()
        print(x.transpose())
        print()
        print()
        print()
    #
    #
    print(e, '->', doubled)
    print()
    print(p)
    print()
    print(q)
    print()
    print(a)
    print()
    print(b)
    print()
    input('?')
    assert r == p
    # if r != p:
    #     print(pre)
    #     print(mid)
    #     print('word =', ' '.join(map(str, e)))
    #     print()
    #     print(p)
    #     print()
    #     print(r)
    #     print()
    #     # for i in range(len(e)):
    #     #     ee = e[:i + 1]
    #     #     a, b = (Word(*reversed(ee)) | Word(*reversed(mid)) | Word(*ee)).eg_insert()
    #     #     print(a)
    #     #     print()
    #     print(t)
    #     print()
    #     input('?\n')


def test_hecke_doubling_random(n):
    e = Permutation.random_fpf_involution_word(n)
    for i in range(n):
        print(i, 'of', n)
        subtest(e[:i])
        print()


def test_hecke_doubling(n):
    assert n % 2 == 0
    g = Permutation.fpf_involutions(n)
    for w in g:
        print(w.cycle_repr())
        print()
        for e in get_fpf_involution_words(tuple(w.oneline)):
            subtest(e)
        print()



from permutations import *
from words import *

def test(n):
    g = Permutation.involutions(n)
    for w in g:
        print(w)
        for e in get_involution_words(tuple(w.oneline)):
            ee = tuple(reversed(e)) + e
            print(''.join(map(str, e)), ':', ''.join(map(str, ee)))
            p, q = Word(*e).involution_insert()
            pp, qq = Word(*ee).hecke_insert()
            print()
            print(p)
            print()
            print(pp)
            print('\n')
            if p.double() != pp:
                print(p.double())
                input('?')
        print()


from permutations import *
from words import *

def test(n):
    def compact(e):
        return ''.join(map(str, e))
    g = Permutation.fpf_involutions(n)
    for w in g:
        print(w)
        for e in get_fpf_involution_words(tuple(w.oneline)):
            if len(e) == 0:
                continue
            a, b = min(e), max(e)
            if any(i not in e and i + 1 not in e for i in range(a, b)):
                continue
            p, q = Word(*e).fpf_insert()
            d = p.fpf_double()
            middle = tuple(i for i in range(min(e) - 1, max(e) + 2) if i % 2 != 0)
            # middle += (max(middle) + 2,)
            ee = tuple(reversed(e)) + middle + e
            pp, qq = Word(*ee).eg_insert()
            # if d != pp:
            #     m = MarkedNumber(middle[-1])
            #     found = set(pp.find(m).mapping.keys()) - set(d.find(m).mapping.keys())
            #     if len(found) != 1:
            #         print(m, 'found:\n', found)
            #         input('??')
            #     for (i, j) in found:
            #         d = d.set(i, j, m)
            ind = [i for i in range(pp.max_row + 1) if (i, i) in pp.mapping and pp.entry(i, i).number % 2 == 0]
            if ind:
                print(compact(e), ':', compact(ee))
                print(ind)
                print()
                print(p)
                print()
                print(pp)
                print('\n')
            for i in ind:
                assert (i, i - 1) in pp.mapping and pp.entry(i, i - 1).number + 1 == pp.entry(i, i).number
                a = pp.entry(i, i).number
                d = pp.set(i, i, a - 1)
                d = pp.set(i, i - 1, a)
                d = pp.set(i - 1, i, a)
            if d.strict_upper_half() != pp.strict_upper_half():
                print(d)
                input('?')
        print()


from pipedreams import *

def print_atomic(n):
    for w in Permutation.involutions(n):
        print(w.cycle_repr())
        print('')
        dreams = [Pipedream.from_word(*e) for e in w.get_pipe_dreams()]
        for d in dreams:
            if d.is_atomic():
                print(d.atomic_part())
                print('')
        print('\n\n')


def generate_involutions(n):
    for perm in Permutation.involutions(n):
        if perm(n) != n:
            Pipedream.save_involution(perm)


def generate_fpf_involutions(n):
    for perm in Permutation.fpf_involutions(n):
        if not (perm(n) == n - 1 and perm(n - 1) == n):
            Pipedream.save_fpf_involution(perm)








from words import *

schurq = [
    {(): 1},                                        # Q_0
    {(1,): 1},                                      # Q_1
    {(1, 1): 2, (2,): 1},                           # Q_2
    {(1, 1, 1): 2, (1, 2): 1, (2, 1): 1},           # Q_21
    {(1, 1, 1): 4, (2, 1): 2, (1, 2): 2, (3,): 1},  # Q_3
    {(3,): 1},
    {(1, 1, 1, 1): 8, (2, 1, 1): 4, (1, 2, 1): 4, (1, 1, 2): 4, (2, 2): 2, (3, 1): 1, (1, 3): 1},
    {(1, 1, 1, 1): 8, (2, 1, 1): 4, (1, 2, 1): 4, (1, 1, 2): 4, (2, 2): 2, (3, 1): 2, (1, 3): 2, (4,): 1},
    {(3, 1): 1, (1, 3): 1, (4,): 1},
]


def simplify(x):
    s = {k: v for k, v in x.items()}
    for q in schurq:
        while all(s.get(k, 0) >= q[k] for k in q):
            s = {k: s[k] - q.get(k, 0) for k in s}
            s = {k: s[k] for k in s if s[k] != 0}
    return s


def span(w, n):
    # if len(w) < n:
    #     for i in range(len(w)):
    #         a = w[i]
    #         yield Word(*(w[:i] + (a,) + w[i:]))
    # for i in range(len(w) - 1):
    #     a, b = w[i], w[i + 1]
    #     if a == b:
    #         yield Word(*(w[:i] + (a,) + w[i + 2:]))
    for i in range(len(w) - 3):
        b, c, bb, a = w[i], w[i + 1], w[i + 2], w[i + 3]
        if a <= b == bb < c:
            yield Word(*(w[:i] + (a, b, c, b) + w[i + 4:]))
        a, b, c, bb = w[i], w[i + 1], w[i + 2], w[i + 3]
        if a <= b == bb < c:
            yield Word(*(w[:i] + (b, c, b, a) + w[i + 4:]))
    for i in range(len(w) - 2):
        b, bb, a = w[i], w[i + 1], w[i + 2]
        if a < b == bb:
            yield Word(*(w[:i] + (b, a, b) + w[i + 3:]))
            yield Word(*(w[:i] + (a, b, b) + w[i + 3:]))
        b, a, bb = w[i], w[i + 1], w[i + 2]
        if a < b == bb:
            yield Word(*(w[:i] + (a, b, b) + w[i + 3:]))
            yield Word(*(w[:i] + (b, b, a) + w[i + 3:]))
        a, b, bb = w[i], w[i + 1], w[i + 2]
        if a < b == bb:
            yield Word(*(w[:i] + (b, b, a) + w[i + 3:]))
            yield Word(*(w[:i] + (b, a, b) + w[i + 3:]))
        #
        c, a, b = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield Word(*(w[:i] + (a, c, b) + w[i + 3:]))
        a, c, b = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield Word(*(w[:i] + (c, a, b) + w[i + 3:]))
        b, c, a = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield Word(*(w[:i] + (b, a, c) + w[i + 3:]))
        b, a, c = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield Word(*(w[:i] + (b, c, a) + w[i + 3:]))


def get_classes(n, k):
    gen = Word.all(n, k)
    words = set(gen)
    while words:
        seed = {words.pop()}
        new = set()
        while seed:
            add = set()
            for w in seed:
                new.add(w)
                add |= set(span(w, n))
            seed = add - new
        q = Vector()
        for w in new:
            q += w.quasisymmetrize(Word.weakly_unimodal_zeta)
        yield new, q
        words = words - new


def simplify(q):
    s = {
        k: v - min({q[t] for t in itertools.permutations(k)})
        for k, v in q.items()
    }
    s = {k: v for k, v in s.items() if v}
    return s


def test(n, k):
    ans = []
    for a, q in get_classes(n, k):
        s = simplify(q)
        if s:
            ans += [(a, s)]
    for a, s in sorted(ans, key=lambda w: tuple(sorted(w[1].keys()))):
        if any(k in w for w in a):
            print(' '.join(map(str, a)), '\n  ', s, '\n')
