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

invol, skew = None, None
invol, skew = test(3, invol, skew)


SignedPermutation(-1, -2, 3, -4, -5)

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
n = 3
w = SignedPermutation.longest_element(n)
words = list(w.get_involution_words())


def tableau(word, n):
    parts = []
    for i in range(n):
        parts += [word[:i + 1]]
        word = word[i + 1:]
    parts = reversed(parts)
    string = ';'.join([','.join([str(i) for i in part]) for part in parts])
    return Tableau.from_string(string)


for e in words:
    print(tableau(e, n))
    print('')


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

extra = []
while b:
    if a[0] == b[0]:
        a = a[1:]
        b = b[1:]
    else:
        extra + =[a[0]]
        a = a[1:]






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

co = coincidences_d(5)
for w, l in co.items():
     print(w,':',l)


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
n = 3; w = EvenSignedPermutation.longest_element(n); w.get_atoms()


len(list(w.get_flattened_involution_words()))


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


ans = [(a, b) for (a, b, c) in invword(n) if not c]
for a, b in ans:
    p = Word(*a).involution_insert()[0]
    q = Word(*b).modified_hecke_insert(verbose=False)[0]
    r = Word(*b).hecke_insert()[0]
    print(p, '\n')
    print(r, '\n')
    print(q, '\n\n')
    print(a, b)


from permutations import *
from words import *

def test(n):
    g = Permutation.fpf_involutions(n)
    for w in g:
        print(w)
        for e in get_fpf_involution_words(tuple(w.oneline)):
            print(''.join(map(str, e)), ':')
            p, q = Word(*e).fpf_insert()
            pp, qq = Word(*e).alt_fpf_insert()
            print()
            print(p)
            print()
            print(pp)
            print('\n')
            if p.fpf_double() != pp:
                print(p.fpf_double())
                input('?')
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
