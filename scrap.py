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


