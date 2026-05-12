from permutations import Permutation
from even import EvenSignedPermutation
import subprocess


def aspan(w, symmetric=True):
    for i in range(len(w) - 2):
        if symmetric:
            c, a, b = w[i], w[i + 1], w[i + 2]
            if a < b < c:
                yield w[:i] + (b, c, a) + w[i + 3:]

        b, c, a = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield w[:i] + (c, a, b) + w[i + 3:]


def dspan(w, symmetric=True):
    for i in range(len(w) - 2):
        if symmetric:
            c, a, b = w[i], w[i + 1], w[i + 2]
            if a < b < c:
                yield w[:i] + (b, c, a) + w[i + 3:]

        b, c, a = w[i], w[i + 1], w[i + 2]
        if a < b < c:
            yield w[:i] + (c, a, b) + w[i + 3:]

    if len(w) >= 3:
        if symmetric:
            c, a, b = -w[0], w[1], w[2]
            if a < b < c:
                yield (-b, c, a) + w[3:]

        b, c, a = -w[0], w[1], w[2]
        if a < b < c:
            yield (-c, a, b) + w[3:]


def get_aclasses(n):
    ans = []
    tuples = {tuple([w(i) for i in range(1, n + 1)]) for w in Permutation.all(n)}
    while tuples:
        ans.append(set())
        add = {tuples.pop()}
        while add:
            nextadd = set()
            for w in add:
                ans[-1].add(w)
                for v in aspan(w):
                    nextadd.add(v)
            add = nextadd - ans[-1]
        tuples -= ans[-1]
    return ans


def get_dclasses(n):
    ans = []
    tuples = {tuple(w) for w in EvenSignedPermutation.all(n)}
    while tuples:
        ans.append(set())
        add = {tuples.pop()}
        while add:
            nextadd = set()
            for w in add:
                ans[-1].add(w)
                for v in dspan(w):
                    nextadd.add(v)
            add = nextadd - ans[-1]
        tuples -= ans[-1]
    return ans


def is_agood(a):
    for w in a:
        for i in range(len(w) - 2):
            c, b, a = w[i], w[i + 1], w[i + 2]
            if c > b > a:
                return False
        if len(w) >= 3:
            c, b, a = -w[0], w[1], w[2]
            if c > b > a:
                return False
    return True


def is_dgood(a):
    for w in a:
        for i in range(len(w) - 2):
            c, b, a = w[i], w[i + 1], w[i + 2]
            if c > b > a:
                return False
        if len(w) >= 3:
            c, b, a = -w[0], w[1], w[2]
            if c > b > a:
                return False
        #for i in range(len(w) - 3):
        #    a, b, d, c = w[i], w[i + 1], -w[i + 2], -w[i + 3]
        #    if 0 < a < b < c < d:
        #        return False
        #if len(w) >= 4:
        #    a, b, d, c = -w[0], w[1], -w[2], -w[3]
        #    if 0 < a < b < c < d:
        #        return False
    return True


def test_aatoms(n=4):
    for a in get_aclasses(n):
        if is_agood(a):
            for w in a:
                sigma = Permutation(w).inverse()
                assert sigma.is_atom()
        else:
            for w in a:
                sigma = Permutation(w).inverse()
                assert not sigma.is_atom()



def test_datoms(n=4):
    for a in get_dclasses(n):
        if is_dgood(a):
            for w in a:
                sigma = EvenSignedPermutation(*w).inverse()
                if not sigma.is_atom():
                    ddraw(a)
                    input('?')
                    break
                #assert sigma.is_atom()
        else:
            for w in a:
                sigma = EvenSignedPermutation(*w).inverse()
                if sigma.is_atom():
                    ddraw(a)
                    input('??')
                    break

def good_aclasses(n):
    a = get_aclasses(n)
    return sorted([x for x in a if is_agood(x)], key=len, reverse=True)


def good_dclasses(n):
    a = get_dclasses(n)
    return sorted([x for x in a if is_dgood(x)], key=len, reverse=True)


def ddraw(a):
    def printer(w):
        ans = str(w)
        #for b in EvenSignedPermutation(*w).get_reduced_words():
        #    ans += '\n' + str(b)
        return ans

    e = {(w, v, False) for w in a for v in dspan(w, False)}
    s = []
    s += ['digraph G {']
    s += ['    overlap=false;']
    s += ['    splines=spline;']
    s += ['    node [shape=box; fontname="courier"; style=filled];']
    for x in a:
        s += ['    "%s" [fillcolor=white];' % printer(x)]
    for x, y, b in e:
        s += ['    "%s" -> "%s" [style="%s"];' % (printer(x), printer(y), 'dotted' if b else 'bold')]
    s += ['}']
    s = '\n'.join(s)

    name = 'temp'
    file = '/Users/emarberg/examples/atoms/'
    dotfile = file + 'dot/DI/' + name + '.dot'
    pngfile = file + 'png/DI/' + name + '.png'
    with open(dotfile, 'w') as f:
        f.write(s)
    subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])
    subprocess.run(["open", pngfile])
