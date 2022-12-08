from stable.partitions import Partition
from stable.polynomials import beta
from stable.vectors import Vector
import itertools


SET_PARTITIONS_CACHE = {0: [()]}
BIG_MULTIPERMUTATIONS_CACHE = {}
WEAK_MULTIPERMUTATIONS_CACHE = {}

PI_P_CACHE = {}
PI_Q_CACHE = {}


def osp_printer(osp):
    return '(' + '|'.join([''.join(map(str, part)) for part in osp]) + ')'


def set_partitions(n):
    assert n >= 0
    if n not in SET_PARTITIONS_CACHE:
        ans = []
        for sp in set_partitions(n - 1):
            ans.append(sp + ((n,),))
            sp = list(sp)
            for i in range(len(sp)):
                sp[i] += (n,)
                ans.append(tuple(sp))
                sp[i] = sp[i][:-1]
        SET_PARTITIONS_CACHE[n] = ans
    for sp in SET_PARTITIONS_CACHE[n]:
        yield sp

def ordered_set_partitions(n):
    for sp in set_partitions(n):
        for osp in itertools.permutations(sp):
            yield osp


def reduce_to_big_multipermutation(osp):
    for part in osp:
        for i in range(len(part) - 1):
            if part[i + 1] - part[i] == 1:
                p = part[i + 1]
                osp = tuple(tuple(j if j < p else j - 1 for j in q if j != p) for q in osp)
                osp = tuple(q for q in osp if q)
                osp, length = reduce_to_big_multipermutation(osp)
                return osp, length + 1
    return osp, 0


def is_weak_multipermutation(osp):
    for a, part in enumerate(osp):
        for i in range(len(part) - 1):
            if part[i + 1] - part[i] == 1:
                p = part[i] - 1
                q = part[i + 1] + 1
                test = set()
                for b in range(a + 1):
                    test |= set(osp[b])
                if p in test or q in test:
                    return False
    return True


def is_big_multipermutation(osp):
    for part in osp:
        for i in range(len(part) - 1):
            if part[i + 1] - part[i] == 1:
                return False
    return True


def big_multipermutations(n):
    assert n >= 0
    if n not in BIG_MULTIPERMUTATIONS_CACHE:
        ans = list(filter(is_big_multipermutation, ordered_set_partitions(n)))
        BIG_MULTIPERMUTATIONS_CACHE[n] = ans
    for osp in BIG_MULTIPERMUTATIONS_CACHE[n]:
        yield osp


def weak_multipermutations(n):
    assert n >= 0
    if n not in WEAK_MULTIPERMUTATIONS_CACHE:
        ans = list(filter(is_weak_multipermutation, ordered_set_partitions(n)))
        WEAK_MULTIPERMUTATIONS_CACHE[n] = ans
    for osp in WEAK_MULTIPERMUTATIONS_CACHE[n]:
        yield osp


def peak_set(a):
    n = sum(map(len, a))
    mapping = {}
    for i, part in enumerate(a):
        for x in part:
            mapping[x] = i
    return {i for i in range(2, n) if mapping[i - 1] < mapping[i] > mapping[i + 1]}


def compute_pi_p(n):
    PI_P_CACHE[n] = {}
    for a in big_multipermutations(n):
        pk = ([0] + sorted(peak_set(a)) + [n]) if n > 0 else [0]
        alpha = tuple(pk[i] - pk[i - 1] for i in range(1, len(pk)))
        PI_P_CACHE[n][alpha] = PI_P_CACHE[n].get(alpha, Vector(printer=osp_printer)) + Vector({a: 1})


def compute_pi_q(n):
    PI_Q_CACHE[n] = {}
    for a in weak_multipermutations(n):
        pk = ([0] + sorted(peak_set(a)) + [n]) if n > 0 else [0]
        alpha = tuple(pk[i] - pk[i - 1] for i in range(1, len(pk)))
        a, reduction = reduce_to_big_multipermutation(a)
        coeff = 2**(len(alpha) - reduction) * beta**reduction
        PI_Q_CACHE[n][alpha] = PI_Q_CACHE[n].get(alpha, Vector(printer=osp_printer)) + Vector({a: coeff})


def pi_p(alpha):
    n = sum(alpha)
    if n not in PI_P_CACHE:
        compute_pi_p(n)
    return PI_P_CACHE[n].get(alpha, 0)


def pi_q(alpha):
    n = sum(alpha)
    if n not in PI_Q_CACHE:
        compute_pi_q(n)
    return PI_Q_CACHE[n].get(alpha, 0)


def test_pi_q_simple():
    assert pi_q((3, 1)) == 4 * pi_p((3, 1)) + 2 * beta * pi_p((2, 1))
        

def test_pi_q(nn=6):
    for n in range(nn + 1):
        print('* testing n=', n)
        for alpha in Partition.compositions(n):
            lhs = pi_q(alpha)
            if lhs == 0:    
                continue
            
            rhs = Vector()
            n = len(alpha)
            for v in range(2 ** n):
                coeff = 2 ** n
                gamma = list(alpha)
                for i in range(n):
                    if v % 2 == 0:
                        coeff = coeff // 2 * beta
                        gamma[i] -= 1
                    v = v // 2
                gamma = tuple(gamma)
                rhs += pi_p(gamma) * coeff
            assert lhs == rhs


def test_inverse_pi_q(nn=6):
    for n in range(nn + 1):
        print('* testing n=', n)
        for alpha in Partition.compositions(n):
            ell = len(alpha)
            lhs = pi_p(alpha) * 2**(n + ell)
            if lhs == 0:
                continue

            rhs = Vector()
            for m in range(sum(alpha) + 1):
                for gamma in Partition.compositions(m):
                    if len(gamma) != ell:
                        continue
                    if not all(gamma[i] <= alpha[i] for i in range(ell)):
                        continue
                    delta = sum(alpha) - sum(gamma)
                    rhs += pi_q(gamma) * 2**sum(gamma) * (-beta)**delta
            assert lhs == rhs