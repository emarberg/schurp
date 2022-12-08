from stable.partitions import Partition
from stable.utils import *
import traceback


def partition_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows):
            yield mu, ()
        n += 1


def strict_partition_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows, strict=True):
            yield mu, ()
        n += 1


def skew_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows):
            for nu in Partition.subpartitions(mu):
                yield mu, nu
        n += 1


def strict_skew_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows, strict=True):
            for nu in Partition.subpartitions(mu, strict=True):
                yield mu, nu
        n += 1


# s, P, Q, g, j, G, J, gp, gq, jp, gq, GP, GQ, gp*gp, gq*gq, GP*GP, GQ*GQ, GS?, skew
data = {
    's': (partition_iterator, s, schur_expansion),
    'P': (strict_partition_iterator, P, P_expansion),
    'Q': (strict_partition_iterator, Q, Q_expansion),
    'j': (partition_iterator, j, j_expansion),
    'g': (partition_iterator, g, g_expansion),
    'G': (partition_iterator, G, G_expansion),
    'jp': (strict_partition_iterator, jp, jp_expansion),
    'jq': (strict_partition_iterator, jq, jq_expansion),
    'gp': (strict_partition_iterator, gp, gp_expansion),
    'gq': (strict_partition_iterator, gq, gq_expansion),
    'GP': (strict_partition_iterator, GP, GP_expansion),
    'GQ': (strict_partition_iterator, GQ, GQ_expansion),
#    'mp_g': (partition_iterator, mp_g, mp_g_expansion),
#    'mp_gp': (strict_partition_iterator, mp_gp, mp_gp_expansion),
#    'mp_gq': (strict_partition_iterator, mp_gq, mp_gq_expansion),
    'mn_G': (partition_iterator, mn_G, mn_G_expansion),
    'mn_GP': (strict_partition_iterator, mn_GP, mn_GP_expansion),
    'mn_GQ': (strict_partition_iterator, mn_GQ, mn_GQ_expansion),
    'skew_G': (skew_iterator, G, None),
    'skew_GP': (strict_skew_iterator, GP, None),
    'skew_GQ': (strict_skew_iterator, GQ, None),
    'ss_skew_G': (skew_iterator, G_doublebar, None),
    'ss_skew_GP': (strict_skew_iterator, GP_doublebar, None),
    'ss_skew_GQ': (strict_skew_iterator, GQ_doublebar, None),
    'skew_g': (skew_iterator, g, None),
    'skew_gp': (strict_skew_iterator, gp, None),
    'skew_gq': (strict_skew_iterator, gq, None),
#    'S': (strict_partition_iterator, S, S_expansion),
#    'gs': (strict_partition_iterator, gs, gs_expansion),
#    'js': (strict_partition_iterator, js, js_expansion),
#    'GS': (strict_partition_iterator, GS, GS_expansion),
#    'skew_GS': (strict_skew_iterator, GS, None),
#    'ss_skew_GS': (strict_skew_iterator, GS_doublebar, None),
#    'skew_gs': (strict_skew_iterator, gs, None),
#    'GS GS': (strict_skew_iterator, lambda n, mu, nu: GS(n, mu) * GS(n, nu), None)
}


def decompose(n, iterator, function, decomp):
    mu, nu = next(iterator)
    f = function(n, mu, nu)
    try:
        expansion = decomp(f)
        return expansion.is_nonnegative()
    except Exception:
        return False


def update(n, results, trials_left):
    print()
    pairs = []
    for (x, y) in results:
        if results[x, y] and not any(x != z != y and results[x, z] and results[z, y] for z in data):
            print('(', x, ')', '-->', '(', y, ')')
            pairs.append((x, y))
    print()
    # for (x, y) in results:
    #     if not results[x, y]:
    #         print('(', x, ')', '-/->', '(', y, ')')
    print()
    print('n =', n, 'trials left:', trials_left)
    print()
    print(repr(pairs))
    print()
    print(len(pairs))
    print()


def test_positivity(nn, trials=1000):
    expected = None #[('s', 'g'), ('s', 'mn_G'), ('P', 's'), ('Q', 'P'), ('j', 's'), ('G', 's'), ('jp', 's'), ('jp', 'gp'), ('jq', 's'), ('jq', 'gp'), ('jq', 'gq'), ('GP', 'G'), ('GQ', 'G'), ('skew_G', 'G'), ('skew_GP', 'GP'), ('skew_GQ', 'GQ'), ('ss_skew_G', 'G'), ('ss_skew_GP', 'GP'), ('ss_skew_GQ', 'GQ'), ('skew_g', 'g'), ('skew_gp', 'gp'), ('skew_gq', 'gq')]
    results = {(x, y): True for x in data for y in data if x != y}
    for n in range(1, nn + 1):
        iterators = {name: val[0](n) for name, val in data.items()}
        for i in range(trials):
            for x in data:
                it = iterators[x]
                _, fn, _ = data[x]
                for y in data:
                    if expected is not None and (x, y) not in expected:
                        results[x, y] = False
                    if x == y or not results[x, y]:
                        continue
                    _, _, dec = data[y]
                    results[x, y] &= decompose(n, it, fn, dec)
            update(n, results, trials - i)
                    


