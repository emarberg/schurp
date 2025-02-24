from stable.partitions import Partition
from stable.utils import *
from stable.polynomials import Polynomial
import traceback
import math


def test_onevar_GP(n=10):
    for mu in Partition.all(n, strict=True):
        print(sum(mu), 'mu =', mu)
        for nu in Partition.subpartitions(mu, strict=True):
            sh = Partition.shifted_shape(mu, nu)
            if any((i + 1, j + 1) in sh for (i, j) in sh):
                expected = Polynomial.zero()
            else:
                mush = Partition.shifted_shape(mu)
                nush = Partition.shifted_shape(nu)
                corners = {(i, j) for (i, j) in nush if (i + 1, j) not in nush and (i, j + 1) not in nush and (i + 1, j + 1) not in sh}

                begs = {(i, j) for (i, j) in sh if (i + 1, j) not in sh and (i, j - 1) not in sh and i != j}
                gaps = {(i, j) for (i, j) in corners if (i == j or (i + 1, j) in sh) and (i, j + 1) in sh}
                exts = {(i, j) for (i, j) in corners if (i + 1, j) in sh or (i, j + 1) in sh or i == j} - gaps
                iles = corners - gaps - exts

                a = len(gaps)
                b = len(begs) - len(gaps)
                c = 2 * len(corners) - len(exts) - len(gaps)

                expected = 2**a * (2 + X(0) * X(1))**b * (1 + X(0) * X(1))**c * X(1)**(sum(mu) - sum(nu))

            actual = GP_doublebar(1, mu, nu).polynomial()
            if actual != expected:
                print()
                print(actual, '==', expected)
                Partition.print_shifted(mu, nu)
                print('mu =', mu, 'nu =', nu, 'begs =', begs, 'gaps =', gaps, 'exts =', exts, 'iles =', iles)
            assert actual == expected


def test_onevar_GQ(n=10):
    for mu in Partition.all(n, strict=True):
        print(sum(mu), 'mu =', mu)
        for nu in Partition.subpartitions(mu, strict=True):
            sh = Partition.shifted_shape(mu, nu)
            if any((i + 1, j + 1) in sh for (i, j) in sh):
                expected = Polynomial.zero()
            else:
                mush = Partition.shifted_shape(mu)
                nush = Partition.shifted_shape(nu)
                corners = {(i, j) for (i, j) in nush if (i + 1, j) not in nush and (i, j + 1) not in nush and (i + 1, j + 1) not in sh}

                begs = {(i, j) for (i, j) in sh if (i + 1, j) not in sh and (i, j - 1) not in sh}
                gaps = {(i, j) for (i, j) in corners if (i + 1, j) in sh and (i, j + 1) in sh}
                exts = {(i, j) for (i, j) in corners if (i + 1, j) in sh or (i, j + 1) in sh} - gaps

                a = len(gaps)
                b = len(begs) - len(gaps)
                c = 2 * len(corners) - len(exts) - len(gaps)

                expected = 2**a * (2 + X(0) * X(1))**b * (1 + X(0) * X(1))**c * X(1)**(sum(mu) - sum(nu))

            actual = GQ_doublebar(1, mu, nu).polynomial()
            if actual != expected:
                print()
                print(actual, '==', expected)
                Partition.print_shifted(mu, nu)
                print('mu =', mu, 'nu =', nu, 'begs =', begs, 'gaps =', gaps, 'exts =', exts, 'iles =', iles)
            assert actual == expected


def test_specialization(n=10, m=5):
    for mu in Partition.all(n, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            for i in range(1, m + 1):
                f = GP(i, mu, nu).polynomial()
                f = f.set_variable(0, -1)
                for j in range(1, i + 1):
                    f = f.set_variable(j, 1)
                if f == 0:
                    print('GP', mu, '/', nu, len(mu), ':', i, ': result =', f)
                assert f in [0, 1]
                #assert f == (1 if i >= len(mu) else 0)

                f = GQ(i, mu, nu).polynomial()
                f = f.set_variable(0, -1)
                for j in range(1, i + 1):
                    f = f.set_variable(j, 1)
                print('GQ', mu, '/', nu, len(mu), ':', i, ': result =', f)
                assert f in [0, 1]
                #assert f == (1 if i >= len(mu) else 0)

            print()


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
                    


