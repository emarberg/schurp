from crystals import (
    AbstractGLCrystal, 
    AbstractQCrystal,
    AbstractPrimedQCrystal
)
from binfty import InfiniteCrystal
from permutations import Permutation
from stable.partitions import Partition
from stable.tableaux import Tableau


def test_demazure_braids(n=4, h=8):
    words = Permutation.longest_element(n).get_reduced_words()
    y = InfiniteCrystal.sqrt_tableaux(n)
    for thresh in range(h + 1):
        boolean = True
        for a in words:
            for b in words:
                if a < b:    
                    boolean &= y.is_demazure_isomorphic(thresh, a, b)
        print('rank =', n, 'threshhold =', thresh, boolean)


def embed_q(n, mu, tab, extended=False):
    if tab is None:
        return None
    try:
        rows = tab.get_rows(unpack=False)
    except TypeError:
        rows = tab.get_rows()
    while len(rows) < len(mu):
        rows.append([])
    for i, m in enumerate(mu):
        rows[i] = m * (n - i,) + tuple(rows[i])
        if extended and len(rows[i]) == m:
            rows[i] = (m - 1) * (n - i,) + (i - n,)
    return tab.from_rows(rows, shifted=True)


def embed_qplus(n, mu, tab):
    return embed_q(n, mu, tab, True)


def test_embed_q(n=3, k=10, cls=AbstractQCrystal):
    embed = embed_q if cls is AbstractQCrystal else embed_qplus
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    partitions = list(filter(Partition.is_strict, partitions))
    for index, lam in enumerate(partitions):
        b = cls.decomposition_tableaux_from_strict_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    try:
                        if u is not None:
                            assert embed(n, mu, u) == b.e_operator_on_decomposition_tableaux(i, embed(n, mu, t))
                    except:
                        print('mu =', mu, 'e_{%s}' % i)
                        print(t)
                        print(u)
                        print(embed(n, mu, t))
                        print(embed(n, mu, u))
                        print('yet:')
                        print(b.e_operator_on_decomposition_tableaux(i, embed(n, mu, t)))
                        return False

                    try:
                        u = b.f_operator(i, t)
                        if u is not None:
                            assert embed(n, mu, u) == b.f_operator_on_decomposition_tableaux(i, embed(n, mu, t))
                    except:
                        print('mu =', mu, 'f_{%s}' % i)
                        print(t)
                        print(u)
                        print(embed(n, mu, t))
                        print(embed(n, mu, u))
                        print('yet:')
                        print(b.f_operator_on_decomposition_tableaux(i, embed(n, mu, t)))
                        return False
    return True


def test_embed_qplus(n=2, k=10):
    # fails
    assert not test_embed_q(n, k, AbstractPrimedQCrystal)


def test_embed_sqrt_q(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    partitions = list(filter(Partition.is_strict, partitions))
    for index, lam in enumerate(partitions):
        b = AbstractQCrystal.sqrtcrystal_from_strict_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_q(n, mu, u) == embed_q(n, mu, t).sqrt_decomposition_e_operator(i)
                    
                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_q(n, mu, u) == embed_q(n, mu, t).sqrt_decomposition_f_operator(i)


def embed_gl(mu, tab):
    if tab is None:
        return None
    rows = tab.get_rows(unpack=False)
    while len(rows) < len(mu):
        rows.append([])
    for i, m in enumerate(mu):
        rows[i] = m * (i + 1,) + tuple(rows[i])
    return Tableau.from_rows(rows)


def test_embed_gl(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    for index, lam in enumerate(partitions):
        b = AbstractGLCrystal.from_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == b.e_operator_on_semistandard_tableaux(i, embed_gl(mu, t))

                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == b.f_operator_on_semistandard_tableaux(i, embed_gl(mu, t))


def test_embed_sqrt_gl(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    for index, lam in enumerate(partitions):
        b = AbstractGLCrystal.sqrtcrystal_from_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == embed_gl(mu, t).sqrt_e_operator(i)
                    
                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == embed_gl(mu, t).sqrt_f_operator(i)


