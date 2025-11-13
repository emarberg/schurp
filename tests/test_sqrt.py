from crystals import(
    AbstractGLCrystal,
)
from permutations import Permutation
from partitions import Partition


def test_sqrt_demazure_tableaux(n=3, k=5):
    partitions = sorted({mu.transpose().tuple() for i in range(k + 1) for mu in Partition.all(i, max_part=n)})
    print(partitions)
    for mu in partitions:
        print('n =', n, 'mu =', mu)
        c = AbstractGLCrystal.sqrtcrystal_from_partition(mu, n)
        hw = c.get_highest_weight_elements()
        for u in hw:
            for w in Permutation.all(n):
                d = {}
                for i in w.get_reduced_words():
                    d[i] = c.demazure(u, *i)
                for i in d:
                    for j in d:
                        if i < j:
                            assert AbstractGLCrystal.find_isomorphism(d[i], d[j], u, u) is not None


def test_sqrt_demazure_words(n=3, k=5):
    print('n =', n, 'k =', k)
    c = AbstractGLCrystal.sqrtcrystal_of_words(k, n)
    hw = c.get_highest_weight_elements()
    for u in hw:
        for w in Permutation.all(n):
            d = {}
            for i in w.get_reduced_words():
                d[i] = c.demazure(u, *i)
            for i in d:
                for j in d:
                    if i < j:
                        assert AbstractGLCrystal.find_isomorphism(d[i], d[j], u, u) is not None