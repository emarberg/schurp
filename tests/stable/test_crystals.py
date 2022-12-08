from stable.crystals import (
    WordCrystalGenerator,
    TableauCrystalGenerator,
    ShiftedCrystalGenerator,
    OrthogonalSetvaluedShiftedCrystalGenerator,
    URTShiftedCrystalGenerator,
    MRTShiftedCrystalGenerator,
    MRT_Symplectic_ShiftedCrystalGenerator
)
from stable.permutations import Permutation
from stable.tableaux import Partition
from stable.insertion import InsertionAlgorithm
import pytest


def test_word_crystal_generator():
    g = WordCrystalGenerator((3, 2, 1), 3, 3)
    print(g.factorizations)
    assert set(g.factorizations) == {
        ((1, 2, 1), (2, 2, 3)),
        ((1, 2, 1), (1, 2, 3)),
        ((1, 2, 1), (1, 1, 3)),
        ((1, 2, 1), (1, 1, 2)),
        ((2, 1, 2), (2, 3, 3)),
        ((2, 1, 2), (1, 3, 3)),
        ((2, 1, 2), (1, 2, 3)),
        ((2, 1, 2), (1, 2, 2))
    }
    print(g.edges)


@pytest.mark.slow
def test_shifted_crystal_generate():
    for mvs in [False]:
        for fpf in [False]:
            ShiftedCrystalGenerator((4, 2), rank=3, excess=0, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((4, 2), rank=3, excess=1, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((4, 2), rank=3, excess=2, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((4, 2), rank=3, excess=3, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((2,), rank=3, excess=1, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((2,), rank=3, excess=2, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((2, 1), rank=3, excess=0, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((2, 1), rank=3, excess=1, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((2, 1), rank=3, excess=2, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3, 1), rank=3, excess=0, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3, 1), rank=3, excess=1, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3, 1), rank=3, excess=2, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3,), rank=3, excess=0, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3,), rank=3, excess=1, multisetvalued=mvs, is_symplectic=fpf).generate()
            ShiftedCrystalGenerator((3,), rank=3, excess=2, multisetvalued=mvs, is_symplectic=fpf).generate()


@pytest.mark.slow
def test_urt_crystal():
    for mu in Partition.all(10, max_part=4, strict=True):
        for excess in range(3):
            URTShiftedCrystalGenerator(mu, rank=3, excess=excess).generate()


@pytest.mark.slow
def test_mrt_crystal():
    for mu in Partition.all(10, max_part=4, strict=True):
        for excess in range(5):
            MRTShiftedCrystalGenerator(mu, rank=3, excess=excess).generate()


@pytest.mark.slow
def test_mrt_sp_crystal():
    for mu in Partition.all(10, max_part=4, strict=True):
        for excess in range(5):
            MRT_Symplectic_ShiftedCrystalGenerator(mu, rank=3, excess=excess).generate()


@pytest.mark.slow
def test_word_crystal_generate():
    WordCrystalGenerator((3, 2, 1), 3, 3).generate()
    WordCrystalGenerator((3, 2, 1), 3, 3, True).generate()

    WordCrystalGenerator((3, 2, 1), 4, 4).generate()
    WordCrystalGenerator((3, 2, 1), 4, 4, True).generate()

    WordCrystalGenerator((3, 2, 1), 5, 5).generate()
    WordCrystalGenerator((3, 2, 1), 5, 5, True).generate()


@pytest.mark.slow
def test_tableau_crystal_generate():
    TableauCrystalGenerator((2,), 3).generate()


@pytest.mark.slow
def test_word_crystal_locality():
    def to_tuple(word, record, k):
        tup = tuple(tuple(word[x] for x in range(len(record)) if record[x] == y) for y in range(1, k + 1))
        return ' | '.join([' '.join([str(a) for a in w]) for w in tup])

    for word in Permutation.hecke_words(5):
        for k in range(len(word)):
            for f in WordCrystalGenerator.get_increasing_factorizations(word, k):
                for op_index in range(1, len(f) - 1):
                    w, i = WordCrystalGenerator.factorization_tuple_to_array(f)
                    p, q = InsertionAlgorithm.hecke(w, i)
                    q = q.f_crystal_operator(op_index)
                    if q is not None:
                        v, j = InsertionAlgorithm.inverse_hecke(p, q)

                        omit = [op_index, op_index + 1]

                        t = {w[x] for x in range(len(w)) if i[x] not in omit and w[x] != v[x]}
                        t1 = tuple(w[x] for x in range(len(w)) if i[x] not in omit)
                        t2 = tuple(v[x] for x in range(len(v)) if j[x] not in omit)

                        if len(t) > 1:
                            print('t =', t)
                            print(op_index)
                            # print(w)
                            # print(i)
                            print(to_tuple(w, i, k))
                            print(to_tuple(v, j, k))
                            print()
                            # print(v)
                            # print(j)
                            print()
                        assert len(t) <= 1


@pytest.mark.slow
def test_decreasing_word_crystal_locality():
    def to_tuple(word, record, k):
        tup = tuple(tuple(word[x] for x in range(len(record)) if record[x] == y) for y in range(1, k + 1))
        return ' | '.join([' '.join([str(a) for a in w]) for w in tup])

    for word in Permutation.hecke_words(5):
        for k in range(len(word)):
            for f in WordCrystalGenerator.get_increasing_factorizations(word, k, decreasing=True):
                for op_index in range(1, len(f) - 1):
                    w, i = WordCrystalGenerator.factorization_tuple_to_array(f)
                    p, q = InsertionAlgorithm.hecke(w, i, multisetvalued=False)
                    q = q.f_crystal_operator(op_index)
                    if q is not None:
                        v, j = InsertionAlgorithm.inverse_hecke(p, q, multisetvalued=False)

                        omit = [op_index, op_index + 1]

                        t = {w[x] for x in range(len(w)) if i[x] not in omit and w[x] != v[x]}
                        t1 = tuple(w[x] for x in range(len(w)) if i[x] not in omit)
                        t2 = tuple(v[x] for x in range(len(v)) if j[x] not in omit)

                        if len(t) > 1:
                            print('t =', t)
                            print(op_index)
                            # print(w)
                            # print(i)
                            print(to_tuple(w, i, k))
                            print(to_tuple(v, j, k))
                            print()
                            # print(v)
                            # print(j)
                            print()
                        assert len(t) <= 1
