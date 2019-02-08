from permutations import Permutation


def test_fpf_atoms():
    z = Permutation.cycles([[1, 8], [2, 6], [3, 10], [4, 7], [5, 9]])
    assert tuple(z.get_min_fpf_atom().inverse().oneline) == (1, 8, 2, 6, 3, 10, 4, 7, 5, 9)
    assert len(z.get_fpf_atoms()) == 8
    assert Permutation(2, 6, 4, 7, 1, 8, 5, 9, 3, 10).inverse() in z.get_fpf_atoms()
