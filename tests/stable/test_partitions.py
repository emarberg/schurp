from stable.tableaux import Partition


def test_transpose():
    mu = ()
    assert Partition.transpose(mu) == mu

    mu = (1,)
    assert Partition.transpose(mu) == mu

    mu = (2,)
    assert Partition.transpose(mu) == (1, 1)

    mu = (5, 3, 3, 2, 2, 1, 1, 1, 1)
    assert Partition.transpose(mu) == (9, 5, 3, 1, 1)


def test_complement():
    nu = ()
    assert Partition.complement(0, nu) == ()
    assert Partition.complement(1, nu) == (1,)
    assert Partition.complement(2, nu) == (2, 1)

    nu = (6, 4, 2, 1)
    assert Partition.complement(6, nu) == (5, 3)


def test_contains():
    mu = (1,)
    nu = ()
    assert Partition.contains(mu, mu)
    assert Partition.contains(mu, nu)
    assert not Partition.contains(nu, mu)

    mu = (3, 2, 1)
    nu = (3, 2)
    assert Partition.contains(mu, mu)
    assert Partition.contains(mu, nu)
    assert not Partition.contains(nu, mu)

    mu = (3, 2, 1)
    nu = (1, 1, 1)
    assert Partition.contains(mu, nu)
    assert not Partition.contains(nu, mu)

    mu = (1, 1, 1, 1)
    nu = (1, 1, 1)
    assert Partition.contains(mu, nu)
    assert not Partition.contains(nu, mu)


def test_covers():
    mu = ()
    assert set(Partition.remove_inner_corners(mu)) == {()}
    assert set(Partition.remove_shifted_inner_corners(mu)) == {()}

    mu = (1,)
    assert set(Partition.remove_inner_corners(mu)) == {(1,), ()}
    assert set(Partition.remove_shifted_inner_corners(mu)) == {(1,), ()}

    mu = (3, 2, 1)
    assert set(Partition.remove_inner_corners(mu)) == {
        (3, 2, 1),
        (2, 2, 1),
        (3, 1, 1),
        (3, 2),
        (2, 1, 1),
        (2, 2),
        (3, 1),
        (2, 1)}
    assert set(Partition.remove_shifted_inner_corners(mu)) == {
        (3, 2),
        (3, 2, 1)
    }


def test_generate_partitions():
    assert set(Partition.generate(0)) == {()}
    assert set(Partition.generate(1)) == {(1,)}
    assert set(Partition.generate(2)) == {(1, 1), (2,)}
    assert set(Partition.generate(3)) == {(1, 1, 1), (2, 1), (3,)}

    assert set(Partition.generate(4)) == {(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)}
    assert set(Partition.generate(4, max_row=0)) == set()
    assert set(Partition.generate(4, max_row=1)) == {(4,)}
    assert set(Partition.generate(4, max_row=2)) == {(2, 2), (3, 1), (4,)}
    assert set(Partition.generate(4, max_row=3)) == {(2, 1, 1), (2, 2), (3, 1), (4,)}
    assert set(Partition.generate(4, max_row=4)) == {(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)}


def test_subpartitions():
    assert set(Partition.subpartitions(())) == {()}
    assert set(Partition.subpartitions([1])) == {(), (1,)}
    assert set(Partition.subpartitions([1, 1])) == {(), (1,), (1, 1)}
    assert set(Partition.subpartitions([2, 1])) == {(), (1,), (1, 1), (2,), (2, 1)}
    assert len(list(Partition.subpartitions([2, 1]))) == 5

    assert set(Partition.subpartitions((), True)) == {()}
    assert set(Partition.subpartitions([1], True)) == {(), (1,)}
    assert set(Partition.subpartitions([1, 1], True)) == {(), (1,)}
    assert set(Partition.subpartitions([2, 1], True)) == {(), (1,), (2,), (2, 1)}
    assert len(list(Partition.subpartitions([2, 1], True))) == 4


def test_find_shifted_corners():
    mu = (3, 2, 1)
    assert Partition.find_shifted_inner_corner(mu, 0) == 3
    assert Partition.find_shifted_inner_corner(mu, 1) is None
    assert Partition.find_shifted_inner_corner(mu, 2) is None

    assert Partition.find_shifted_outer_corner(mu, 0) is None
    assert Partition.find_shifted_outer_corner(mu, 1) is None
    assert Partition.find_shifted_outer_corner(mu, 2) is None
    assert Partition.find_shifted_outer_corner(mu, 3) == 1


def test_decompose_shifted_shape_by_vertical_strips():
    f = Partition.decompose_shifted_shape_by_vertical_strips
    assert set(f(())) == {()}
    assert set(f((1,))) == {(), (1,)}
    assert set(f((2, 1))) == {(1,), (2,), (2, 1)}
    assert set(f((3, 1))) == {(2,), (2, 1), (3,), (3, 1)}


def test_decompose_shifted_shape_by_rims():
    f = Partition.decompose_shifted_shape_by_rims
    assert set(f(())) == {()}
    assert set(f((1,))) == {(), (1,)}
    assert set(f((2, 1))) == {(1,), (2,), (2, 1)}
    assert set(f((3, 1))) == {(3, 1), (3,), (2, 1), (2,), (1,)}
