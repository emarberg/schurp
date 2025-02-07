from stable.tableaux import Tableau, MarkedReversePlanePartition
from stable.partitions import Partition


def test_sqrt_signature_word():
    t = Tableau({
        (1, 1): 1, 
        (1, 2): (1, 2),
        (1, 3): 2,
        (1, 4): 1,
        (1, 5): (1, 2),
        (1, 6): 1,
        (1, 7): (1, 2),
        (1, 8): 2,
        (1, 9): (1, 2),
    })
    print(t)
    print(t.sqrt_f_operator(1))


def test_from_string(m=4):
    for n in range(m):
        for mu in Partition.generate(n, strict=True):
            for t in Tableau.semistandard_shifted_marked(n, mu):
                u = Tableau(str(t))
                if u != t:
                    print(t)
                    print(u)
                    print(t.boxes, u.boxes)
                    print()
                assert Tableau(str(t)) == t


def test_row_reading_word():
    t = Tableau()
    assert t.row_reading_word() == tuple()
    assert t == Tableau.from_row_reading_word(t.row_reading_word())

    t = Tableau({
        (1, 1): 1, (1, 2): 3, (1, 3): 4,
        (2, 1): 2
    })
    assert t.row_reading_word() == (2, 1, 3, 4)
    u = Tableau.from_row_reading_word(t.row_reading_word())
    assert t == u

def test_standardize():
    t = Tableau({
        (1, 1): 1, (1, 2): (-3, 2), (1, 3): (4, 5), (1, 4): (-7, -6, 7), (1, 5): (9, 10),
        (2, 2): (4, 4), (2, 3): (6, -7), (2, 4): -8,
        (3, 3): 8
    })
    assert t.standardize() == Tableau({
        (1, 1): 1, (1, 2): (-3, 2), (1, 3): (6, 7), (1, 4): (-10, -8, 12), (1, 5): (15, 16),
        (2, 2): (4, 5), (2, 3): (9, -11), (2, 4): -13,
        (3, 3): 14
    })


def test_crystal_reading_word():
    t = Tableau({
        (1, 1): 1, (1, 2): (3, 2), (1, 3): (3, 4, 5), (1, 4): (7, 6, 7), (1, 5): (9, 10),
        (2, 2): (4, 4), (2, 3): (6, 7), (2, 4): 8,
        (3, 3): 8
    })
    print(t)
    print(t.crystal_reading_word())
    assert t.is_semistandard()
    assert t.crystal_reading_word() == (
        (1,),
        (4, 2, 3, 4),
        (8, 6, 3, 4, 5, 7),
        (8, 6, 7, 7),
        (9, 10)
    )
    assert t == Tableau.from_crystal_reading_word(t.crystal_reading_word())


def test_crystal_operators():
    t = Tableau({
        (1, 1): (1, 2, 4),
        (2, 1): (5, 5, 6, 6),
        (3, 1): (7, 7, 7),
    })
    print(t)
    print(t.e_crystal_operator(4))
    assert t.e_crystal_operator(4) == Tableau({
        (1, 1): (1, 2, 4, 4),
        (2, 1): (5, 6, 6),
        (3, 1): (7, 7, 7),
    })
    assert t.f_crystal_operator(4) is None

    t = Tableau({
        (1, 1): 1, (1, 2): (1, 2, 3),
        (2, 1): 3
    })
    u = Tableau({
        (1, 1): (1, 2), (1, 2): (2, 3),
        (2, 1): 3
    })
    assert t.f_crystal_operator(1, multisetvalued=False) == u
    assert u.e_crystal_operator(1, multisetvalued=False) == t

    t = Tableau({
        (1, 1): 1, (1, 2): (1, 2),
        (2, 1): 2
    })
    u = Tableau({
        (1, 1): 1, (1, 2): (1, 3),
        (2, 1): 2
    })
    assert t.f_crystal_operator(2, multisetvalued=False) == u
    assert t.e_crystal_operator(2, multisetvalued=False) is None
    assert u.e_crystal_operator(2, multisetvalued=False) == t


def test_horizontal_strips():
    mu = tuple()
    nu = ()
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        (tuple(), {(1, 1), (1, 2)}, []),
        ((1,), {(1, 2)}, [(1, 1)]),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (1, 1)
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        ((1,), {(2, 1)}, []),
        ((1, 1), set(), [(2, 1)]),
    ]

    mu = (2, 1)
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        ((1,), {(1, 2), (2, 1)}, []),
        ((2,), {(2, 1)}, [(1, 2)]),
        ((1, 1), {(1, 2)}, [(2, 1)]),
        ((2, 1), set(), [(1, 2), (2, 1)]),
    ]

    nu = (1, 1)
    assert list(Tableau._horizontal_strips(mu, nu)) == [
        ((1, 1), {(1, 2)}, []),
        ((2, 1), set(), [(1, 2)]),
    ]


def test_shifted_horizontal_strips():
    mu = tuple()
    nu = ()
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        (tuple(), {(1, 1), (1, 2)}, []),
        ((1,), {(1, 2)}, [(1, 1)]),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (2, 1)
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        ((2,), {(2, 2)}, []),
        ((2, 1), set(), [(2, 2)]),
    ]

    mu = (3, 1)
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        ((2,), {(2, 2), (1, 3)}, []),
        ((3,), {(2, 2)}, [(1, 3)]),
        ((2, 1), {(1, 3)}, [(2, 2)]),
        ((3, 1), set(), [(1, 3), (2, 2)]),
    ]

    nu = (3,)
    assert list(Tableau._shifted_horizontal_strips(mu, nu)) == [
        ((3,), {(2, 2)}, []),
        ((3, 1), set(), [(2, 2)]),
    ]


def test_shifted_vertical_strips():
    mu = ()
    nu = ()
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        ((1,), {(1, 2)}, []),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (2, 1)
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        ((1,), {(1, 2), (2, 2)}, []),
        ((2,), {(2, 2)}, [(1, 2)]),
        ((2, 1), set(), [(2, 2)]),
    ]

    mu = (3, 1)
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        ((2,), {(2, 2), (1, 3)}, []),
        ((3,), {(2, 2)}, [(1, 3)]),
        ((2, 1), {(1, 3)}, [(2, 2)]),
        ((3, 1), set(), [(1, 3), (2, 2)]),
    ]

    nu = (1, 1)
    assert list(Tableau._shifted_vertical_strips(mu, nu)) == [
        ((2, 1), {(1, 3)}, []),
        ((3, 1), set(), [(1, 3)]),
    ]


def test_rpp_horizontal_strips():
    mu = ()
    nu = ()
    assert Tableau._rpp_horizontal_strips(mu, nu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._rpp_horizontal_strips(mu, nu) == [
        ((), {(1, 1)}),
        ((1,), set()),
    ]

    mu = (1, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu, nu)} == {
        (),
        (1, 1),
        (1,),
    }

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu, nu)} == {
        (2,),
        (1, 1),
        (1,),
        (),
        (2, 1),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu, nu)} == {
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
        (1, 1),
        (3, 1),
    }

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._rpp_vertical_strips(mu, nu)} == {
        (2,),
        (1, 1),
        (1,),
        (),
        (2, 1),
    }

    nu = (1, 1)
    assert {nu for nu, _ in Tableau._rpp_vertical_strips(mu, nu)} == {
        (1, 1),
        (2, 1),
    }


def test_shifted_rpp_horizontal_strips():
    mu = ()
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), {(1, 1)}),
        ((1,), set())
    ]

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (2, 1),
        (2,),
        (1,),
        (),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (3, 1),
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
    }


def test_shifted_rpp_vertical_strips():
    mu = ()
    assert Tableau._shifted_rpp_vertical_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._shifted_rpp_vertical_strips(mu) == [
        ((), {(1, 1)}),
        ((1,), set())
    ]

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_vertical_strips(mu)} == {
        (2, 1),
        (2,),
        (1,),
        (),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_vertical_strips(mu)} == {
        (3, 1),
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
    }


def test_semistandard():
    mu = ()
    assert Tableau.semistandard(0, mu) == {Tableau()}
    assert Tableau.semistandard(1, mu) == {Tableau()}
    assert Tableau.semistandard(2, mu) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard(1, mu) == {Tableau({(1, 1): 1})}
    assert Tableau.semistandard(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
    }

    mu = (1, 1)
    assert Tableau.semistandard(2, mu) == {
        Tableau({(1, 1): 1, (2, 1): 2})
    }

    mu = (2,)
    assert Tableau.semistandard(2, mu) == {
        Tableau({(1, 1): 1, (1, 2): 1}),
        Tableau({(1, 1): 1, (1, 2): 2}),
        Tableau({(1, 1): 2, (1, 2): 2}),
    }

    mu = (1,)
    assert Tableau.semistandard(3, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): 3}),
    }


def test_semistandard_setvalued():
    mu = ()
    assert Tableau.semistandard_setvalued(0, mu) == {Tableau()}
    assert Tableau.semistandard_setvalued(1, mu) == {Tableau()}
    assert Tableau.semistandard_setvalued(2, mu) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_setvalued(1, mu) == {Tableau({(1, 1): 1})}
    assert Tableau.semistandard_setvalued(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): (1, 2)})
    }

    mu = (1, 1)
    assert Tableau.semistandard_setvalued(2, mu) == {
        Tableau({(1, 1): 1, (2, 1): 2})
    }

    mu = (2,)
    assert Tableau.semistandard_setvalued(2, mu) == {
        Tableau({(1, 1): 1, (1, 2): 1}),
        Tableau({(1, 1): 1, (1, 2): 2}),
        Tableau({(1, 1): 2, (1, 2): 2}),
        Tableau({(1, 1): 1, (1, 2): (1, 2)}),
        Tableau({(1, 1): (1, 2), (1, 2): 2}),
    }


def test_semistandard_marked():
    mu = ()
    assert Tableau.semistandard_marked(0, mu) == {Tableau()}
    assert Tableau.semistandard_marked(1, mu) == {Tableau()}
    assert Tableau.semistandard_marked(2, mu) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_marked(1, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
    }
    assert Tableau.semistandard_marked(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
    }

    mu = (3, 1)
    print(Tableau.semistandard_marked(1, mu))
    assert Tableau.semistandard_marked(1, mu) == {
        Tableau({(1, 1): -1, (2, 1): -1, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 1): 1, (1, 2): 1, (1, 3): 1}),
    }
    mu = (3, 2)
    assert Tableau.semistandard_marked(1, mu) == set()


def test_semistandard_marked_setvalued():
    mu = ()
    assert Tableau.semistandard_marked_setvalued(0, mu) == {Tableau()}
    assert Tableau.semistandard_marked_setvalued(1, mu) == {Tableau()}
    assert Tableau.semistandard_marked_setvalued(2, mu) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_marked_setvalued(1, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)})
    }
    assert Tableau.semistandard_marked_setvalued(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
        Tableau({(1, 1): (-2, 2)}),
        Tableau({(1, 1): (1, 2)}),
        Tableau({(1, 1): (1, -2)}),
        Tableau({(1, 1): (1, -2, 2)}),
        Tableau({(1, 1): (-1, 2)}),
        Tableau({(1, 1): (-1, -2)}),
        Tableau({(1, 1): (-1, -2, 2)}),
        Tableau({(1, 1): (-1, 1, 2)}),
        Tableau({(1, 1): (-1, 1, -2)}),
        Tableau({(1, 1): (-1, 1, -2, 2)})
    }


def test_semistandard_shifted_marked():
    mu = ()
    assert Tableau.semistandard_shifted_marked(0, mu) == {Tableau()}
    assert Tableau.semistandard_shifted_marked(1, mu) == {Tableau()}
    assert Tableau.semistandard_shifted_marked(2, mu) == {Tableau()}

    mu = (1,)

    assert Tableau.semistandard_shifted_marked(1, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
    }
    assert Tableau.semistandard_shifted_marked(1, mu, diagonal_primes=False) == {
        Tableau({(1, 1): 1}),
    }

    assert Tableau.semistandard_shifted_marked(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
    }
    assert Tableau.semistandard_shifted_marked(2, mu, diagonal_primes=False) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
    }

    mu = (3, 1)
    assert Tableau.semistandard_shifted_marked(1, mu) == set()
    assert Tableau.semistandard_shifted_marked(1, mu, diagonal_primes=False) == set()

    assert Tableau.semistandard_shifted_marked(2, mu) == {
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): -2, (1, 3): 2}),
    }
    assert Tableau.semistandard_shifted_marked(2, mu, diagonal_primes=False) == {
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
    }

    mu = (3, 2)
    assert Tableau.semistandard_shifted_marked(1, mu) == set()
    assert Tableau.semistandard_shifted_marked(1, mu, diagonal_primes=False) == set()


def test_semistandard_shifted_marked_setvalued():
    mu = ()
    assert Tableau.semistandard_shifted_marked_setvalued(0, mu) == {Tableau()}
    assert Tableau.semistandard_shifted_marked_setvalued(1, mu) == {Tableau()}
    assert Tableau.semistandard_shifted_marked_setvalued(2, mu) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_shifted_marked_setvalued(1, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)})
    }
    assert Tableau.semistandard_shifted_marked_setvalued(2, mu) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
        Tableau({(1, 1): (-2, 2)}),
        Tableau({(1, 1): (1, 2)}),
        Tableau({(1, 1): (1, -2)}),
        Tableau({(1, 1): (1, -2, 2)}),
        Tableau({(1, 1): (-1, 2)}),
        Tableau({(1, 1): (-1, -2)}),
        Tableau({(1, 1): (-1, -2, 2)}),
        Tableau({(1, 1): (-1, 1, 2)}),
        Tableau({(1, 1): (-1, 1, -2)}),
        Tableau({(1, 1): (-1, 1, -2, 2)})
    }


def test_semistandard_marked_rpp():
    mu = (1,)
    print(Tableau.semistandard_marked_rpp(1, mu, diagonal_nonprimes=True))
    print()
    assert Tableau.semistandard_marked_rpp(1, mu, diagonal_nonprimes=True) == {
        MarkedReversePlanePartition({(1, 1): 1}),
        MarkedReversePlanePartition({(1, 1): -1}),
    }
    print(Tableau.semistandard_marked_rpp(1, mu, diagonal_nonprimes=False))
    print()
    assert Tableau.semistandard_marked_rpp(1, mu, diagonal_nonprimes=False) == {
        MarkedReversePlanePartition({(1, 1): -1}),
    }


def test_skew_semistandard():
    mu = ()
    nu = ()
    assert Tableau.semistandard(0, mu, nu) == {Tableau()}
    assert Tableau.semistandard(1, mu, nu) == {Tableau()}
    assert Tableau.semistandard(2, mu, nu) == {Tableau()}

    mu = (1, 1)
    nu = (1,)
    assert Tableau.semistandard(2, mu, nu) == {
        Tableau({(2, 1): 1}),
        Tableau({(2, 1): 2})
    }

    mu = (2,)
    nu = (1,)
    assert Tableau.semistandard(2, mu, nu) == {
        Tableau({(1, 2): 1}),
        Tableau({(1, 2): 2}),
    }

    mu = (2, 1)
    nu = (1,)
    assert Tableau.semistandard(2, mu, nu) == {
        Tableau({(1, 2): 1, (2, 1): 1}),
        Tableau({(1, 2): 1, (2, 1): 2}),
        Tableau({(1, 2): 2, (2, 1): 1}),
        Tableau({(1, 2): 2, (2, 1): 2}),
    }

    mu = (2, 2)
    nu = (1,)
    assert Tableau.semistandard(3, mu, nu) == {
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 2}),
        Tableau({(1, 2): 1, (2, 1): 2, (2, 2): 2}),
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 3}),
        Tableau({(1, 2): 1, (2, 1): 2, (2, 2): 3}),
        Tableau({(1, 2): 1, (2, 1): 3, (2, 2): 3}),
        Tableau({(1, 2): 2, (2, 1): 1, (2, 2): 3}),
        Tableau({(1, 2): 2, (2, 1): 2, (2, 2): 3}),
        Tableau({(1, 2): 2, (2, 1): 3, (2, 2): 3}),
    }

    mu = (2, 2, 1)
    nu = (1,)
    assert Tableau.semistandard(3, mu, nu) == {
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 2, (3, 1): 3}),
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 2, (3, 1): 2}),
        Tableau({(1, 2): 1, (2, 1): 2, (2, 2): 2, (3, 1): 3}),
        #
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 3, (3, 1): 3}),
        Tableau({(1, 2): 1, (2, 1): 1, (2, 2): 3, (3, 1): 2}),
        Tableau({(1, 2): 1, (2, 1): 2, (2, 2): 3, (3, 1): 3}),
        #
        Tableau({(1, 2): 2, (2, 1): 1, (2, 2): 3, (3, 1): 3}),
        Tableau({(1, 2): 2, (2, 1): 1, (2, 2): 3, (3, 1): 2}),
        Tableau({(1, 2): 2, (2, 1): 2, (2, 2): 3, (3, 1): 3}),
    }


def test_skew_semistandard_setvalued():
    n = 2
    mu = (2,)
    nu = (1,)
    tabs = Tableau.semistandard_setvalued(n, mu, nu)
    print(sorted(tabs))
    assert tabs == {
        Tableau({(1, 2): 1}),
        Tableau({(1, 2): 2}),
        Tableau({(1, 2): (1, 2)})
    }

    n = 2
    mu = (1, 1)
    nu = (1,)
    tabs = Tableau.semistandard_setvalued(n, mu, nu)
    print(sorted(tabs))
    assert tabs == {
        Tableau({(2, 1): 1}),
        Tableau({(2, 1): 2}),
        Tableau({(2, 1): (1, 2)})
    }


def test_skew_semistandard_marked_setvalued():
    for n in [1, 2, 3]:
        mu = (2, 1)
        nu = (1,)
        tabs = Tableau.semistandard_marked_setvalued(n, mu, nu)
        assert len(tabs) == (2**(2 * n) - 1)**2


def test_descents():
    t = Tableau({
        (1, 1): 1, (1, 2): -2, (1, 3): 3, (1, 4): -4, (1, 5): 5, (1, 6): -10,
        (2, 2): 6, (2, 3): -7, (2, 4): -9, (2, 5): -12,
        (3, 3): 8, (3, 4): -11,
    })
    assert t.descent_set() == {1, 3, 5, 6, 8, 9, 11}
