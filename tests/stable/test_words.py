from stable.words import Word


def test_crystal_operators():
    w = (2, 1, 1)
    assert Word.e_crystal_operator(1, w) is None
    assert Word.f_crystal_operator(1, w) == (2, 1, 2)
    assert Word.e_crystal_operator(2, w) is None
    assert Word.f_crystal_operator(2, w) == (3, 1, 1)

    w = (2, 1, 2)
    assert Word.e_crystal_operator(1, w) == (2, 1, 1)
    assert Word.f_crystal_operator(1, w) is None
    assert Word.e_crystal_operator(2, w) is None
    assert Word.f_crystal_operator(2, w) == (2, 1, 3)

    w = (2, 1, 3)
    assert Word.e_crystal_operator(1, w) is None
    assert Word.f_crystal_operator(1, w) is None
    assert Word.e_crystal_operator(2, w) == (2, 1, 2)
    assert Word.f_crystal_operator(2, w) == (3, 1, 3)

    w = (3, 1, 3)
    assert Word.e_crystal_operator(1, w) is None
    assert Word.f_crystal_operator(1, w) == (3, 2, 3)
    assert Word.e_crystal_operator(2, w) == (2, 1, 3)
    assert Word.f_crystal_operator(2, w) is None

    w = (3, 1, 1)
    assert Word.e_crystal_operator(1, w) is None
    assert Word.f_crystal_operator(1, w) == (3, 1, 2)
    assert Word.e_crystal_operator(2, w) == (2, 1, 1)
    assert Word.f_crystal_operator(2, w) is None
