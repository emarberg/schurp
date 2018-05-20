from numbers import MarkedNumber


def test_init():
    try:
        MarkedNumber('a')
    except:
        pass
    else:
        raise Exception()


def test_eq():
    a = MarkedNumber(-3)
    b = MarkedNumber(3)
    c = -3

    try:
        a == c
    except:
        pass
    else:
        raise Exception()

    assert a != b


def test_lt():
    assert MarkedNumber(-1) < MarkedNumber(1) < MarkedNumber(-2) < MarkedNumber(2)


def test_repr():
    assert str(MarkedNumber(-1)) == "1'"
    assert str(MarkedNumber(1)) == "1"


def test_sorted():
    a = MarkedNumber(-3)
    b = MarkedNumber(2)
    c = MarkedNumber(-6)
    d = MarkedNumber(11)
    e = MarkedNumber(-11)
    x = [a, b, c, d, e]
    assert list(sorted(x)) == [b, a, c, e, d]


def test_hash():
    assert hash(MarkedNumber(-3)) == -3
    assert hash(MarkedNumber(2)) == 2
