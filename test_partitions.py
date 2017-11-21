from partitions import StrictPartition


def test_init():
    StrictPartition(2, 1)
    StrictPartition()


def test_init_failure():
    try:
        StrictPartition(1, 1)
    except:
        pass
    else:
        raise Exception()

    try:
        StrictPartition(1, 2)
    except:
        pass
    else:
        raise Exception()


def test_repr():
    assert str(StrictPartition(2, 1)) == '* * \n  * '
    assert str(StrictPartition()) == ''


def test_shape():
    p = StrictPartition(2, 1)
    assert p.shape == {(1, 1), (1, 2), (2, 2)}


def test_row():
    p = StrictPartition(2, 1)
    assert p.row(1) == ((1, 1), (1, 2))
    assert p.row(2) == ((2, 2),)
    assert p.row(3) == ()
    assert p.row(0) == ()
    assert p.row(-1) == ()


def test_column():
    p = StrictPartition(2, 1)
    assert p.column(1) == ((1, 1),)
    assert p.column(2) == ((1, 2), (2, 2))
    assert p.column(3) == ()
    assert p.column(0) == ()
    assert p.column(-1) == ()


def test_num_rows():
    p = StrictPartition(2, 1)
    assert p.num_rows == 2

    p = StrictPartition()
    assert p.num_rows == 0


def test_num_columns():
    p = StrictPartition(2, 1)
    assert p.num_columns == 2

    p = StrictPartition()
    assert p.num_columns == 0
