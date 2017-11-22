from partitions import StrictPartition, Shape


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
    assert str(StrictPartition(2, 1)) == '* *\n  *'
    assert str(StrictPartition()) == ''


def test_shape():
    p = StrictPartition(2, 1)
    assert p.shape.positions == {(1, 1), (1, 2), (2, 2)}


def test_row():
    s = StrictPartition(2, 1).shape
    assert s.row(1) == {(1, 1), (1, 2)}
    assert s.ordered_row(1) == sorted({(1, 1), (1, 2)})
    assert s.row(2) == {(2, 2)}
    assert s.row(3) == set()
    assert s.row(0) == set()
    assert s.row(-1) == set()
    assert s.max_row == 2


def test_column():
    s = StrictPartition(2, 1).shape
    assert s.column(1) == {(1, 1)}
    assert s.column(2) == {(1, 2), (2, 2)}
    assert s.ordered_column(2) == sorted({(1, 2), (2, 2)})
    assert s.column(3) == set()
    assert s.column(0) == set()
    assert s.column(-1) == set()
    assert s.max_column == 2


def test_empty_shape():
    s = Shape()
    assert s.max_row == 0
    assert s.max_column == 0
    assert s.positions == set()
