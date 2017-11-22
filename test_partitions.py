from partitions import Partition, StrictPartition, Shape


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


def test_corners():
    s = Shape()
    assert s.corners() == set()

    s = Partition(4, 4, 3, 2).shape
    assert s.corners() == {(2, 4), (3, 3), (4, 2)}


def test_horizontal_border_strips():
    s = Shape()
    assert s.horizontal_border_strips() == set()

    s = Partition(3).shape
    assert s.horizontal_border_strips() == {((1, 3),), ((1, 2), (1, 3)), ((1, 1), (1, 2), (1, 3))}

    s = Partition(3, 3, 3, 2).shape
    assert s.horizontal_border_strips() == {
        ((3, 3),),
        ((3, 3), (4, 2)),
        ((3, 3), (4, 1), (4, 2)),
        ((4, 2),),
        ((4, 1), (4, 2))
    }
