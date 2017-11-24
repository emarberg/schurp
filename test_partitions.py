from partitions import Partition, StrictPartition, Shape


def test_call():
    p = Partition(4, 3)
    assert p(-1) == 0
    assert p(0) == 0
    assert p(1) == 4
    assert p(2) == 3
    assert p(3) == 0
    assert p(4) == 0


def test_pieri():
    p = Partition(4, 3)

    a = Partition(5, 3)
    b = Partition(4, 4)
    c = Partition(4, 3, 1)
    assert p.pieri(1) == {a: 1, b: 1, c: 1}

    a = Partition(5, 4)
    b = Partition(5, 3, 1)
    c = Partition(4, 4, 1)
    d = Partition(4, 3, 2)
    e = Partition(6, 3)
    assert p.pieri(2) == {a: 1, b: 1, c: 1, d: 1, e: 1}


def test_strict_pieri():
    p = StrictPartition(4, 3)

    a = StrictPartition(5, 3)
    b = StrictPartition(4, 3, 1)
    assert p.pieri(1) == {a: 1, b: 1}

    a = StrictPartition(5, 3, 1)
    b = StrictPartition(6, 3)
    c = StrictPartition(4, 3, 2)
    assert p.pieri(2) == {a: 2, b: 1, c: 1}


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
    assert s.row(1) == Shape({(1, 1), (1, 2)})
    assert s.row(2) == Shape({(2, 2)})
    assert s.row(3) == Shape()
    assert s.row(0) == Shape()
    assert s.row(-1) == Shape()
    assert s.max_row == 2

def test_column():
    s = StrictPartition(2, 1).shape
    assert s.column(1) == Shape({(1, 1)})
    assert s.column(2) == Shape({(1, 2), (2, 2)})
    assert s.column(3) == Shape()
    assert s.column(0) == Shape()
    assert s.column(-1) == Shape()
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
    assert s.horizontal_border_strips() == {
        Shape({(1, 3)}),
        Shape({(1, 2), (1, 3)}),
        Shape({(1, 1), (1, 2), (1, 3)})
    }

    s = Partition(3, 3, 3, 2).shape
    assert s.horizontal_border_strips() == {
        Shape({(3, 3)}),
        Shape({(3, 3), (4, 2)}),
        Shape({(3, 3), (4, 1), (4, 2)}),
        Shape({(4, 2)}),
        Shape({(4, 1), (4, 2)})
    }


def test_vertical_border_strips():
    s = Shape()
    assert s.vertical_border_strips() == set()

    s = Partition(1, 1, 1).shape
    assert s.vertical_border_strips() == {
        Shape({(3, 1)}),
        Shape({(2, 1), (3, 1)}),
        Shape({(1, 1), (2, 1), (3, 1)})
    }

    s = Partition(4, 4, 3).shape
    assert s.vertical_border_strips() == {
        Shape({(3, 3)}),
        Shape({(3, 3), (2, 4)}),
        Shape({(3, 3), (1, 4), (2, 4)}),
        Shape({(2, 4)}),
        Shape({(1, 4), (2, 4)})
    }
