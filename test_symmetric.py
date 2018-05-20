from symmetric import (
    SchurQ, Vector
)


def test_schurq_multiplication():
    assert SchurQ(1) * SchurQ() == Vector({SchurQ(1): 1})
    assert SchurQ(1) * SchurQ(1) == 2 * SchurQ(2)
    assert SchurQ(1) * SchurQ(2) == SchurQ(2, 1) + 2 * SchurQ(3)
    assert SchurQ(1) * SchurQ(3) == SchurQ(3, 1) + 2 * SchurQ(4)
    assert SchurQ(1) * SchurQ(4) == SchurQ(4, 1) + 2 * SchurQ(5)
    assert SchurQ(1) * SchurQ(5) == SchurQ(5, 1) + 2 * SchurQ(6)

    assert SchurQ(2) * SchurQ(2) == 2 * SchurQ(3, 1) + 2 * SchurQ(4)
    assert SchurQ(2) * SchurQ(3) == SchurQ(3, 2) + 2 * SchurQ(4, 1) + 2 * SchurQ(5)
    assert SchurQ(2) * SchurQ(4) == SchurQ(4, 2) + 2 * SchurQ(5, 1) + 2 * SchurQ(6)
    assert SchurQ(2) * SchurQ(5) == SchurQ(5, 2) + 2 * SchurQ(6, 1) + 2 * SchurQ(7)

    assert SchurQ(3) * SchurQ(3) == 2 * SchurQ(4, 2) + 2 * SchurQ(5, 1) + 2 * SchurQ(6)
    assert SchurQ(3) * SchurQ(4) == SchurQ(4, 3) + 2 * SchurQ(5, 2) + 2 * SchurQ(6, 1) + 2 * SchurQ(7)
    assert SchurQ(3) * SchurQ(5) == SchurQ(5, 3) + 2 * SchurQ(6, 2) + 2 * SchurQ(7, 1) + 2 * SchurQ(8)

    assert SchurQ(4) * SchurQ(4) == 2 * SchurQ(5, 3) + 2 * SchurQ(6, 2) + 2 * SchurQ(7, 1) + 2 * SchurQ(8)
    assert SchurQ(4) * SchurQ(5) == SchurQ(5, 4) + 2 * SchurQ(6, 3) + 2 * SchurQ(7, 2) + 2 * SchurQ(8, 1) + 2 * SchurQ(9)

    assert SchurQ(5) * SchurQ(5) == 2 * SchurQ(6, 4) + 2 * SchurQ(7, 3) + 2 * SchurQ(8, 2) + 2 * SchurQ(9, 1) + 2 * SchurQ(10)
