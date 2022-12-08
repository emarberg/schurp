from stable.dualops import sh_horizontal_strip, sh_vertical_strip
from stable.partitions import Partition


partitions = list(Partition.all(30, strict=True))

def test_simple():
    op = sh_horizontal_strip(0, 3)
    assert op(()) == (4,)
    assert op((1,)) is None
    assert op((6,)) == (6, 4)

    op = sh_horizontal_strip(1, 2)
    assert op((4, 1)) == (4, 3)

    op = sh_vertical_strip(0, 0)
    assert op(()) == (1,)
    assert op((1,)) is None
    assert op((2,)) == (2, 1)

    op = sh_vertical_strip(0, 2)
    assert op((2, 1)) == (3, 2, 1)

    op = sh_vertical_strip(0, 1)
    assert op((9, 8, 5, 4)) is None


def getnonzero_hv(a=4, b=2, g=partitions, condition=None):
    for aa in range(a + 1):
        for bb in range(b + 1):
            for cc in range(a + 1):
                for dd in range(b + 1):
                    boolean = True if condition is None else condition(aa, aa + bb, cc, cc + dd)
                    if not boolean:
                        continue
                    h = sh_horizontal_strip(aa, aa + bb)
                    v = sh_vertical_strip(cc, cc + dd)
                    val = {mu for mu in g if h(v(mu)) is not None}
                    if val:
                        print('h[%s %s] and v[%s %s]' % (aa, aa + bb, cc, cc + dd))
                        yield ((aa, aa + bb, cc, cc + dd), val)


def getnonzero_vh(a=4, b=2, g=partitions, condition=None):
    for aa in range(a + 1):
        for bb in range(b + 1):
            for cc in range(a + 1):
                for dd in range(b + 1):
                    boolean = True if condition is None else condition(aa, aa + bb, cc, cc + dd)
                    if not boolean:
                        continue
                    h = sh_horizontal_strip(aa, aa + bb)
                    v = sh_vertical_strip(cc, cc + dd)
                    val = {mu for mu in g if v(h(mu)) is not None}
                    if val:
                        print('h[%s %s] and v[%s %s]' % (aa, aa + bb, cc, cc + dd))
                        yield ((aa, aa + bb, cc, cc + dd), val)


def separated(a, b, c, d):
    return max(a, b) + 1 < min(c, d) or max(c, d) + 1 < min(a, b)


def adjacent(a, b, c, d):
    return max(a, b) + 1 == min(c, d) or max(c, d) + 1 == min(a, b)


def intersect(a, b, c, d):
    one = set(range(a, b + 1))
    two = set(range(c, d + 1))
    return len(one & two) > 0


def intersect_with_zero_condition(a, b, c, d):
    return intersect(a, b, c, d) and len([x for x in [a, b, c, d] if x != 0]) != 1


def printex(mapping, both=0):
    for a, b, c, d in mapping:
        h = sh_horizontal_strip(a, b)
        v = sh_vertical_strip(c, d)
        for mu in mapping[a, b, c, d]:
            if both >= 0:
                print('h[%s %s] v[%s %s]' % (a, b, c, d))
                print(Partition.printable(mu, shifted=True))
                print()
                print(Partition.printable(v(mu), shifted=True))
                print()
                print(Partition.printable(h(v(mu)), shifted=True))
                print()
                print()
                print()
            if both <= 0:
                print('v[%s %s] h[%s %s]' % (c, d, a, b))
                print(Partition.printable(mu, shifted=True))
                print()
                print(Partition.printable(h(mu), shifted=True))
                print()
                print(Partition.printable(v(h(mu)), shifted=True))
                print()
                print()
                print()
            input()


def getnoncommuting(a=4, b=2, g=partitions, condition=None):
    for aa in range(a + 1):
        for bb in range(b + 1):
            for cc in range(a + 1):
                for dd in range(b + 1):
                    boolean = True if condition is None else condition(aa, aa + bb, cc, cc + dd)
                    if not boolean:
                        continue
                    h = sh_horizontal_strip(aa, aa + bb)
                    v = sh_vertical_strip(cc, cc + dd)
                    val = {mu for mu in g if h(v(mu)) != v(h(mu))}
                    if val:
                        print('h[%s %s] and v[%s %s]' % (aa, aa + bb, cc, cc + dd))
                        yield ((aa, aa + bb, cc, cc + dd), val)


def getnoncommuting_h(a=4, b=2, g=partitions, condition=None):
    for aa in range(a + 1):
        for bb in range(b + 1):
            for cc in range(a + 1):
                for dd in range(b + 1):
                    boolean = True if condition is None else condition(aa, aa + bb, cc, cc + dd)
                    if not boolean:
                        continue
                    h1 = sh_horizontal_strip(aa, aa + bb)
                    h2 = sh_horizontal_strip(cc, cc + dd)
                    val = {mu for mu in g if h1(h2(mu)) != h2(h1(mu))}
                    if val:
                        print('h[%s %s] and h[%s %s]' % (aa, aa + bb, cc, cc + dd))
                        yield ((aa, aa + bb, cc, cc + dd), val)


def getnoncommuting_v(a=4, b=2, g=partitions, condition=None):
    for aa in range(a + 1):
        for bb in range(b + 1):
            for cc in range(a + 1):
                for dd in range(b + 1):
                    boolean = True if condition is None else condition(aa, aa + bb, cc, cc + dd)
                    if not boolean:
                        continue
                    v1 = sh_vertical_strip(aa, aa + bb)
                    v2 = sh_vertical_strip(cc, cc + dd)
                    val = {mu for mu in g if v1(v2(mu)) != v2(v1(mu))}
                    if val:
                        print('v[%s %s] and v[%s %s]' % (aa, aa + bb, cc, cc + dd))
                        yield ((aa, aa + bb, cc, cc + dd), val)


def op_equal(a1, b1, c1, d1, a2, b2, c2, d2, g=partitions):
    h1 = sh_horizontal_strip(a1, b1)
    v1 = sh_vertical_strip(c1, d1)
    h2 = sh_horizontal_strip(a2, b2)
    v2 = sh_vertical_strip(c2, d2)

    ans0 = not any(h1(v1(mu)) != v2(h2(mu)) for mu in g)
    ans1 = any(h1(v1(mu)) is not None for mu in g)
    ans2 = any(v2(h2(mu)) is not None for mu in g)
    return ans0, ans1, ans2


def find_op_equal(a=4, b=4, g=partitions, condition=adjacent):
    patterns = [
        (aa, aa + bb, cc, cc + dd)
        for aa in range(a + 1)
        for bb in range(b + 1)
        for cc in range(a + 1)
        for dd in range(b + 1)
        if condition is None or condition(aa, aa + bb, cc, cc + dd)
    ]
    for a1, b1, c1, d1 in patterns:
        for a2, b2, c2, d2 in patterns:
            ans0, ans1, ans2 = op_equal(a1, b1, c1, d1, a2, b2, c2, d2, g)
            if ans0:
                print("h[%s %s] v[%s %s] == v[%s %s] h[%s %s]" % (a1, b1, c1, d1, a2, b2, c2, d2), ans1, ans2)
                assert b1 + d1 - c1 - a1 == b2 + d2 - c2 - a2
                # yield (a1, b1, c1, d1, a2, b2, c2, d2)


def test_intersect_zero_hv(a=5, b=5):
    assert not any(getnonzero_hv(a, b, condition=intersect))


def test_intersect_zero_vh(a=5, b=5):
    assert not any(getnonzero_vh(a, b, condition=intersect_with_zero_condition))


def test_separated_commuting(a=5, b=5):
    assert not any(getnoncommuting(a, b, condition=separated))


def test_separated_commuting_h(a=5, b=5):
    assert not any(getnoncommuting_h(a, b, condition=separated))


def test_separated_commuting_v(a=5, b=5):
    assert not any(getnoncommuting_v(a, b, condition=separated))


def test_adjacent(a=5, b=5):
    assert not any(getnoncommuting_v(a, b, condition=separated))
