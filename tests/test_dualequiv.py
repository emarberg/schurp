from partitions import StrictPartition
from tableaux import Tableau


def assaf_dual_equivalence_op(tab, index):
    def find(t, val):
        for i, j in t:
            v = t.get(i, j)
            if abs(v) == val:
                return i, j, v, (j - i) if v.number > 0 else (i - j)

    i1, j1, v1, d1 = find(tab, index - 1)
    i2, j2, v2, d2 = find(tab, index)
    i3, j3, v3, d3 = find(tab, index + 1)

    if d1 <= d2 <= d3 or d3 <= d2 <= d1:
        return tab
    
    if d2 <= d1 <= d3 or d3 <= d1 <= d2:
        if d2 == d1:
            return tab.set(i3, j3, -v3)
        if d1 == d3:
            return tab.set(i2, j2, -v2)
        if abs(abs(d3) - abs(d2)) == 1:
            return tab.set(i2, j2, -v2).set(i3, j3, -v3)
        return tab.set(i2, j2, v2.number // abs(v2) * abs(v3)).set(i3, j3, v3.number // abs(v3) * abs(v2))
        
    if d1 <= d3 <= d2 or d2 <= d3 <= d1:
        if d2 == d3:
            return tab.set(i1, j1, -v1)
        if d1 == d3:
            return tab.set(i2, j2, -v2)
        if abs(abs(d2) - abs(d1)) == 1:
            return tab.set(i2, j2, -v2).set(i1, j1, -v1)
        return tab.set(i2, j2, v2.number // abs(v2) * abs(v1)).set(i1, j1, v1.number // abs(v1) * abs(v2))


def test_simple():
    p = Tableau.from_string("1,2',4;,3")
    q = Tableau.from_string("1,2,4;,3")
    r = Tableau.from_string("1,2,3;,4")

    assert p.dual_equivalence_operator(2) == p
    assert p.dual_equivalence_operator(1) == q

    assert q.dual_equivalence_operator(1) == p
    assert q.dual_equivalence_operator(2) == r

    assert r.dual_equivalence_operator(1) == r
    assert r.dual_equivalence_operator(2) == q

    p = Tableau.from_string("1,2,4';,3")
    q = Tableau.from_string("1,2',4';,3")
    r = Tableau.from_string("1,2',3';,4")

    assert p.dual_equivalence_operator(2) == p
    assert p.dual_equivalence_operator(1) == q

    assert q.dual_equivalence_operator(1) == p
    assert q.dual_equivalence_operator(2) == r

    assert r.dual_equivalence_operator(1) == r
    assert r.dual_equivalence_operator(2) == q

    p = Tableau.from_string("1,2',3;,4")
    q = Tableau.from_string("1,2,3';,4")

    assert p.dual_equivalence_operator(1) == q
    assert p.dual_equivalence_operator(2) == q

    assert q.dual_equivalence_operator(1) == p
    assert q.dual_equivalence_operator(2) == p


def test_simple_assaf():
    p = Tableau.from_string("1,2',4;,3")
    q = Tableau.from_string("1,2,4;,3")
    r = Tableau.from_string("1,2,3;,4")

    assert assaf_dual_equivalence_op(p, 3) == p
    assert assaf_dual_equivalence_op(p, 2) == q

    assert assaf_dual_equivalence_op(q, 2) == p
    assert assaf_dual_equivalence_op(q, 3) == r

    assert assaf_dual_equivalence_op(r, 2) == r
    assert assaf_dual_equivalence_op(r, 3) == q

    p = Tableau.from_string("1,2,4';,3")
    q = Tableau.from_string("1,2',4';,3")
    r = Tableau.from_string("1,2',3';,4")

    assert assaf_dual_equivalence_op(p, 3) == p
    assert assaf_dual_equivalence_op(p, 2) == q

    assert assaf_dual_equivalence_op(q, 2) == p
    assert assaf_dual_equivalence_op(q, 3) == r

    assert assaf_dual_equivalence_op(r, 2) == r
    assert assaf_dual_equivalence_op(r, 3) == q

    p = Tableau.from_string("1,2',3;,4")
    q = Tableau.from_string("1,2,3';,4")

    assert assaf_dual_equivalence_op(p, 2) == q
    assert assaf_dual_equivalence_op(p, 3) == q

    assert assaf_dual_equivalence_op(q, 2) == p
    assert assaf_dual_equivalence_op(q, 3) == p


def test_all(n=8):
    for m in range(n + 1):
        print('m =', m)
        for mu in StrictPartition.all(m):
            for t in Tableau.get_standard_shifted(mu):
                for i in range(1, m - 1):
                    assert t.dual_equivalence_operator(i) == assaf_dual_equivalence_op(t, i + 1)