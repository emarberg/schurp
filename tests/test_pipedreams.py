from permutations import Permutation
from pipedreams import Pipedream
from schubert import InvSchubert, FPFSchubert


def test_bottom_pipedream():
    w = Permutation(3, 1, 4, 6, 5, 2)
    p = w.get_bottom_pipe_dream()
    assert p.crossings == {(1, 1), (1, 2), (3, 1), (4, 1), (4, 2), (5, 1)}


def test_ladder_moves():
    w = Permutation(1, 4, 3, 2)
    p = w.get_bottom_pipe_dream()

    q = Pipedream({(1, 2), (2, 1), (2, 2)})
    r = Pipedream({(1, 3), (2, 1), (3, 1)})
    print(p)
    print(q)
    print(r)
    print()
    for x in p.ladder_moves():
        print(x)
    assert set(p.ladder_moves()) == {q, r}

    test = set(p.upper_ladder_interval())
    assert test == set(w.get_pipe_dreams())
    expected = {
        p,
        q,
        r,
        Pipedream({(1, 2), (1, 3), (3, 1)}),
        Pipedream({(1, 2), (1, 3), (2, 2)}),
    }
    assert test == expected


def test_ladder_moves_span():
    n = 5
    for w in Permutation.all(n):
        a = set(w.get_pipe_dreams())
        b = {
            Pipedream.from_word(*dream)
            for word in w.get_reduced_words()
            for dream in w._get_pipe_dreams_helper(word)
        }
        assert a == b


def test_involutions():
    n = 6
    for i, perm in enumerate(Permutation.involutions(n)):
        print('. . .', i, perm)
        dreams = perm.get_involution_pipe_dreams()
        test_s = sum([dream.inv_monomial() for dream in dreams])
        s = InvSchubert.get(perm)
        assert s == test_s


def test_fpf_involutions():
    n = 6
    for i, perm in enumerate(Permutation.fpf_involutions(n)):
        print('. . .', i, perm)
        dreams = perm.get_fpf_involution_pipe_dreams()
        test_s = sum([dream.fpf_monomial() for dream in dreams])
        s = FPFSchubert.get(perm)
        assert s == test_s


def test_involution_ladder_moves():
    w = Permutation(1, 2, 3, 6, 5, 4)
    p = w.get_bottom_involution_pipe_dream()
    assert p == Pipedream({(5, 1), (4, 1)})

    q = Pipedream({(4, 1), (4, 2)})
    print(p)
    print()
    for d in p.involution_ladder_moves():
        print(d)
        print()
    assert q in p.involution_ladder_moves()


def test_involution_ladder_moves_span():
    n = 6
    for w in Permutation.involutions(n):
        a = set(w.get_involution_pipe_dreams())
        b = set(w.get_bottom_involution_pipe_dream().upper_involution_ladder_interval())
        if a != b:
            print(w)
            for p in a & b:
                print(p)
                print()
            print('Extra:\n')
            for p in b - a:
                print(p)
                print()
            print('Missing:\n')
            for p in a - b:
                print(p)
                print()
            print()
        assert b.issubset(a)
        assert a == b

        a = set(w.get_involution_pipe_dreams(extended=True))
        b = set(w.get_bottom_involution_pipe_dream().upper_involution_ladder_interval(extended=True))


def test_fpf_ladder_moves():
    p = Pipedream({(3, 2), (4, 2)})
    q = Pipedream({(3, 1), (3, 2)})
    print(p)
    print()
    for d in p.fpf_involution_ladder_moves():
        print(d)
        print()
    assert q in p.fpf_involution_ladder_moves()

    w = Permutation(2, 1, 6, 5, 4, 3)
    p = w.get_bottom_fpf_pipe_dream()
    assert p == Pipedream({(4, 1), (5, 1)})

    q = Pipedream({(3, 2), (5, 1)})
    print(p)
    print()
    for d in p.fpf_involution_ladder_moves():
        print(d)
        print()
    assert {q} == set(p.fpf_involution_ladder_moves())


def test_fpf_ladder_moves_span():
    n = 6
    for w in Permutation.fpf_involutions(n):
        a = set(w.get_fpf_involution_pipe_dreams())
        b = set(w.get_bottom_fpf_pipe_dream().upper_fpf_involution_ladder_interval())
        if a != b:
            print(w)
            for p in a & b:
                print(p)
                print()
            print('Extra:\n')
            for p in b - a:
                print(p)
                print()
            print('Missing:\n')
            for p in a - b:
                print(p)
                print()
            print()
        assert b.issubset(a)
        assert a == b
