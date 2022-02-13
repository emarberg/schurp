from pipedreams import *
from permutations import Permutation

# TODO

def test_get_tile():
    pass

  
def test_is_blank():
    pass


def test_get_blank_tiles():
    pass


def test_get_minimal_blank_tile():
    pass


def test_get_pipe(n=5):
    for w in Permutation.all(n):
        for b in BumplessPipedream.from_permutation(w, n):
            for wi in range(1, n):
                j = n
                i = w(wi)
                print(w)
                print(b)
                x = b.get_pipe(wi, j, 'H')
                print(wi, j, '-->', i, '=?=', x)
                print()
                assert x == i

def test_modify_column_move_rectangle():
    # x = 1
    # y = 2 
    # x_prime = 4
    # p = 2
    # Example on P.5, Gao & Huang
    # test case for step 2: (x = 1, y = 2, x_prime = 4, p = 2)
    # t={(1, 3): '┌', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', (2, 2): '┌', (4, 2): '┌', (5, 2): '│', (6, 2): '│', (7, 2): '│', (2, 3): '┼', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', (3, 5): '┌', (3, 6): '─', (3, 7): '─', (5, 4): '─', (6, 4): '┌', (7, 4): '│', (4, 6): '┌', (4, 7): '─', (5, 6): '┼', (6, 6): '│', (7, 6): '┼', (5, 3): '┌', (5, 5): '┼', (5, 7): '─', (6, 3): '│', (7, 3): '│', (6, 7): '┌', (7, 7): '┼', (7, 5): '┌', (4, 3): '┘', (3, 3): '│', (3, 2): '┘', (6, 5): '┘', (4, 5): '│'}

    # t2={(1, 2): '┌', (1, 3): '─', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', (2, 2): '│', (3, 2): '┼', (4, 2): '│', (5, 2): '│', (6, 2): '│', (7, 2): '│', (2, 3): '┌', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', (3, 5): '┌', (3, 6): '─', (3, 7): '─', (5, 4): '─', (6, 4): '┌', (7, 4): '│', (4, 6): '┌', (4, 7): '─', (5, 6): '┼', (6, 6): '│', (7, 6): '┼', (5, 3): '┌', (5, 5): '┼', (5, 7): '─', (6, 3): '│', (7, 3): '│', (6, 7): '┌', (7, 7): '┼', (7, 5): '┌', (3, 3): '┘', (6, 5): '┘', (4, 5): '│'}
    
    t={(1, 4): '─', (1, 5): '─', (1, 6): '─', (3, 2): '─', (5, 2): '┌', (6, 2): '│', (2, 5): '┌', (2, 6): '─', (4, 1): '│', (5, 1): '│', (6, 1): '│', (3, 6): '┌', (4, 6): '┼', (5, 6): '┼', (6, 6): '┼', (4, 5): '┌', (5, 5): '┼', (6, 5): '┼', (5, 4): '┌', (6, 4): '┼', (6, 3): '┌', (3, 5): '┘', (3, 3): '┼', (3, 4): '─', (5, 3): '┘', (2, 3): '│', (4, 3): '│', (1, 3): '┌', (3, 1): '┌'}
    (x,y) = (2,4)
    x_prime = 3
    r = x

    D = SymmetricBumplessPipedream(t)
    print('Input:' )
    print(D)
    print('After step 2: ', D.modify_column_move_rectangle(x, y))
    # assert False
    # E, a, r = D.delta()





def test_delta():
    
    # Example on P.6, Gao & Huang test_modify_column_move_rectangle(x = 1, y = 2, x_prime = 4, p = 2)
    t = {(1, 3): '┌', (1, 4): '─', (1, 5): '─', (2, 2): '┌', (4, 2): '┌', (5, 2): '│', (2, 3): '┼', (2, 4): '─', (2, 5): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (3, 5): '┌', (4, 5): '┼', (5, 5): '┼', (4, 4): '┌', (5, 4): '┼', (5, 3): '┌', (4, 3): '┘', (3, 3): '│', (3, 2): '┘'}
    # tiles of delta(D)
    t2 = {(1, 2): '┌', (1, 3): '─', (1, 4): '─', (1, 5): '─', 
    (2, 2): '│', (3, 2): '┼', (4, 2): '│', (5, 2): '│', 
    (2, 3): '┌', (2, 4): '─', (2, 5): '─', 
    (3, 1): '┌', (4, 1): '│', (5, 1): '│', 
    (4, 3): '┌', (3, 5): '┌', (4, 4): '─', (5, 4): '┌', (4, 5): '┼', (5, 5): '┼', 
    (5, 3): '│', (3, 3): '┘'}
    #### a = 3, r = 1 ####
    
    D = BumplessPipedream(t)
    print('Input:' )
    print(D)
    E, a, r = D.delta()

    print("we expected output:")
    expected = BumplessPipedream(t2)
    print(expected)

    print()
    print("we got output:")
    print(E)

    assert E == expected
    assert a == 3
    assert r == 1
    # assert False


def test_get_sequence_simple(n=5):
    for z in Permutation.all(n):
        S = BumplessPipedream.from_permutation(z)
        for x in S:
            print("\n\n\nnew case:")
            x.get_sequence()
    
    
def test_get_gao_huang_pipedream(n=4):
    for z in Permutation.all(n):
        print("\n\n\nnew case:", z)
        bumpless = BumplessPipedream.from_permutation(z) # this is set
        pipedreams = z.get_pipe_dreams() # this is an iterator
        pipedreams = set(pipedreams)
        image = {x.get_gao_huang_pipedream() for x in bumpless}
        print("image:", image)
        print("expected:", pipedreams)
        assert image == pipedreams


def test_gao_huang_symmetry(n=6):
    # fails
    succeeds = True
    for z in Permutation.fpf_involutions(n):
        print("\n\n\nnew case:", z)
        bumpless = BumplessPipedream.from_permutation(z) # this is set
        image = {x.get_gao_huang_pipedream() for x in bumpless if x.is_symmetric()}
        print("preimage:", {x for x in bumpless if x.is_symmetric()})
        print("image:", image)
        print("desired:", {x for x in z.get_pipe_dreams() if x.is_symmetric()})
        b = all({x.is_symmetric() for x in image})
        print("succeeds:", b)
        succeeds = succeeds and b
    assert not succeeds 


def test_get_symmetric_pipedream(n=8):
    for z in Permutation.fpf_involutions(n):
        print("\n\n\nnew case:", z)
        bumpless = SymmetricBumplessPipedream.from_fpf_involution(z) # this is set
        symmetric = set(bumpless)
        pipedreams = z.get_fpf_involution_pipe_dreams() # this is an iterator
        pipedreams = set(pipedreams)
        image = {x.get_symmetric_pipedream() for x in symmetric}
        print("image:", len(image))
        print("expected:", len(pipedreams))
        assert len(image) == len(pipedreams)
        assert image == pipedreams


def test_s():
    w = Permutation(2,1,6,5,4,3)
    b = list(SymmetricBumplessPipedream.from_fpf_involution(w))
    for xx in b:
        D, a, r = xx.symmetric_delta()
        print('(a, r) =', (a, r))
        print(10 * "\n")

    

def test_symmetric_delta(n=6):
    # t = {(1, 3): '┌', (1, 4): '─', (1, 5): '─', (1, 6): '─', (3, 2): '─', (4, 2): '┌', (5, 2): '│', (6, 2): '│', (2, 4): '┌', (2, 5): '─', (2, 6): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (3, 6): '┌', (4, 6): '┼', (5, 6): '┼', (6, 6): '┼', (4, 5): '┌', (5, 5): '┼', (6, 5): '┼', (5, 4): '┌', (6, 4): '┼', (6, 3): '┌', (3, 4): '┘', (3, 3): '┼', (4, 3): '┘', (2, 3): '│'}
    # t = {(1, 4): '┌', (1, 5): '─', (1, 6): '─', (3, 2): '┌', (5, 2): '┌', (6, 2): '│', (2, 5): '┌', (2, 6): '─', (4, 1): '┌', (5, 1): '│', (6, 1): '│', (3, 6): '┌', (4, 6): '┼', (5, 6): '┼', (6, 6): '┼', (4, 5): '┌', (5, 5): '┼', (6, 5): '┼', (5, 4): '┌', (6, 4): '┼', (6, 3): '┌', (3, 5): '┘', (3, 3): '┼', (3, 4): '─', (5, 3): '┘', (2, 3): '┌', (4, 3): '│', (4, 2): '┘', (2, 4): '┘'}
    # X = SymmetricBumplessPipedream(t)
    for z in Permutation.fpf_involutions(n):
        S = SymmetricBumplessPipedream.from_fpf_involution(z)
        for X in S:
            X.symmetric_get_sequence()
            # input("\n\n\n(press enter)\n")
    # assert False

    
