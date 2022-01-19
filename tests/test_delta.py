from pipedreams import BumplessPipedream
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
                i = w.inverse()(wi)
                print(w)
                print(b)
                x = b.get_pipe(wi, j, 'H')
                print(wi, j, '-->', i, '=?=', x)
                print()
                assert x == i

def test_modify_column_move_rectangle():
    x = 1
    y = 2 
    x_prime = 4
    p = 2
    # Example on P.5, Gao & Huang
    # test case for step 2: (x = 1, y = 2, x_prime = 4, p = 2)
    t={(1, 3): '┌', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', (2, 2): '┌', (4, 2): '┌', (5, 2): '│', (6, 2): '│', (7, 2): '│', (2, 3): '┼', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', (3, 5): '┌', (3, 6): '─', (3, 7): '─', (5, 4): '─', (6, 4): '┌', (7, 4): '│', (4, 6): '┌', (4, 7): '─', (5, 6): '┼', (6, 6): '│', (7, 6): '┼', (5, 3): '┌', (5, 5): '┼', (5, 7): '─', (6, 3): '│', (7, 3): '│', (6, 7): '┌', (7, 7): '┼', (7, 5): '┌', (4, 3): '┘', (3, 3): '│', (3, 2): '┘', (6, 5): '┘', (4, 5): '│'}

    t2={(1, 2): '┌', (1, 3): '─', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', (2, 2): '│', (3, 2): '┼', (4, 2): '│', (5, 2): '│', (6, 2): '│', (7, 2): '│', (2, 3): '┌', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', (3, 5): '┌', (3, 6): '─', (3, 7): '─', (5, 4): '─', (6, 4): '┌', (7, 4): '│', (4, 6): '┌', (4, 7): '─', (5, 6): '┼', (6, 6): '│', (7, 6): '┼', (5, 3): '┌', (5, 5): '┼', (5, 7): '─', (6, 3): '│', (7, 3): '│', (6, 7): '┌', (7, 7): '┼', (7, 5): '┌', (3, 3): '┘', (6, 5): '┘', (4, 5): '│'}
    
    r = x
    # test case for step 3: (x = 4, y = 4, x_prime = 6, p = 5)
    # t={(1, 2): '┌', (1, 3): '─', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', 
    #     (2, 2): '│', (3, 2): '┼', (4, 2): '│', (5, 2): '│', (6, 2): '│', (7, 2): '│', 
    #     (2, 3): '┌', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', 
    #     (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', 
    #     (3, 5): '┌', (3, 6): '─', (3, 7): '─', 
    #     (5, 4): '─', (6, 4): '┌', (7, 4): '│', 
    #     (4, 6): '┌', (4, 7): '─', 
    #     (5, 6): '┘', (6, 6): '─', (7, 6): '┌', 
    #     (5, 3): '┌', (5, 5): '┼', (5, 7): '┌', 
    #     (6, 3): '│', (7, 3): '│', 
    #     (6, 7): '┼', (7, 7): '┼', 
    #     (7, 5): '│', 
    #     (3, 3): '┘', 
    #     (6, 5): '┼', 
    #     (4, 5): '│'}

    D = BumplessPipedream(t)
    print('Input:' )
    print(D)
    E, a, r = D.delta()

    # assert E == BumplessPipedream(t2)
    # assert a == 4
    #assert r == 1

    # tiles = D.tiles.copy()
    # if p != y + 1:
    #     # pass # step 2

    #     # Find z and z_prime
    #     z = x + 1
    #     while D.get_tile(z, y + 1) != D.P_TILE and D.get_tile(z, y) == D.C_TILE:
    #         z += 1

    #     z_prime = z + 1
    #     while D.get_tile(z_prime, y) != D.J_TILE:
    #         z_prime += 1
        
    #     assert D.get_pipe(z,y) == D.get_pipe(z_prime,y)
        
    #     # Column moves
    #     del tiles[(x_prime, y + 1)]
            
    #     tiles[(x, y)] = D.C_TILE

    #     # for i in range(z+1,z_prime):
    #     #     tiles[(i, y + 1)] = self.V_TILE

    #     tiles[(z, y + 1)] = D.C_TILE
    #     tiles[(z_prime, y + 1)] = D.J_TILE

    #     for i in range(x + 1, x_prime):
    #         if (i < z and D.get_tile(i, y + 1) == D.P_TILE ) or i == z_prime:
    #             tiles[(i, y)] = D.P_TILE
    #         else:    
    #             tiles[(i, y)] = D.V_TILE

    #     if D.get_tile(x, y + 1) == D.V_TILE:
    #         tiles[(x, y + 1)] = D.J_TILE
    #     elif D.get_tile(x, y + 1) == D.C_TILE:
    #         tiles[(x, y + 1)] = D.H_TILE

    #     if D.get_tile(x_prime, y) == D.H_TILE:
    #         tiles[(x_prime, y)] = D.J_TILE
    #     elif D.get_tile(x_prime, y) == D.C_TILE:
    #         tiles[(x_prime, y)] = D.V_TILE

    #     # assign new labled blank tile
    #     (x, y) = (x_prime, y + 1)
    #     t2={(1, 2): '┌', (1, 3): '─', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', (2, 2): '│', (3, 2): '┼', (4, 2): '│', (5, 2): '│', (6, 2): '│', (7, 2): '│', (2, 3): '┌', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', (3, 5): '┌', (3, 6): '─', (3, 7): '─', (5, 4): '─', (6, 4): '┌', (7, 4): '│', (4, 6): '┌', (4, 7): '─', (5, 6): '┼', (6, 6): '│', (7, 6): '┼', (5, 3): '┌', (5, 5): '┼', (5, 7): '─', (6, 3): '│', (7, 3): '│', (6, 7): '┌', (7, 7): '┼', (7, 5): '┌', (3, 3): '┘', (6, 5): '┘', (4, 5): '│'}
    #     print('Expected pipe dream after step 2: ')
    #     print(BumplessPipedream(t2, 7))
    #     print('Output: ')
    #     print(BumplessPipedream(tiles, D.n))
    #     assert BumplessPipedream(t2, 7) == BumplessPipedream(tiles, D.n)
    #     # assert 1==0
    # else:
    #     tiles[(x, y)] = tiles[(x_prime, y + 1)] = D.C_TILE
    #     if w.get_tile(x, y + 1) == D.V_TILE:
    #         tiles[(x, y + 1)] = D.J_TILE
    #     elif w.get_tile(x, y + 1) == D.C_TILE:
    #         tiles[(x, y + 1)] = D.H_TILE

    #     for i in range(x + 1, x_prime + 1):
    #         if tiles[(i, y + 1)] == D.P_TILE:
    #             tiles[(i, y + 1)] = D.H_TILE
    #             tiles[(i, y)] = D.P_TILE
    #             tiles[(i + 1, y + 1)] = D.C_TILE
    #         else:
    #             tiles[(i, y)] = D.V_TILE

    #     global a
    #     a = y
    #     print((a, r))
    #     # print('Expected pipe dream after step 3: ')
    #     # print(BumplessPipedream(t2, 7))
    #     print('Output: ')
    #     print(BumplessPipedream(tiles, D.n))
    #     assert False



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

    
    # Example on P.5, Gao & Huang
    # test case for step 2: (x = 1, y = 2, x_prime = 4, p = 2)
    # t = {(1, 3): '┌', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─',
    #  (2, 2): '┌', (4, 2): '┌', (5, 2): '│', (6, 2): '│', (7, 2): '│',
    #   (2, 3): '┼', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', 
    #   (3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', 
    #   (3, 5): '┌', (3, 6): '─', (3, 7): '─',
    #    (5, 4): '─', (6, 4): '┌', (7, 4): '│', 
    #    (4, 6): '┌', (4, 7): '─', 
    #    (5, 6): '┘', (6, 6): '─', (7, 6): '┌', 
    #    (5, 3): '┌', (5, 5): '┼', (5, 7): '┌', 
    #    (6, 3): '│', (7, 3): '│', 
    #    (6, 7): '┼', (7, 7): '┼', 
    #    (7, 5): '│', 
    #    (4, 3): '┘', (3, 3): '│', 
    #    (3, 2): '┘', 
    #    (6, 5): '┼', (4, 5): '│'}
    # tiles of delta(D)
    # t2 = {(3, 1): '┌', (4, 1): '│', (5, 1): '│', (6, 1): '│', (7, 1): '│', 
    # (2, 2): '│', (2, 3): '┌', (2, 4): '─', (2, 5): '─', (2, 6): '─', (2, 7): '─', 
    # (3, 2): '┼', (4, 2): '│', (5, 2): '│', (6, 2): '│', (7, 2): '│', 
    # (1, 3): '─', (1, 4): '─', (1, 5): '─', (1, 6): '─', (1, 7): '─', 
    # (5, 3): '┌', (6, 3): '│', (7, 3): '│', 
    # (4, 4): '┌', (4, 6): '┌', (4, 7): '─', 
    # (5, 4): '┼', (6, 4): '│', (7, 4): '│', 
    # (3, 5): '┌', (3, 6): '─', (3, 7): '─', 
    # (6, 5): '┌', (7, 5): '│',
    #  (6, 6): '─', (6, 7): '┼', 
    #  (7, 6): '┌', (7, 7): '┼', 
    #  (5, 7): '┌', (5, 6): '┘', 
    #  (5, 5): '─', (4, 5): '┘',
    #  (3, 3): '┘', (1, 2): '┌'}
    ### a = 4, r = 1 ####
    
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

    # print('Input: ' )
    # print(D)

    
    # (x, y) = D.get_minimal_blank_tile()
    # r = x
    # while True:
    # # step 1
    #     while D.is_blank(x, y + 1):
    #         y = y + 1
    #     p = D.get_pipe(x, y + 1)
    #     print(p)

    # # step 2
    #     if p == y+1:
    #         break

    #     # Find the location of J_TILE in pipe p
    #     x_prime = x + 1
    #     while D.get_tile(x_prime, y + 1) != E.J_TILE:
    #         x_prime += 1
            
    #     D = D.modify_column_move_rectangle(x, y, x_prime, p)
    #     x, y = x_prime, y + 1

    # # step 3
    # x_prime = x + 1
    # while D.get_tile(x_prime, y + 1) != E.C_TILE or D.get_pipe(x_prime, y + 1, 'H') != y:
    #     x_prime += 1
    # D = D.modify_column_move_rectangle(x, y, x_prime, p)

    # print('Expected: ')
    # print(BumplessPipedream(t2))

    # print('Output: ')
    # print(D)


def test_get_sequence(n=5):
    pass
    # First example on P.6, Gao & Huang 
    # t = {(1, 3): '┌', (1, 4): '─', (1, 5): '─', (2, 2): '┌', (4, 2): '┌', (5, 2): '│', (2, 3): '┼', (2, 4): '─', (2, 5): '─', (3, 1): '┌', (4, 1): '│', (5, 1): '│', (3, 5): '┌', (4, 5): '┼', (5, 5): '┼', (4, 4): '┌', (5, 4): '┼', (5, 3): '┌', (4, 3): '┘', (3, 3): '│', (3, 2): '┘'}
    

    # Second example on P.6, Gao & Huang 
    # t={(1, 3): '┌', (1, 4): '─', (1, 5): '─', (2, 2): '┌', (3, 2): '│', (4, 2): '┼', (5, 2): '│', (2, 4): '┌', (2, 5): '─', (4, 1): '┌', (5, 1): '│', (3, 5): '┌', (4, 5): '┼', (5, 5): '┼', (4, 4): '┌', (5, 4): '┼', (5, 3): '┌', (3, 4): '┘', (3, 3): '┌', (2, 3): '┘', (4, 3): '┘'}
    # D = BumplessPipedream(t)
    # print('Input:' )
    # D.get_sequence()

    check = True
    count = 0
    for w1 in Permutation.all(n):
        if w1 == Permutation():
            continue
        for b in BumplessPipedream.from_permutation(w1, n):
            # for i in (1, b.n + 1):
            #     if (i,i) in b.diagram:
            #         count += 1
            # if count == b.n:
            #     break
            w = b.get_sequence()
            print(w)
            assert w == w1
    
    

