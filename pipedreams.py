from schubert import InvSchubert, FPFSchubert, Grothendieck, AltInvGrothendieck
from polynomials import x as x_var, one as one_var, y as y_var
import subprocess
import os
from permutations import *


class BumplessPipedream:

    J_TILE = '┘'
    C_TILE = '┌'
    P_TILE = '┼'
    H_TILE = '─'
    V_TILE = '│'
    B_TILE = '•'

    E_TILE = '*'

    TILES = [J_TILE, C_TILE, P_TILE, H_TILE, V_TILE, B_TILE, E_TILE]

    def __init__(self, bends, n=None):
        self.bends = {p: t for (p, t) in bends.items() if t in [self.J_TILE, self.C_TILE, self.E_TILE]}
        self.n = n if n else max([0] + [max(p) for p in bends])
        
        self.tiles = self.bends.copy()
        for column in range(1, self.n + 1):
            south = True
            i = self.n
            j = column
            while i > 0 and j <= self.n:
                if south:
                    if (i, j) not in self.tiles:
                        self.tiles[i, j] = self.V_TILE
                        i -= 1
                    elif self.tiles[i, j] == self.H_TILE:
                        self.tiles[i, j] = self.P_TILE
                        i -= 1
                    elif self.tiles[i, j] in [self.C_TILE, self.E_TILE]:
                        j += 1
                        south = False
                    else:
                        raise Exception
                else:
                    if (i, j) not in self.tiles:
                        self.tiles[i, j] = self.H_TILE
                        j += 1
                    elif self.tiles[i, j] == self.V_TILE:
                        self.tiles[i, j] = self.P_TILE
                        j += 1
                    elif self.tiles[i, j] in [self.J_TILE, self.E_TILE]:
                        i -= 1
                        south = True
                    else:
                        raise Exception

    def __repr__(self):
        ans = []
        for i in range(1, self.n + 1):
            row = []
            for j in range(1, self.n + 1):
                t = self.get_tile(i, j)
                row += [t]
            ans += [''.join(row)]
        return '\n' + '\n'.join(ans) + '\n'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        assert type(other) == type(self)
        return str(self) == str(other)

    def kweight(self):
        ans = x_var(0)**0
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if self.get_tile(i, j) == self.B_TILE:
                    ans *= Grothendieck.beta * (x_var(i) + y_var(j) + Grothendieck.beta * x_var(i) * y_var(j))
                if self.get_tile(i, j) == self.J_TILE:
                    ans *= (1 + Grothendieck.beta * x_var(i)) * (1 + Grothendieck.beta * y_var(j))
        return ans

    def weight(self):
        ans = x_var(0)**0
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if self.get_tile(i, j) == self.B_TILE:
                    ans *= (x_var(i) - y_var(j))
        return ans

    def fpf_weight(self):
        ans = x_var(0)**0
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if i > j and self.get_tile(i, j) == self.B_TILE:
                    ans *= (x_var(i) + x_var(j))
        return ans

    def inv_weight(self):
        ans = x_var(0) ** 0
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if self.get_tile(i, j) == self.B_TILE:
                    if i > j:
                        ans *= (x_var(i) + x_var(j))
                    elif i == j:
                        ans *= x_var(i)
        return ans

    def inv_kweight(self):
        beta = AltInvGrothendieck.beta
        ans = x_var(0)**0
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if self.get_tile(i, j) == self.B_TILE:
                    if i > j:
                        ans *= beta * (x_var(i) + x_var(j) + beta * x_var(i) * x_var(j))
                    elif i == j:
                        ans *= beta * x_var(i)
                if self.get_tile(i, j) == self.J_TILE:
                    if i > j:
                        ans *= (1 + beta * x_var(i)) * (1 + beta * x_var(j))
                    elif i == j:
                        ans *= 1 + beta * x_var(i)
        return ans

    @classmethod
    def rothe(cls, w, n=None):
        n = n if n else w.rank
        bends = {}
        for i in range(1, n + 1):
            j = w(i)
            bends[(i, j)] = cls.C_TILE
        return cls(bends, n)

    def kdroop(self, i, j, a, b, strict=True):
        ans = self.droop(i, j, a, b, strict=strict)
        if ans is not None:
            return ans

        p = [p for p in range(j + 1, b) if self.get_tile(a, p) == self.C_TILE]
        q = [q for q in range(i + 1, a) if self.get_tile(q, b) == self.C_TILE]

        if self.get_tile(i, j) != self.C_TILE:
            return None
        if self.get_tile(a, b) != self.J_TILE:
            return None
        if len(p) == 0 and len(q) == 0:
            return None
            
        c_count = 0
        j_count = 0
        for x in range(i, a + 1):
            for y in range(j, b + 1):
                if not strict and (x, y) in [(i, b), (a, j)] and self.get_tile(x, y) == self.J_TILE:
                    continue
                t = self.get_tile(x, y)
                if t == self.C_TILE:
                    c_count += 1
                if t == self.J_TILE:
                    j_count += 1
        if c_count != 2 or j_count != 1:
            return None

        bends = self.bends.copy()
        del bends[(i, j)]

        if p:
            p = p[0]
            bends[i, p] = self.C_TILE
            del bends[(a, p)]
            if self.get_tile(a, j) == self.J_TILE:
                del bends[(a, j)]
            else:
                bends[(a, j)] = self.C_TILE
        else:
            q = q[0]
            bends[q, j] = self.C_TILE
            del bends[(q, b)]
            if self.get_tile(i, b) == self.J_TILE:
                del bends[(i, b)]
            else:
                bends[(i, b)] = self.C_TILE

        return BumplessPipedream(bends, self.n)

    def droop(self, i, j, a, b, strict=True):
        if self.get_tile(i, j) != self.C_TILE:
            return None
        if self.get_tile(a, b) != self.B_TILE:
                return None
        for x in range(i, a + 1):
            for y in range(j, b + 1):
                if (x, y) == (i, j):
                    continue
                if not strict and (x, y) in [(i, b), (a, j)] and self.get_tile(x, y) == self.J_TILE:
                    continue
                if self.get_tile(x, y) in [self.C_TILE, self.J_TILE]:
                    return None

        bends = self.bends.copy()

        del bends[(i, j)]
        bends[(a, b)] = self.J_TILE
        
        if self.get_tile(i, b) == self.J_TILE:
            del bends[(i, b)]
        else:
            bends[(i, b)] = self.C_TILE

        if self.get_tile(a, j) == self.J_TILE:
            del bends[(a, j)]
        else:
            bends[(a, j)] = self.C_TILE

        return BumplessPipedream(bends, self.n)

    def droops(self, strict=True, ktheoretic=False):
        for (i, j) in self.tiles:
            for a in range(i + 1, self.n + 1):
                for b in range(j + 1, self.n + 1):
                    bpd = (self.kdroop if ktheoretic else self.droop)(i, j, a, b, strict=strict)
                    if bpd is not None:
                        yield bpd

    def symmetric_droops(self, strict=True, ktheoretic=False):
        for (i, j) in self.tiles:
            if i < j:
                continue
            for a in range(i + 1, self.n + 1):
                for b in range(j + 1, self.n + 1):
                    bpd = (self.symmetric_kdroop if ktheoretic else self.symmetric_droop)(i, j, a, b, strict=strict)
                    if bpd is not None:
                        yield bpd

    def symmetric_droop(self, i, j, a, b, strict=True):
        if i == j and a == b:
            ans = self.droop(i, j, a, b, strict=strict)
        elif i == j and a != b:
            ans = self.droop(i, j, a, b, strict=strict)
            if ans is not None:
                ans = ans.droop(i, b, b, a, strict=strict) if a > b else ans.droop(a, j, b, a, strict=strict)
        else:
            ans = self.droop(i, j, a, b, strict=False)
            if ans is not None:
                ans = ans.kdroop(j, i, b, a, strict=strict)
        return ans

    def symmetric_kdroop(self, i, j, a, b, strict=True):
        ans = self.symmetric_droop(i, j, a, b, strict=strict)
        if ans is not None:
            return ans

        if i == j and a == b:
            ans = self.kdroop(i, j, a, b, strict=strict)
        elif i == j and a != b:
            ans = self.kdroop(i, j, a, b, strict=strict)
            if ans is not None:
                ans = ans.kdroop(i, b, b, a, strict=strict) if a > b else ans.kdroop(a, j, b, a, strict=strict)
        else:
            ans = self.kdroop(i, j, a, b, strict=False)
            if ans is not None:
                ans = ans.kdroop(j, i, b, a, strict=strict)
        return ans

    @classmethod
    def from_involution(cls, w, n=None, reduced=True, strict=True):
        assert w == w.inverse()
        ans = set()
        seed = {cls.rothe(w, n)}
        while seed:
            new_seed = set()
            for bpd in seed:
                ans.add(bpd)
                new_seed |= set(bpd.symmetric_droops(strict=strict, ktheoretic=not reduced))
            seed = new_seed - ans
        return ans

    @classmethod
    def from_permutation(cls, w, n=None, reduced=True, strict=True):
        ans = set()
        seed = {cls.rothe(w, n)}
        while seed:
            new_seed = set()
            for bpd in seed:
                ans.add(bpd)
                new_seed |= set(bpd.droops(strict=strict, ktheoretic=not reduced))
            seed = new_seed - ans
        return ans

    @classmethod
    def transpose_tile(cls, t):
        assert t in cls.TILES
        if t == cls.V_TILE:
            return cls.H_TILE
        if t == cls.H_TILE:
            return cls.V_TILE
        return t

    def transpose(self):
        tiles = {}
        for (i, j) in self.tiles:
            t = self.tiles[(i, j)]
            u = self.transpose_tile(t)
            tiles[(j, i)] = u
        return self.__class__(tiles, self.n)

    def is_symmetric(self):
        return self == self.transpose()

    @classmethod
    def from_fpf_involution_slow(cls, z, n=None):
        ans = {SymmetricBumplessPipedream(bpd.tiles, bpd.n) for bpd in cls.from_permutation(z, n) if bpd.is_symmetric()}
        return ans

    # def follow(self, i, j):
    #     if i == self.n:
    #         return j
    #     if (i, j) not in self.tiles:
    #         return 0
    #     t = self.tiles[(i, j)]
    #     if t in [self.C_TILE, self.V_TILE]:
    #         return self.follow(i + 1, j)
    #     if t in [self.J_TILE, self.H_TILE]:
    #         return self.follow(i, j - 1)
    #     return 0

    def modify_column_move_rectangle(self, x, y):
        # do column move on self and determine the second step from P.4

        # Find the location of J_TILE in pipe p
        x_prime = x + 1
        while self.get_tile(x_prime, y + 1) != self.J_TILE:
            x_prime += 1

        tiles = self.tiles.copy()

        # Column moves
        del tiles[(x_prime, y + 1)]
        tiles[(x, y)] = self.C_TILE

        # Modifying first two tiles
        if self.get_tile(x, y + 1) == self.V_TILE:
            tiles[(x, y + 1)] = self.J_TILE
        elif self.get_tile(x, y + 1) == self.C_TILE:
            tiles[(x, y + 1)] = self.H_TILE
        ######### Bugs in here
        # Modifying last two tiles
        if self.get_tile(x_prime, y) == self.H_TILE:
            tiles[(x_prime, y)] = self.J_TILE
        elif self.get_tile(x_prime, y) == self.C_TILE:
            tiles[(x_prime, y)] = self.V_TILE

        for i in range(x + 1, x_prime):
            if self.get_tile(i,y) == self.H_TILE and self.get_tile(i, y + 1) == self.P_TILE:
                tiles[(i, y)] = self.P_TILE
                tiles[(i, y + 1)] = self.H_TILE
            elif self.get_tile(i,y) == self.B_TILE and self.get_tile(i, y + 1) == self.V_TILE:
                tiles[(i, y)] = self.V_TILE
                tiles[(i, y + 1)] = self.B_TILE

        # step 2
        # 2a: find out the J_TILE in column y+1 and the C_TIle in column y

        # Find z 
        z_values = [
            z for z in range(x + 1 , x_prime) 
            if self.get_tile(z, y + 1) == self.P_TILE and self.get_tile(z , y) == self.C_TILE
        ]

        for z in z_values:
            # Find z_prime
            z_prime = z + 1
            while self.get_tile(z_prime, y) != self.J_TILE:
                z_prime += 1

            assert self.get_pipe(z,y) == self.get_pipe(z_prime,y)
            tiles[(z, y + 1)] = self.C_TILE
            tiles[(z, y)] = self.V_TILE
            tiles[(z_prime, y + 1)] = self.J_TILE
            tiles[(z_prime, y)] = self.P_TILE            

        # assign new labled blank tile
        (x, y) = (x_prime, y + 1)

        # Test
        return BumplessPipedream(tiles, self.n), x_prime

    def modify_column_move_rectangle_step_three(self, x, y):
        # do column move on self and determine the third step from P.4

        x_prime = x + 1
        while self.get_tile(x_prime, y + 1) != self.P_TILE or self.get_pipe(x_prime, y + 1, 'H') != y:
            x_prime += 1

        tiles = self.tiles.copy()
        
        # step 3
        tiles[(x, y)] = self.C_TILE
        tiles[(x_prime, y)] = self.V_TILE
        tiles[(x_prime, y + 1)] = self.C_TILE   # ┌
        # Modifying the tile next to the labled tile
        if self.get_tile(x, y + 1) == self.V_TILE:
            tiles[(x, y + 1)] = self.J_TILE    # ┘
        elif self.get_tile(x, y + 1) == self.C_TILE:    # ┌
            tiles[(x, y + 1)] = self.H_TILE
        
        # Modifying the tile between row x+1 and x_prime
        for i in range(x + 1, x_prime):
            # ┌┼ becomes │┌
            if self.get_tile(i, y) == self.C_TILE and self.get_tile(i, y + 1) == self.P_TILE:
                tiles[(i, y)] = self.V_TILE
                tiles[(i, y + 1)] = self.C_TILE
            # ┘│ becomes ┼┘   
            elif self.get_tile(i, y) == self.J_TILE and self.get_tile(i, y + 1) == self.V_TILE:
                tiles[(i, y)] = self.P_TILE
                tiles[(i, y + 1)] = self.J_TILE
            else:
                tiles[(i, y)] = self.get_tile(i, y + 1)
                tiles[(i, y + 1)] = self.get_tile(i, y)

        return BumplessPipedream(tiles, self.n)
                    
    def delta(self):
        D = self
        (x, y) = D.get_minimal_blank_tile()
        r = x
        while True:
            # step 1
            while D.is_blank(x, y + 1):
                y = y + 1
            p = D.get_pipe(x, y + 1)

            # step 2
            if p == y+1:
                break

            D, x_prime = D.modify_column_move_rectangle(x, y)
            x, y = x_prime, y + 1
            # print("After step 2: ", D)

        # step 3
        a = y

        D = D.modify_column_move_rectangle_step_three(x, y)
        # print("After step 3: ", D)
        
        # step 4
        return D, a, r

    def symmetric_delta(self, verbose=True):
        D = self
        (x, y) = D.get_minimal_nondiagonal_blank_tile()
        assert x % 2 != 0
        r = 1 + (x - 1) // 2
        # if D.diagram == {(i,i):(i,i) for i in range(1, D.n) if i % 2 == 1}:
        #    return D, a, r
        if verbose:
            print(D)
        while True:
            # step 1
            while D.is_blank(x, y + 1):
                if x == y + 1:
                    x = x - 1
                y = y + 1
            p = D.get_pipe(x, y + 1)

            # step 2
            if p == y+1:
                break

            D, x_prime = D.modify_column_move_rectangle(x, y)
            if verbose:
                print("After step 2: ", D)
            D = D.symmetric_modify_column_move_rectangle(x, y)

            x, y = x_prime, y + 1
            
            if verbose:
                print("After symmetric step 2: ", D)
                print("Symmetric?", D.is_symmetric())
                print('x =', x, 'y =', y, 'x\' =', x_prime)

        # step 3
        a = y

        D = D.modify_column_move_rectangle_step_three(x, y)
        if verbose:
            print("After step 3: ", D)
        D = D.symmetric_modify_column_move_rectangle_step_three(x, y)
        if verbose:
            print("After symmetric step 3: ", D)
            print("Symmetric?", D.is_symmetric())
        # assert D.is_symmetric()
        
        # step 4
        return D, a, r

    def get_sequence(self):
        D = self
        ans = []

        while D.has_blank_tiles():
            print(D)
            D,a,r = D.delta()
            print('↓' + str((a,r)))
            ans.append(a)
        
        # for m in range(1, self.n + 1):
        #     oneline.append(self.get_pipe(m, self.n))
        print(D)

        w = Permutation.from_word(*ans)
        print('w = ',w)
        # print("Corresponding pipe dream: ",self.get_gao_huang_pipedream())

    def symmetric_get_sequence(self):
        D = self
        ans = []

        while D.has_nondiagonal_blank_tiles(): # and D.diagram != {(i,i):(i,i) for i in range(D.n) if i % 2 == 1}:
            print(D)
            D,a,r = D.symmetric_delta()
            print('↓' + str((a,r)))
            ans.append(a)
        
        # for m in range(1, self.n + 1):
        #     oneline.append(self.get_pipe(m, self.n))
        print(D)

        w = Permutation.from_word(*ans)
        print('w = ',w)

    def get_gao_huang_pipedream(self):
        # return the pipe dream after applying delta
        D = self
        crossings = []
        while D.has_blank_tiles():
            D, a, r = D.delta()
            crossings.append((r, a - (r - 1)))
        return Pipedream(crossings)

    def get_symmetric_pipedream(self, verbose=False):
        assert self.is_symmetric()
        D = self
        crossings = []
        while D.has_nondiagonal_blank_tiles():
            if verbose:
                print(D)
            D, a, r = D.symmetric_delta(verbose)
            assert a != r
            if verbose:
                print('↓' + str((a,r)))
            crossings.append((a - (r - 1), r))
        if verbose:
            print(D)
            print()
            print()
            print(crossings)
        return Pipedream(crossings)

    def symmetric_modify_column_move_rectangle(self, x, y):
        return self.transpose().modify_column_move_rectangle(x, y)[0].transpose()

    def symmetric_modify_column_move_rectangle_step_three(self, x, y):
        return self.transpose().modify_column_move_rectangle_step_three(x, y).transpose()

    def get_permutation(self):
        oneline = [self.get_pipe(i, self.n, 'H') for i in range(1, self.n + 1)]
        return Permutation(*oneline)

    def get_pipe(self, i, j, direction=None):
        # returns the column index of the position on the bottom side where the pipe enters the n-by-n grid
        assert direction in [None, 'V', 'H']
        if i == self.n + 1:
            return j
        
        t = self.get_tile(i, j)
        assert t != self.B_TILE

        if t == self.E_TILE:
            assert direction is not None
            if direction == 'V':
                return self.get_pipe(i, j-1, 'H')
            else:
                return self.get_pipe(i+1, j, 'V')

        if t == self.P_TILE:
            assert direction is not None
            if direction == 'V':
                return self.get_pipe(i+1, j, direction)
            else:
                return self.get_pipe(i, j-1, direction)

        if t == self.J_TILE: # ┘
            return self.get_pipe(i, j-1, 'H')
        elif t == self.C_TILE: # ┌
            return self.get_pipe(i+1, j, 'V')
        elif t == self.H_TILE: # ─
            return self.get_pipe(i, j-1, 'H')
        elif t == self.V_TILE: # │
            return self.get_pipe(i+1, j, 'V')

    def has_blank_tiles(self):
        return len(self.get_blank_tiles()) > 0

    def has_nondiagonal_blank_tiles(self):
        return not all(i == j for (i, j) in self.get_blank_tiles())

    def get_tile(self, i, j):
        assert 1 <= i <= self.n and 1 <= j <= self.n
        return self.tiles.get((i, j), self.B_TILE)

    def is_blank(self, i, j):
        return self.get_tile(i,j) == self.B_TILE

    def get_blank_tiles(self):
        # returns a list with the position of the rightmost tile
        # the list is order such that [row1: rightmost tiles ... leftmost tiles, row2:rightmost tiles ... leftmost tiles,... row n:rightmost tiles ... leftmost tiles] 
        return [(i,j) for i in range(1, self.n + 1) for j in reversed(range(1, self.n + 1)) if self.is_blank(i,j)]

        # the dictionary has order such that {row i: [rightmost tiles (i, y) ... leftmost tiles (i, y0)]}
        # return {i: [(i,j) for i in range(1, self.n + 1) for j in reversed(range(1, self.n + 1)) if self.is_blank(i, j)]}

    def get_minimal_nondiagonal_blank_tile(self):
        return [(i, j) for (i, j) in self.get_blank_tiles() if i != j][0]

    def get_minimal_blank_tile(self):
        # intentially will cause an error if there are no blank tiles
        # find the first row containing contiguous blank tiles
        return self.get_blank_tiles()[0]
        
        # for i in range(1, self.n + 1):
        #     if len(self.get_blank_tiles()[i]) > 1:
        #         return self.get_blank_tiles()[i][0]

    def copy(self):
        return self.__class__(self.tiles, self.n)



class SymmetricBumplessPipedream(BumplessPipedream):

    @classmethod
    def from_fpf_involution(cls, z, n=None):
        ans = set()
        seed = {cls.rothe(z, n)}
        while seed:
            new_seed = set()
            for bpd in seed:
                ans.add(bpd)
                new_seed |= set(bpd.symmetric_droops())
            seed = new_seed - ans
        return ans

    def symmetric_droops(self):
        for (i, j) in self.tiles:
            if i <= j:
                continue
            for a in range(i + 1, self.n + 1):
                for b in range(j + 1, self.n + 1):
                    bpd = self.droop(i, j, a, b)
                    if bpd is not None:
                        bpd = bpd.droop(j, i, b, a)
                        if bpd is not None:
                            yield SymmetricBumplessPipedream(bpd.tiles, bpd.n)

    def __repr__(self):
        ans = []
        for i in range(1, self.n + 1):
            row = []
            for j in range(1, self.n + 1):
                t = self.get_tile(i, j)
                if i > j:
                    t = '\033[91m' + t + '\033[0m'
                # if i == j:
                #    t = '┐' if t == self.P_TILE else ' '
                # elif i < j:
                #    t = ' '
                # if i > j and t == self.B_TILE:
                #    t = str(len({self.follow(i, a) for a in range(1, j)} - {0}))
                row += [t]
            ans += [''.join(row)]
        return '\n' + '\n'.join(ans) + '\n'


class Pipedream:

    DIRECTORY = '/Users/emarberg/examples/pipedreams/examples/'

    @classmethod
    def from_word(cls, *words):
        crossings = []
        for i, word in enumerate(words):
            assert all(word[j] > word[j + 1] for j in range(len(word) - 1))
            assert all(w - i >= 0 for w in word)
            for w in word:
                crossings.append((i + 1, w - i))
        return Pipedream(crossings)

    def __init__(self, crossings):
        self.crossings = {(i, j) for (i, j) in crossings}
        self.n = max([i + j for (i, j) in self.crossings]) if self.crossings else 1

    def count_diagonal(self):
        return len([i for (i, j) in self.crossings if i == j])

    def __iter__(self):
        return self.crossings.__iter__()

    def transpose(self):
        return Pipedream({(j, i) for (i, j) in self.crossings})

    def reflect(self):
        return Pipedream(self.crossings | {(j, i) for (i, j) in self.crossings})

    def __repr__(self):
        s = []
        for i in range(1, self.n + 1):
            s += [(self.n + 1 - i) * ['.']]
        for (i, j) in self.crossings:
            s[i - 1][j - 1] = '+'
        return '\n' + '\n'.join(' '.join(row) for row in s) + '\n'

    def __hash__(self):
        return hash(tuple(sorted(self.crossings)))

    def __eq__(self, other):
        assert type(other) == Pipedream
        return self.crossings == other.crossings

    def fpf_involution_ladder_moves(self, extended=False):
        cr = self.crossings
        for (i, j) in cr:
            if (i, j + 1) in cr:
                continue
            x = i - 1
            while (x, j) in cr and (x, j + 1) in cr:
                x = x - 1
            if x == 0:
                continue
            if j > 1 and (x, j - 1) not in cr and (x, j) in cr and (x, j + 1) not in cr and (x, j + 2) not in cr:
                #     j j+1
                # x . + .
                #   . + +
                #   . + +
                #   . + +
                # i . + .
                p = any((x - d, j - 2 + d) in self.crossings for d in range(1, x))
                q = any((x - d, j - 1 + d) in self.crossings for d in range(1, x))
                r = any((x - d, j + d) in self.crossings for d in range(1, x))
                s = any((x - d, j + 1 + d) in self.crossings for d in range(1, x))
                t = any((x - d, j + 2 + d) in self.crossings for d in range(1, x))
                if not (p or q or r or s or t):
                    yield Pipedream((self.crossings - {(i, j)}) | {(x, j - 1)})

            # if (x - 1, j) not in cr and (x - 1, j + 1) in cr and (x, j) not in cr and (x, j + 1) not in cr:
            #     x = x - 1
            #     #   j j+1
            #     #   . . .
            #     # x . + .
            #     #   . . .
            #     #   + +
            #     #   + +
            #     #   + +
            #     # i +
            #     p = any((x - d, j - 1 + d) in self.crossings for d in range(1, x))
            #     q = any((x - d, j + d) in self.crossings for d in range(1, x))
            #     r = any((x - d, j + 1 + d) in self.crossings for d in range(1, x))
            #     s = any((x + 1 - d, j + 1 + d) in self.crossings for d in range(1, x + 1))
            #     t = any((x + 2 - d, j + 1 + d) in self.crossings for d in range(1, x + 2))
            #     if not (p or q or r or s or t):
            #         yield Pipedream((self.crossings - {(i, j)}) | {(x, j)})
            #     x = x + 1

            if (x, j) not in cr and (extended or x > j + 1):
                yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

    def upper_fpf_involution_ladder_interval(self, extended=False):
        level = {self}
        while level:
            new_level = set()
            for dream in level:
                yield dream
                for e in dream.fpf_involution_ladder_moves(extended):
                    new_level.add(e)
            level = new_level

    def involution_ladder_moves(self, extended=False):
        for (i, j) in self.crossings:
            if (i, j + 1) in self.crossings:
                continue
            x = i - 1
            while (x, j) in self.crossings and (x, j + 1) in self.crossings:
                x = x - 1
            if x == 0:
                continue
            if (x, j) in self.crossings and (x, j + 1) not in self.crossings and (x, j + 2) not in self.crossings:
                #   j j+1
                #   . .
                # x + . .
                #   + +
                #   + +
                #   + +
                # i +
                p = any((x - d, j - 1 + d) in self.crossings for d in range(1, x))
                q = any((x - d, j + d) in self.crossings for d in range(1, x))
                r = any((x - d, j + 1 + d) in self.crossings for d in range(1, x))
                s = any((x - d, j + 2 + d) in self.crossings for d in range(1, x))
                if not (p or q or r or s) and (extended or x >= j + 1):
                    yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

            if (x, j) not in self.crossings and (extended or x >= j + 1):
                yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

    # def involution_downladder_moves(self, extended=False):
    #     for (i, j) in self.crossings:
    #         j = j - 1
    #         if j == 0:
    #             continue
    #         x = i + 1
    #         while (x, j) in self.crossings and (x, j + 1) in self.crossings:
    #             x = x + 1
    #         if (x, j) not in self.crossings and (x, j + 1) not in self.crossings and (i, j) in self.crossings and (i, j + 2) not in self.crossings:
    #             #   j j+1
    #             #   . .
    #             # i + +
    #             #   + +
    #             #   + +
    #             #   + +
    #             # x . .
    #             p = any((i - d, j - 1 + d) in self.crossings for d in range(1, i))
    #             q = any((i - d, j + d) in self.crossings for d in range(1, i))
    #             r = any((i - d, j + 1 + d) in self.crossings for d in range(1, i))
    #             s = any((i - d, j + 2 + d) in self.crossings for d in range(1, i))
    #             if not (p or q or r or s):
    #                 yield Pipedream((self.crossings - {(i, j + 1)}) | {(x, j)})

    #         if (i, j) not in self.crossings and (x, j + 1) not in self.crossings:
    #             yield Pipedream((self.crossings - {(i, j + 1)}) | {(x, j)})

    # def lower_involution_ladder_interval(self, extended=False):
    #     level = {self}
    #     while level:
    #         new_level = set()
    #         for dream in level:
    #             yield dream
    #             for e in dream.involution_downladder_moves(extended):
    #                 new_level.add(e)
    #         level = new_level

    def upper_involution_ladder_interval(self, extended=False):
        level = {self}
        while level:
            new_level = set()
            for dream in level:
                yield dream
                for e in dream.involution_ladder_moves(extended):
                    new_level.add(e)
            level = new_level

    def involution_chute_moves(self):
        for (i, j) in self.crossings:
            x = j - 1
            while (i, x) in self.crossings and (i + 1, x) in self.crossings:
                x = x - 1

            if not (x == 0 or (i, x) in self.crossings):
                if (i + 1, x) not in self.crossings and (i, x) not in self.crossings and (i - 1, j) not in self.crossings:
                    #   x       j
                    # i . + + + +
                    #   . + + + +
                    #
                    p = any((i + 2 - d, j + d) in self.crossings for d in range(1, i + 2))
                    q = any((i + 1 - d, j + d) in self.crossings for d in range(1, i + 1))
                    r = any((i - d, j + d) in self.crossings for d in range(1, i))
                    s = any((i - 1 - d, j + d) in self.crossings for d in range(1, i - 1))
                    if not (p or q or r or s):
                        yield Pipedream((self.crossings - {(i, j)}) | {(i + 1, x)})

                if (i + 1, j) not in self.crossings and (i, x) not in self.crossings and (i + 1, x) not in self.crossings:
                    yield Pipedream((self.crossings - {(i, j)}) | {(i + 1, x)})

            if i == 1 or (i - 1, j) in self.crossings:
                continue

            x = j + 1
            while (i - 1, x) in self.crossings and (i, x) in self.crossings:
                x = x + 1

            if (i, x) in self.crossings and (i - 1, x) not in self.crossings and (i - 2, x) not in self.crossings:
                #   j       x
                #   . + + + .
                # i + + + + +
                #
                p = any((i + 1 - d, x + d) in self.crossings for d in range(1, i + 1))
                q = any((i - d, x + d) in self.crossings for d in range(1, i))
                r = any((i - 1 - d, x + d) in self.crossings for d in range(1, i - 1))
                s = any((i - 2 - d, x + d) in self.crossings for d in range(1, i - 2))
                if not (p or q or r or s):
                    yield Pipedream((self.crossings - {(i, j)}) | {(i - 1, x)})

            if (i, x) not in self.crossings and (i - 1, x) not in self.crossings:
                yield Pipedream((self.crossings - {(i, j)}) | {(i - 1, x)})

    def involution_chute_span(self):
        level = {self}
        seen = set()
        while level:
            new_level = set()
            for dream in level:
                yield dream
                seen.add(dream)
                for e in dream.involution_chute_moves():
                    if e not in seen:
                        new_level.add(e)
            level = new_level

    def ladder_moves(self):
        for (i, j) in self.crossings:
            if (i, j + 1) in self.crossings:
                continue
            x = i - 1
            while (x, j) in self.crossings and (x, j + 1) in self.crossings:
                x = x - 1
            if x == 0 or (x, j) in self.crossings or (x, j + 1) in self.crossings:
                continue
            yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

    def upper_ladder_interval(self):
        level = {self}
        while level:
            new_level = set()
            for dream in level:
                yield dream
                for e in dream.ladder_moves():
                    new_level.add(e)
            level = new_level

    def is_symmetric(self):
        return all((j, i) in self.crossings for (i, j) in self.crossings)

    def lower_part(self):
        return Pipedream({(i, j) for (i, j) in self.crossings if i >= j})

    def strict_lower_part(self):
        return Pipedream({(i, j) for (i, j) in self.crossings if i > j})

    def weight(self):
        ans = []
        for (i, j) in self.crossings:
            while i - 1 >= len(ans):
                ans += [0]
            ans[i - 1] += 1
        return tuple(ans)

    def monomial(self):
        ans = one_var()
        for (i, j) in self.crossings:
            ans *= x_var(i)
        return ans

    def inv_monomial(self):
        ans = one_var()
        for (i, j) in self.crossings:
            if i == j:
                ans *= x_var(i)
            elif i > j:
                ans *= x_var(i) + x_var(j)
        return ans

    def fpf_monomial(self):
        ans = one_var()
        for (i, j) in self.crossings:
            if i > j:
                ans *= x_var(i) + x_var(j)
        return ans

    def tikz_simple(self):
        s = []
        s += ['\\begin{center}']
        s += ['\\begin{tikzpicture}[scale=0.5]']
        for i in range(1, self.n + 1):
            i_ = self.n - i
            s += ['\\node at (0, %s) {%s};' % (i_, i)]
            s += ['\\node at (%s, %s) {%s};' % (i, self.n, i)]
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if i + j > self.n + 1:
                    continue
                j_ = self.n - j
                if (j, i) not in self.crossings:
                    s += ['\\node at (%s, %s) {$\\cdot$};' % (i, j_)]
                if (j, i) in self.crossings:
                    s += ['\\node at (%s, %s) {$+$};' % (i, j_)]
        s += ['\\end{tikzpicture}']
        s += ['\\end{center}']
        return '\n'.join(s)

    def tikz(self):
        s = []
        s += ['\\begin{center}']
        s += ['\\begin{tikzpicture}']
        for i in range(1, self.n + 1):
            i_ = self.n - i
            s += ['\\node at (0, %s) {%s};' % (i_, i)]
            s += ['\\node at (%s, %s) {%s};' % (i, self.n, i)]
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if i + j > self.n + 1:
                    continue
                j_ = self.n - j
                if (j, i) not in self.crossings:
                    # s += ['\\draw [domain=270:360] plot ({%s-0.5+0.5* cos(\\x)}, {%s+0.5+0.5*sin(\\x)});' % (i, j_)]
                    s += ['\\draw (%s-0.5,%s) -- (%s,%s+0.5);' % (i, j_, i, j_)]
                if (j, i) not in self.crossings and i + j <= self.n:
                    # s += ['\\draw [domain=90:180] plot ({%s+0.5+0.5*cos(\\x)}, {%s-0.5+0.5*sin(\\x)});' % (i, j_)]
                    s += ['\\draw (%s+0.5,%s) -- (%s,%s-0.5);' % (i, j_, i, j_)]
                if (j, i) in self.crossings:
                    s += [
                        '\\draw (%s-0.5,%s) -- (%s+0.5,%s);' % (i, j_, i, j_),
                        '\\draw (%s,%s-0.5) -- (%s,%s+0.5);' % (i, j_, i, j_),
                    ]
        s += ['\\end{tikzpicture}']
        s += ['\\end{center}']
        return '\n'.join(s)

    def word(self):
        word = []
        for i in range(1, self.n + 1):
            for j in range(self.n + 1 - i, 0, -1):
                if (i, j) in self.crossings:
                    word += [j + i - 1]
        return tuple(word)

    def words(self):
        words = []
        for i in range(1, self.n + 1):
            word = []
            for j in range(self.n + 1 - i, 0, -1):
                if (i, j) in self.crossings:
                    word += [j + i - 1]
            words += [tuple(word)]
        return tuple(words)

    def decreasing_words(self, n=None):
        n = self.n if n is None else n
        words = []
        for i in range(1, self.n + 1):
            word = []
            for j in range(self.n + 1 - i, 0, -1):
                if (i, j) in self.crossings:
                    word += [n - (j + i - 1)]
            words += [tuple(word)]
        return tuple(words)

    def reverse_row_reading_words(self):
        words = []
        for i in range(self.n, 0, -1):
            word = []
            for j in range(1, self.n + 2 - i):
                if (i, j) in self.crossings:
                    word += [j + i - 1]
            words += [tuple(word)]
        return tuple(words)

    def unimodal_reading_words(self):
        words = []
        for i in range(2, self.n + 1):
            word = []
            j = 0
            while i + j <= self.n:
                a = i + j
                b = j + 1
                if (a, b) in self.crossings:
                    word += [a + b - 1]
                j += 1
            words += [tuple(word)]
        return tuple(words)

    def column_reading_words(self):
        words = []
        for j in range(self.n, 0, -1):
            word = []
            for i in range(1, self.n + 1):
                if (i, j) in self.crossings:
                    word += [i + j - 1]
            words += [tuple(word)]
        return tuple(words)

    def filename_prefix(self):
        word = []
        for i in range(1, self.n + 1):
            subword = []
            for j in range(self.n + 1 - i, 0, -1):
                if (i, j) in self.crossings:
                    subword += [str(j + i - 1)]
            word += [','.join(subword)]
        while word and word[-1] == '':
            word = word[:-1]
        return ':'.join(word)

    def filename(self):
        return self.filename_prefix() + '.tex'

    def save(self, directory=None):
        if directory is None:
            directory = self.DIRECTORY
        lines = [
            '\\documentclass[11pt]{article}',
            '\\usepackage{tikz}',
            '\\begin{document}',
            '\\title{}',
            '\\date{}',
            '\\maketitle',
            self.tikz(),
            '\\end{document}',
        ]
        file = '\n'.join(lines)
        with open(self.filename(), 'w') as f:
            f.write(file)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(["pdflatex", self.filename()], stdout=devnull)
            subprocess.run(["mv", self.filename(), directory], stdout=devnull)
            subprocess.run(["mv", self.filename_prefix() + '.pdf', directory], stdout=devnull)
            subprocess.run(["rm", self.filename_prefix() + '.log'], stdout=devnull)
            subprocess.run(["rm", self.filename_prefix() + '.aux'], stdout=devnull)

    @classmethod
    def save_involution(cls, perm):
        if len(perm) == 0:
            return
        dreams = perm.get_involution_pipe_dreams()
        test_s = sum([dream.inv_monomial() for dream in dreams])
        s = InvSchubert.get(perm)

        filename = ''.join(map(str, perm.oneline))
        directory = cls.DIRECTORY + 'inv/' + (test_s != s) * '_' + str(len(dreams)) + '_' + filename + '/'
        with open(os.devnull, 'w') as devnull:
            subprocess.run(["mkdir", directory], stdout=devnull)
        lines = [
            '\\documentclass[11pt]{article}',
            '\\usepackage{tikz}',
            '\\usepackage{amsfonts}'
            '\\begin{document}',
            '\\title{}',
            '\\date{}',
            '\\maketitle',
        ]
        lines += [e.tikz() for e in dreams]
        lines += ['']
        lines += ['\\[w = %s\\]' % perm.cycle_repr()]

        lines += ['\\[\\hat D(w) = \\{%s\\}\\]' % ', '.join([str(ij) for ij in perm.rothe_diagram() if ij[0] >= ij[1]])]
        lines += ['\\[\\hat\\mathfrak{S}_w = %s\\]' % s]
        lines += ['\\[\\hat\\mathfrak{S}_w - \\mathfrak{P}_w = %s\\]' % (s - test_s)]
        lines += ['\\end{document}']
        file = '\n'.join(lines)
        with open(filename + '.tex', 'w') as f:
            f.write(file)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(["pdflatex", filename + '.tex'], stdout=devnull)
            subprocess.run(["mv", filename + '.tex', directory], stdout=devnull)
            subprocess.run(["mv", filename + '.pdf', directory], stdout=devnull)
            subprocess.run(["rm", filename + '.log'], stdout=devnull)
            subprocess.run(["rm", filename + '.aux'], stdout=devnull)

    @classmethod
    def save_fpf_involution(cls, perm):
        if len(perm) == 0:
            return
        dreams = perm.get_fpf_pipe_dreams()
        test_s = sum([dream.fpf_monomial() for dream in dreams])
        s = FPFSchubert.get(perm)

        filename = ''.join(map(str, perm.oneline))
        directory = cls.DIRECTORY + 'fpf/' + (s != test_s) * '_' + str(len(dreams)) + '_' + filename + '/'
        with open(os.devnull, 'w') as devnull:
            subprocess.run(["mkdir", directory], stdout=devnull)
        lines = [
            '\\documentclass[11pt]{article}',
            '\\usepackage{tikz}',
            '\\usepackage{amsfonts}'
            '\\begin{document}',
            '\\title{}',
            '\\date{}',
            '\\maketitle',
        ]
        lines += [e.tikz_simple() for e in dreams]
        lines += ['']
        lines += ['\\[w = %s\\]' % perm.cycle_repr()]
        lines += ['\\[\\hat D_{FPF}(w) = \\{%s\\}\\]' % ', '.join([str(ij) for ij in perm.rothe_diagram() if ij[0] > ij[1]])]
        lines += ['\\[\\hat\\mathfrak{S}^{FPF}_w = %s\\]' % s]
        lines += ['\\[\\hat\\mathfrak{S}^{FPF}_w - \\mathfrak{P}_w = %s\\]' % (s - test_s)]
        lines += ['\\end{document}']
        file = '\n'.join(lines)
        with open(filename + '.tex', 'w') as f:
            f.write(file)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(["pdflatex", filename + '.tex'], stdout=devnull)
            subprocess.run(["mv", filename + '.tex', directory], stdout=devnull)
            subprocess.run(["mv", filename + '.pdf', directory], stdout=devnull)
            subprocess.run(["rm", filename + '.log'], stdout=devnull)
            subprocess.run(["rm", filename + '.aux'], stdout=devnull)
