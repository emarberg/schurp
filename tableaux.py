from cached import cached_value
from partitions import Shape, Partition, StrictPartition
from marked import MarkedNumber
import random
from collections import defaultdict


STANDARD_CACHE = {}
SEMISTANDARD_CACHE = {}
STANDARD_SHIFTED_MARKED_CACHE = {}
SEMISTANDARD_MARKED_RPP_CACHE = {}

HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_RPP_HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_RPP_VERTICAL_STRIPS_CACHE = {}

# for French notation
FRENCH = True


class Tableau:
    def __init__(self, dictionary=None):
        if dictionary is None:
            dictionary = {}

        self.mapping = {}
        for key, value in dictionary.items():
            if type(value) == int:
                value = MarkedNumber(value)
            self.mapping[key] = value
        self.max_row = max({i for i, j in self.mapping} | {0})
        self.max_column = max({j for i, j in self.mapping} | {0})

    def __iter__(self):
        return self.mapping.__iter__()

    def __getitem__(self, item):
        return self.mapping.get(item, None)

    def tex(self):
        rows = []
        for i in range(1, self.max_row + 1):
            row = []
            for j in range(1, self.max_column + 1):
                v = self.entry(i, j)
                row += [('*(white) ' + str(v)) if v is not None else '\\none']
            rows += [' & '.join(row)]
        return '$\\colorbox{lightgray!50}{\\begin{ytableau}' + ' \\\\ '.join(reversed(rows)) + '\\end{ytableau}}$'

    def is_shifted_k_flagged(self, k):
        for (i, j) in self:
            entry = self.entry(i, j)
            if entry.is_primed() and abs(entry) > j + k:
                return False
            if not entry.is_primed() and abs(entry) > i + k:
                return False
        return True

    def is_k_flagged(self, k):
        for (i, j) in self:
            entry = self.entry(i, j)
            assert not entry.is_primed()
            if abs(entry) > i + k:
                return False
        return True

    def unprime_diagonal(self):
        return Tableau({box: -val if box[0] == box[1] and val.is_primed() else val for (box, val) in self.mapping.items()})

    def unprime_two_diagonals(self):
        return Tableau({box: -val if box[1] <= box[0] + 1 and val.is_primed() else val for (box, val) in self.mapping.items()})

    def is_togglable(self, box):
        val = self.mapping[box]
        return not any(abs(v) == abs(val) and (b[0] > box[0] or b[1] < box[1]) for (b, v) in self.mapping.items())
        # return not any(abs(v) == abs(val) and b[0] >= box[0] and b[1] <= box[1] and b != box for (b, v) in self.mapping.items())

    def unprime_togglable(self):
        return Tableau({
            box: -val
            if val.is_primed() and
            self.is_togglable(box)
            # ((box[0] + 1, box[1]) not in self or abs(self.entry(box[0] + 1, box[1])) != abs(val)) and
            # ((box[0], box[1] - 1) not in self or abs(self.entry(box[0], box[1] - 1)) != abs(val))
            else val for (box, val) in self.mapping.items()})

    def unprime(self):
        return Tableau({box: -val if val.is_primed() else val for (box, val) in self.mapping.items()})

    @classmethod
    def shifted_from_rows(cls, rows):
        return cls.from_rows(rows, True)

    @classmethod
    def from_rows(cls, rows, shifted=False):
        mapping = {}
        for i in range(len(rows)):
            for j in range(len(rows[i])):
                mapping[(i + 1, j + 1 + (i if shifted else 0))] = rows[i][j]
        return Tableau(mapping)

    def standardize(self):
        values = sorted([(i, j, self[(i, j)]) for (i, j) in self], key=lambda x: (x[2], -x[0], x[1]) if not x[2].is_marked() else (x[2], -x[1], x[0]))
        mapping = {}
        for e, triple in enumerate(values):
            i, j, v = triple
            mapping[(i, j)] = e + 1 if not v.is_marked() else -e - 1
        return Tableau(mapping)

    def destandardize(self, alpha=None):
        if alpha is None:
            des = [0] + sorted(self.descent_set()) + [len(self)]
            alpha = [des[i] - des[i - 1] for i in range(1, len(des))]
        while alpha and alpha[-1] == 0:
            alpha = alpha[:-1]
        partialsums = [0]
        for a in alpha:
            partialsums.append(partialsums[-1] + a)
        des = self.descent_set()
        assert all(a in partialsums for a in des)
        assert sum(alpha) == len(self)
        conversion = {}
        for i in range(1, len(partialsums)):
            for a in range(partialsums[i - 1] + 1, partialsums[i] + 1):
                conversion[MarkedNumber(a)] = i
                conversion[MarkedNumber(-a)] = -i
        return Tableau({ij: conversion[self[ij]] for ij in self})

    def dual_equivalence_operator(self, index, jndex=None):
        jndex = index if jndex is None else jndex
        ans = self
        for i in range(index, jndex + 1) if index <= jndex else range(index, jndex - 1, -1):
            ans = ans._dual_equivalence_operator(i)
        return ans

    def _dual_equivalence_operator(self, index):
        tab = self

        def find(t, a):
            for i, j in t:
                if abs(t.get(i, j)) == a:
                    return i, j

        def locations(t, x):
            a, b, c = None, None, None
            w = t.shifted_reading_word()
            for i in range(len(w)):
                a = i if abs(w[i]) == x else a
                b = i if abs(w[i]) == (x + 1) else b
                c = i if abs(w[i]) == (x + 2) else c
            return a, b, c

        if index == 0:
            i, j = find(tab, 2)
            return tab.set(i, j, -tab.get(i, j))

        a, b, c = locations(tab, index)

        if a is None or b is None or c is None or a < b < c or c < b < a:
            return tab

        i1, j1 = find(tab, index)
        i2, j2 = find(tab, index + 1)
        i3, j3 = find(tab, index + 2)

        x1 = tab.get(i1, j1)
        x2 = tab.get(i2, j2)
        x3 = tab.get(i3, j3)

        if i1 == j1 and i3 == j3 and x1.number * x3.number < 0:
            return tab.set(i1, j1, -x1).set(i2, j2, -x2).set(i3, j3, -x3)

        if b < a < c or c < a < b:
            i1, j1, i2, j2 = i2, j2, i3, j3

        x1, x2 = tab.get(i1, j1), tab.get(i2, j2)

        if i1 == i2 or j1 == j2:
            tab = tab.set(i1, j1, -x1 if i1 != j1 else x1)
            tab = tab.set(i2, j2, -x2 if i2 != j2 else x2)
        else:
            tab = tab.set(i1, j1, x1 - 1 if x1.number < 0 else x1 + 1)
            tab = tab.set(i2, j2, x2 + 1 if x2.number < 0 else x2 - 1)

        return tab

    def shifted_crystal_word(self):
        rows = [
            (i, j, self.mapping[(i, j)])
            for (i, j) in sorted(self.mapping, key=lambda x:(-x[0], x[1]))
        ]
        cols = [
            (i, j, self.mapping[(i, j)])
            for (i, j) in sorted(self.mapping, key=lambda x:(-x[1], x[0]))
        ]

        word, positions = [], []

        # for i, j, v in cols:
        #     if v.is_marked():
        #         word.append(v)
        #         positions.append((i, j))
        # for i, j, v in rows:
        #     if not v.is_marked():
        #         word.append(v)
        #         positions.append((i, j))

        # return word, positions

        a, b = 0, 0
        for t in range(max(self.max_row, self.max_column), 0, -1):
            while a < len(cols) and cols[a][1] == t:
                (i, j, v) = cols[a]
                if v.is_marked():
                    word.append(v)
                    positions.append((i, j))
                a += 1
            while b < len(rows) and rows[b][0] == t:
                (i, j, v) = rows[b]
                if not v.is_marked():
                    word.append(v)
                    positions.append((i, j))
                b += 1
        assert len(word) == len(positions) == len(self)
        return word, positions

    def shifted_crystal_s(self, index):
        weight = {i + 1: a for i, a in enumerate(self.weight())}
        k = weight.get(index, 0) - weight.get(index + 1, 0)
        if k == 0:
            return self
        ans = self
        if k < 0:
            for _ in range(-k):
                ans = ans.shifted_crystal_e(index)
        else:
            for _ in range(k):
                ans = ans.shifted_crystal_f(index)
        return ans

    def extended_shifted_crystal_e0(self, i):
        ans = self
        for j in range(i - 1, 0, -1):
            ans = ans.shifted_crystal_s(j)
        ans = ans.shifted_crystal_e(0)
        if ans is not None:
            for j in range(1, i):
                ans = ans.shifted_crystal_s(j)
        return ans

    def extended_shifted_crystal_f0(self, i):
        ans = self
        for j in range(i - 1, 0, -1):
            ans = ans.shifted_crystal_s(j)
        ans = ans.shifted_crystal_f(0)
        if ans is not None:
            for j in range(1, i):
                ans = ans.shifted_crystal_s(j)
        return ans

    def shifted_crystal_e(self, index, verbose=False):
        if index == 0:
            if (1, 1) not in self or self[(1, 1)] != MarkedNumber(-1):
                return None
            return self.set(1, 1, 1)

        if index == -1:
            twos = [cell for cell in self if abs(self[cell]) == 2 and cell[0] == 1]
            if len(twos) == 0:
                return None
            a, b = min(twos)
            if b > 1 and self[(a, b)] == MarkedNumber(2):
                return None
            v = MarkedNumber(-1) if (b == 1 and self[(a, b)] == MarkedNumber(-2)) else MarkedNumber(1)
            return self.set(a, b, v)

        word, positions = self.shifted_crystal_word()

        p, queue = None, []
        for i in reversed(range(len(word))):
            v = abs(word[i])
            if v == index:
                queue.append(i)
            elif v == index + 1 and queue:
                queue.pop()
            elif v == index + 1:
                p = i

        # print('p:', word, p, 'index:', index)

        if p is None:
            return

        x, (a, b) = word[p], positions[p]
        y = self[(a - 1, b)]  # could be None
        z = self[(a, b - 1)]  # could be None

        # z x
        # ? y

        if verbose:
            print('x:', x, '(a, b):', (a, b), '\n')

        if not x.is_marked():
            # z cannot be i+1 since x is first unpaired
            assert z is None or z != MarkedNumber(index + 1)
            # ? cannot be i since x unpaired
            assert (a - 1, b - 1) not in self or self[(a - 1, b - 1)] < MarkedNumber(index)

            if z is not None and z.number == -index - 1:
                # (i+1)' i+1 ->  i (i+1)'
                # ?      y   ->  ? y
                #
                # ? cannot be i by unpairedness
                assert (a - 1, b - 1) not in self or self[(a - 1, b - 1)] != MarkedNumber(index)
                # y if nonempty is at most (i+1)'
                assert (a - 1, b) not in self or self[(a - 1, b)] <= MarkedNumber(-index - 1)
                # z=(i+1)' cannot be on main diagonal by firstness
                assert a != b - 1
                if verbose:
                    print('\n* case R1(a)\n')
                return self.set(a, b, z).set(a, b - 1, index)
            if y is None or y <= MarkedNumber(-index):
                # z i+1 -> z i
                # ? y   -> ? y
                #
                # z cannot be (i+1)' as not in previous case or i+1 since x is first unpaired
                # so z <= i
                assert z is None or z <= MarkedNumber(index)
                #
                if verbose:
                    print('\n* case R1(b)\n')
                return self.set(a, b, index)

            # z x = z i+1 -> z i+1' = z H -> z H~
            # ? y = ? y   -> ? y    = ? y -> ? y
            #
            # z <= i as z cannot be i+1' or i+1 by firstness+previous cases
            assert z is None or z <= MarkedNumber(index)
            # y must be i' or i or (i+1)' by previous case
            assert y in [MarkedNumber(index), MarkedNumber(-index - 1)]
            #
            # if y=i then the column [ i+1 i ]^T
            # in positions [ x y ]^T must be repeated rightwards
            # until there is a column [ i+1 i+1' ]^T by unpairedness
            #
            # change x from i+1 to i+1' then let H the lower right (i+1)-ribbon
            # starting at x. Change to i the first i+1' in H that can be
            # changed to i to form H~. This i+1' might occur in the position of x.
            #
            # some such i+1' in H must exist by unpairedness

            ans = self.set(a, b, -index - 1)
            rx, ry = a, b
            while True:
                if ans[(rx, ry)].is_marked() and ((rx - 1, ry) not in ans or ans[(rx - 1, ry)] not in [MarkedNumber(index), MarkedNumber(-index - 1)]):
                    ans = ans.set(rx, ry, index)
                    break
                if (rx - 1, ry) in ans and abs(ans[(rx - 1, ry)]) == index + 1:
                    rx -= 1
                elif (rx, ry + 1) in ans and abs(ans[(rx, ry + 1)]) == index + 1:
                    ry += 1
                else:
                    assert False

            if verbose:
                print('\n* case R1(c)\n  (rx, ry) =', rx, ry, '\n  (a, b) =', a, b, '\n')
                # Assaf-Oguz error: "southwest" -> "southeast"
            return ans

        # z x
        # ? y

        else:
            # y cannot be (i+1)' since x is first unpaired
            assert y is None or y != MarkedNumber(-index - 1)
            # ? cannot be i' or i
            assert (a - 1, b - 1) not in self or abs(self[(a - 1, b - 1)]) < index

            if y is not None and y.number == index:
                # z (i+1)' ->  z i
                # ? i      ->  ? i'
                #
                # z <= i by semistandardness
                # ? cannot be i' by unpairedness
                assert (a - 1, b - 1) not in self or self[(a - 1, b - 1)] != MarkedNumber(-index)
                # x cannot be on the diagonal by firstness+unpairedness (?)
                assert a != b
                #
                if verbose:
                    print('\n* case R2(a)\n')
                return self.set(a, b, index).set(a - 1, b, -index)
            if z is None or abs(z) < index:
                # z (i+1)' ->  z i'
                # ? y      ->  ? y
                #
                # y cannot be (i+1)'
                # y cannot be i as not in previous case
                # so y <= i'
                assert y is None or y <= MarkedNumber(-index)
                #
                if verbose:
                    print('\n* case R2(b)\n')
                return self.set(a, b, -index)

            # H (i+1)' -> H~ i
            # ? y      -> ?  y
            #
            # y <= i' since y cannot be (i+1)' or i by previous cases
            assert y is None or y <= MarkedNumber(-index)
            # z must be i or i' by semistandardness and not in previous case
            assert abs(z) == index
            # thus H is a nonempty i-ribbon
            #
            # If H ends on diagonal then H~ = H, unless next diagonal
            # box has opposite prime and contains i+1 or i+1'.
            #
            # Otherwise H ends in i off the diagonal by unpairedness,
            # and H~ is formed by add a prime to this entry.

            rx, ry = a, b - 1
            while True:
                if self[(rx, ry)].is_marked():
                    if (rx + 1, ry) not in self or abs(self[(rx + 1, ry)]) != index:
                        break
                    rx += 1
                else:
                    if (rx, ry - 1) not in self or abs(self[(rx, ry - 1)]) != index:
                        break
                    ry -= 1

            if verbose:
                print('\n* case R2(c)\n  (rx, ry) =', rx, ry, '\n  (a, b) =', a, b, '\n')
                # Assaf-Oguz error:
                #   the right side of the picture in Figure 18 is wrong, has an extra cell
                # Assaf-Oguz error:
                #   "northeastern" -> "northwestern"
                # Assaf-Oguz error:
                #   x should change to i (as shown in Figure 18) not i' (as stated)

            ans = self
            if rx != ry:
                assert not self[(rx, ry)].is_marked()
                ans = ans.set(rx, ry, -index)
            elif rx == ry and self[(rx, ry)].is_marked() and self[(rx + 1, ry + 1)] == MarkedNumber(index + 1):
                ans = ans.set(rx, ry, index)
                ans = ans.set(rx + 1, ry + 1, -index - 1)
            elif rx == ry and not self[(rx, ry)].is_marked() and self[(rx + 1, ry + 1)] == MarkedNumber(-index - 1):
                ans = ans.set(rx, ry, -index)
                ans = ans.set(rx + 1, ry + 1, index + 1)

            return ans.set(a, b, index)

    def shifted_crystal_last_unpaired_box(self, index):
        word, positions = self.shifted_crystal_word()

        p, queue = None, []
        for i in range(len(word)):
            v = abs(word[i])
            if v == index + 1:
                queue.append(i)
            elif v == index and queue:
                queue.pop()
            elif v == index:
                p = i

        if p is None:
            return None

        x, (a, b) = word[p], positions[p]
        return x, (a, b)

    def shifted_crystal_f(self, index, verbose=False):
        if index == 0:
            if (1, 1) not in self or self[(1, 1)] != MarkedNumber(1):
                return None
            return self.set(1, 1, -1)

        if index == -1:
            ones = [cell for cell in self if abs(self[cell]) == 1]
            if len(ones) == 0:
                return None
            if len([cell for cell in self if self[cell] == MarkedNumber(-2) and cell[0] == 1]) > 0:
                return None
            a, b = max(ones)
            assert a == 1
            v = MarkedNumber(2) if (b == 1 and self[(a, b)] == MarkedNumber(1)) else MarkedNumber(-2)
            return self.set(a, b, v)

        unpaired = self.shifted_crystal_last_unpaired_box(index)
        if unpaired is None:
            return
        else:
            x, (a, b) = unpaired
        y = self[(a + 1, b)]  # could be None
        z = self[(a, b + 1)]  # could be None

        if verbose:
            print('x:', x, '(a, b):', (a, b), '\n')

        # y
        # x z

        if not x.is_marked():
            # z cannot be i since x is last unpaired
            assert z is None or z != MarkedNumber(index)
            assert (a + 1, b + 1) not in self or self[(a + 1, b + 1)] > MarkedNumber(index + 1)

            if z is not None and z.number == -index - 1:
                # ? empty     ->  ?      empty
                # i (i+1)'    ->  (i+1)' i+1
                #
                # i cannot be on main diagonal ny unpairedness
                assert a != b
                #
                # ? must be empty or (i+1)' or i+1
                assert (a + 1, b) not in self or abs(self[(a + 1, b)]) >= index + 1
                # ??:=empty must be empty because:
                assert (a + 1, b + 1) not in self or abs(self[(a + 1, b + 1)]) > index + 1
                #
                # (1) if ? is empty then ?? must also be empty
                #
                # (2) if ? is not empty then each i to the
                # left of x is below (i+1)' or i+1 and
                # so ?? must be empty as otherwise ?? would
                # have to be i+1 so i in x would be matched
                if verbose:
                    print('\n* case L1(a)\n')
                return self.set(a, b, z).set(a, b + 1, -z)

            if y is None or abs(y) > index + 1:
                # empty empty   ->  empty empty
                # i     ?       ->  i+1   ?
                #
                # ? cannot be i' or i or, by previous case, (i+1)'
                assert y is None or y not in [MarkedNumber(-index), MarkedNumber(index), MarkedNumber(-index - 1)]
                if verbose:
                    print('\n* case L1(b)\n')
                return self.set(a, b, index + 1)

            # H empty   ->  H~     empty
            # i ?       ->  (i+1)' ?
            #
            # ? is i+1 or empty
            assert z is None or z == MarkedNumber(index + 1) or abs(z) > index + 1
            # H is a nonempty (i+1)-ribbon
            # the empty position must be empty by (1)-(2) above
            assert (a + 1, b + 1) not in self or abs(self[(a + 1, b + 1)]) > index + 1
            #
            # If H ends on diagonal then H~ = H, unless previous diagonal
            # box has opposite prime and contains i or i'.
            #
            # Otherwise H ends in (i+1)' off the diagonal by unpairedness,
            # and H~ is formed by removing the prime from this entry.

            rx, ry = a + 1, b
            while True:
                if self[(rx, ry)].is_marked():
                    if (rx + 1, ry) not in self or abs(self[(rx + 1, ry)]) != index + 1:
                        break
                    rx += 1
                else:
                    if (rx, ry - 1) not in self or abs(self[(rx, ry - 1)]) != index + 1:
                        break
                    ry -= 1

            ans = self
            if rx != ry:
                assert self[(rx, ry)].is_marked()
                ans = ans.set(rx, ry, index + 1)
            elif rx == ry and self[(rx, ry)].is_marked() and self[(rx - 1, ry - 1)] == MarkedNumber(index):
                ans = ans.set(rx, ry, index + 1)
                ans = ans.set(rx - 1, ry - 1, -index)
            elif rx == ry and not self[(rx, ry)].is_marked() and self[(rx - 1, ry - 1)] == MarkedNumber(-index):
                ans = ans.set(rx, ry, -index - 1)
                ans = ans.set(rx - 1, ry - 1, index)

            if verbose:
                print('\n* case L1(c)\n')
                # Assaf-Oguz error:
                #   "northeastern" -> "northwestern"
            return ans.set(a, b, -index - 1)

        else:
            # y cannot be i' since x is last unpaired
            assert y is None or y != MarkedNumber(-index)
            assert (a + 1, b + 1) not in self or self[(a + 1, b + 1)] >= MarkedNumber(index + 1)

            if y is not None and y.number == index:
                # i  ?   ->  (i+1)' ?
                # i' ??  ->  i      ??
                #
                # ? cannot be (i+1)' by unpairedness
                assert (a + 1, b + 1) not in self or self[(a + 1, b + 1)] != MarkedNumber(-index - 1)
                if verbose:
                    print('\n* case L2(a)\n')
                return self.set(a, b, index).set(a + 1, b, -index - 1)
            if z is None or z.number == index + 1 or abs(z) >= index + 2:
                # ?  empty  ->  ?      empty
                # i' ??     ->  (i+1)' ??
                #
                # ? cannot be i' by lastness or i by previous case
                assert y is None or y not in [MarkedNumber(index), MarkedNumber(-index)]
                # ?? must be empty or i+1 by assumption
                if verbose:
                    print('\n* case L2(b)\n')
                return self.set(a, b, -index - 1)

            # ?  ??  ->  ? ??  -> ?  ??
            # i' ??? ->  H ??? -> H~ ???
            #
            assert a != b
            # ? cannot be i' or i by lastness+previous cases
            assert y is None or abs(y) != index
            # ?? cannot be (i+1)' by unpairedness
            assert (a + 1, b + 1) not in self or self[(a + 1, b + 1)] != MarkedNumber(-index - 1)
            # ??? must be i or (i+1)' by previous case
            assert z in [MarkedNumber(index), MarkedNumber(-index - 1)]
            #
            # if ??? is (i+1)' then the row i' (i+1)'
            # in positions x z must be repeated below
            # until there is a row i' i, by unpairedness
            #
            # change the i' in x to i then let H the lower right i-rim hook
            # starting at x, then change the first i in H that can be
            # changed to (i+1)' to form H~
            #
            # some such i in H must exist by unpairedness

            ans = self
            rx, ry = a, b
            while True:
                if not self[(rx, ry)].is_marked() and ((rx, ry + 1) not in self or self[(rx, ry + 1)] not in [MarkedNumber(index), MarkedNumber(-index - 1)]):
                    ans = ans.set(rx, ry, -index - 1)
                    break
                if (rx - 1, ry) in self and abs(self[(rx - 1, ry)]) == index:
                    rx -= 1
                elif (rx, ry + 1) in self and abs(self[(rx, ry + 1)]) == index:
                    ry += 1
                else:
                    assert False

            if verbose:
                print('\n* case L2(c)\n')
                # Assaf-Oguz error:
                #   "southwest" -> "southeast"
            return ans.set(a, b, index)

    def restrict(self, n):
        n = MarkedNumber(n) if type(n) == int else n
        assert type(n) == MarkedNumber
        return Tableau({k: v for k, v in self.mapping.items() if v <= n})

    def durfee(self):
        i = 0
        while (i + 1, i + 1) in self.mapping:
            i += 1
        return i

    @classmethod
    def from_composition(cls, alpha):
        ans = Tableau()
        for i, a in enumerate(alpha):
            for col in range(1, a + 1):
                row = 1
                while (row, col) in ans:
                    row += 1
                ans = ans.set(row, col, i + 1)
        return ans

    def weight(self, as_dict=False):
        ans = {}
        for v in self.mapping.values():
            v = v.weight()
            if v not in ans:
                ans[v] = 0
            ans[v] += 1
        if as_dict:
            return ans
        assert all(type(m) == int and m > 0 for m in ans)
        m = max(ans) if ans else 0
        alpha = m * [0]
        for v in ans:
            alpha[v - 1] = ans[v]
        return tuple(alpha)

    def shape(self):
        return Shape(self.mapping.keys())

    def composition(self):
        return tuple(len(_) for _ in self.get_rows())

    def partition(self):
        rows = defaultdict(int)
        for i, j in self.mapping:
            rows[i] += 1
        return Partition(*sorted(rows.values(), reverse=True))

    @classmethod
    def shifted_highest_weight(cls, mu):
        mapping = {}
        for i in range(len(mu)):
            for j in range(mu[i]):
                mapping[i + 1, i + j + 1] = i + 1
        return cls(mapping)

    @classmethod
    def shifted_lowest_weight(cls, mu, rank=None, diagonal_primes=False):
        rank = len(mu) if rank is None else rank
        if sum(mu) == 0:
            return Tableau()
        else:
            boxes = {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(mu[i])}
            # find start
            j = 1
            while (1, j + 1) in boxes:
                j += 1
            i = 1
            # find outer ribbon
            ribbon = []
            while True:
                if (i + 1, j) in boxes:
                    ribbon.append((i, j, -1))
                    i += 1
                elif (i, j - 1) in boxes:
                    ribbon.append((i, j, 1))
                    j -= 1
                else:
                    ribbon.append((i, j, -1 if diagonal_primes else 1))
                    break
            # trim outer ribbon
            nu = list(mu)
            for i, j, _ in ribbon:
                nu[i - 1] -= 1
            nu = Partition.trim(nu)
            # get inner part of tableau and add outer ribbon
            mapping = cls.shifted_lowest_weight(nu, rank - 1, diagonal_primes).mapping
            for (i, j, sign) in ribbon:
                mapping[i, j] = sign * rank
            return cls(mapping)

    def count_diagonal_cells(self):
        return len([(i, j) for i, j in self.mapping if i == j])

    @classmethod
    def from_string(cls, string):
        def mark(i):
            i = i.strip()
            if i == "":
                return None
            if i.endswith("'"):
                return MarkedNumber(-int(i[:-1]))
            else:
                return MarkedNumber(int(i))
        rows = [[mark(i) for i in row.split(',')] for row in string.split(';')]
        dictionary = {
            (i + 1, j + 1): rows[i][j]
            for i in range(len(rows)) for j in range(len(rows[i]))
            if rows[i][j]
        }
        return Tableau(dictionary)

    @classmethod
    def from_partition(cls, mu):
        mapping = {(i, j): 1 for (i, j) in mu.shape}
        return Tableau(mapping)

    def add(self, i, j, v):
        return self.set(i, j, v)

    def set(self, i, j, v):
        mapping = self.mapping.copy()
        mapping[(i, j)] = v
        return Tableau(mapping)

    def toggle(self):
        subtableau = self.find(MarkedNumber(-2))
        if subtableau:
            assert len(subtableau) == 1
            i, j = max(subtableau.cells())
            return self.set(i, j, MarkedNumber(1))

        subtableau = self.find(MarkedNumber(1))
        if subtableau:
            i, j = max(subtableau.cells())
            if i == j:
                return self.set(i, j, MarkedNumber(2))
            else:
                return self.set(i, j, MarkedNumber(-2))

        subtableau = self.find(MarkedNumber(2))
        if subtableau:
            i, j = min(subtableau.mapping)
            return self.set(i, j, MarkedNumber(1))

        return self

    def find(self, *args):
        args = [MarkedNumber(v) if type(v) == int else v for v in args]
        return Tableau({key: value for key, value in self.mapping.items() if value in args})

    def shift(self):
        return Tableau({(i, i + j - 1): self.entry(i, j) for (i, j) in self.mapping})

    def translate_left(self, steps=1):
        return Tableau({(i, j - steps): self.entry(i, j) for i, j in self.mapping})

    def transpose(self):
        return Tableau({(j, i): self.entry(i, j) for i, j in self.mapping})

    def clean_mapping(self):
        return {(i, j): self.entry(i, j).number for (i, j) in self.mapping}

    def double(self, shift=False):
        assert self.is_shifted()
        mapping = {(i, j + shift): self.entry(i, j) for i, j in self.mapping}
        for i, j in self.mapping:
            if (j, j) not in mapping or shift:
                mapping[(j, i)] = self.entry(i, j)
        return Tableau(mapping)

    def fpf_double(self):
        assert self.is_shifted()
        mapping = {(i, j + 1): self.entry(i, j) for i, j in self.mapping}
        for i, j in self.mapping:
            mapping[(j + 1, i)] = self.entry(i, j)
        offset = 0 if (1, 1) not in self.mapping else (self.entry(1, 1).number - 2) // 2
        for i, j in list(mapping.keys()):
            if (i, i + 1) in mapping or (i, i - 1) in mapping:
                mapping[(i, i)] = 2 * (i + offset) - 1
        return Tableau(mapping)

    def halve(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i <= j}
        return Tableau(mapping)

    def lower_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i <= j}
        return Tableau(mapping)

    def strict_lower_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i < j}
        return Tableau(mapping)

    def upper_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i >= j}
        return Tableau(mapping)

    def strict_upper_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i > j}
        return Tableau(mapping)

    def maximum(self):
        if self.mapping:
            return max(self.entries())

    def cells(self):
        return self.mapping.keys()

    def entries(self):
        return self.mapping.values()

    def get(self, i, j):
        return self.entry(i, j)

    def entry(self, i, j):
        return self.mapping.get((i, j), None)

    def pop(self, i, j):
        assert (i, j) in self
        mapping = self.mapping.copy()
        del mapping[(i, j)]
        return self.entry(i, j), Tableau(mapping)

    def get_main_diagonal(self):
        return tuple(self.entry(i, i) for i in range(1, self.max_row + 1) if (i, i) in self.mapping)

    def get_row(self, i):
        columns = sorted([j for (i_, j) in self.mapping if i == i_])
        return tuple(self.entry(i, j) for j in columns)

    def get_column(self, j):
        rows = sorted([i for (i, j_) in self.mapping if j == j_])
        return tuple(self.entry(i, j) for i in rows)

    def get_hook(self, j, k):
        ans = {(x, y) for (x, y) in self if x == j and k <= y}
        ans |= {(x, y) for (x, y) in self if y == k and j < x}
        return ans

    def get_shifted_hook(self, j, k):
        assert self.is_shifted()
        ans = self.get_hook(j, k)
        if (k, k) in ans:
            ans |= {(x, y) for (x, y) in self if x == k + 1}
        return ans

    def hooks(self):
        return Tableau({(i, j): len(self.get_hook(i, j)) for (i, j) in self})

    def shifted_hooks(self):
        assert self.is_shifted()
        return Tableau({(i, j): len(self.get_shifted_hook(i, j)) for (i, j) in self})

    @classmethod
    def random_shifted(cls, mu):
        mu = StrictPartition(*mu.parts[:])
        nu = defaultdict(int)
        for i, j in mu.shape:
            nu[j] += 1
        nu = dict(nu)
        t = Tableau.from_partition(mu)

        for n in range(abs(mu), 0, -1):
            m = random.randint(1, abs(mu))
            i = 1
            while m > mu(i) and i <= len(mu):
                m -= mu(i)
                i += 1
            j = m + i - 1

            while True:
                row = mu(i) - j + i - 1
                col = nu[j] - i
                res = mu(j + 1)
                if row == 0 and col == 0 and res == 0:
                    break
                m = random.randint(1, row + col + res)
                if m <= row:
                    j += m
                elif row + col < m:
                    i = j + 1
                    j = j + m - row - col
                else:
                    i += m - row

            t.mapping[(i, j)] = MarkedNumber(-n if i != j and random.random() < 0.5 else n)
            mu.decrement(i)
            nu[j] -= 1
            if n % 1000 == 0:
                print(n, '. . .')
        return t

    @classmethod
    def inverse_fpf_insertion(cls, p, q):
        # q is assumed to be semistandard shifted with no diagonal primes
        order = sorted(q, key=lambda x: (q[x], x[0] if q[x].is_primed() else x[1]))
        alpha = []
        for i, j in q:
            v = abs(q[i, j])
            while len(alpha) < v:
                alpha.append(0)
            alpha[v - 1] += 1
        mapping = {}
        for i, (a, b) in enumerate(order):
            mapping[a, b] = -(i + 1) if q[a, b].is_primed() else (i + 1)
        qstandard = Tableau(mapping)
        word = cls.inverse_fpf(p, qstandard)
        ans = []
        for a in alpha:
            ans.append(word[:a])
            word = word[a:]
        return tuple(ans)


    @classmethod
    def inverse_fpf(cls, p, q):
        # q is assumed to be standard shifted with no diagonal primes
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}

        order = len(q) * [0]
        for (i, j) in q:
            order[abs(q[(i, j)]) - 1] = (i, j, q[(i, j)] < 0)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j, signed = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]
            if signed:
                for col in range(j - 1, 0, -1):
                    row = i
                    while p.get((row + 1, col), a) < a:
                        row += 1
                    if row == col:
                        i, j = col, col
                        if a % 2 == 0:
                            a, p[(i, j)] = p[(i, j)], a
                        else:
                            a -= 2
                        break
                    elif (row + 1, col) in p and p[(row + 1, col)] == a:
                        a = p[(row, col)]
                    else:
                        a, p[(row, col)] = p[(row, col)], a
            for row in range(i - 1, 0, -1):
                col = j
                while p.get((row, col + 1), a + 1) <= a:
                    col += 1
                if p[(row, col)] == a:
                    a = p[(row, col - 1)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def inverse_involution_insertion(cls, p, q):
        # q is assumed to be semistandard shifted with no diagonal primes
        order = sorted(q, key=lambda x: (q[x], x[0] if q[x].is_primed() else x[1]))
        alpha = []
        for i, j in q:
            v = abs(q[i, j])
            while len(alpha) < v:
                alpha.append(0)
            alpha[v - 1] += 1
        mapping = {}
        for i, (a, b) in enumerate(order):
            mapping[a, b] = -(i + 1) if q[a, b].is_primed() else (i + 1)
        qstandard = Tableau(mapping)
        word = cls.inverse_inv(p, qstandard)
        ans = []
        for a in alpha:
            ans.append(word[:a])
            word = word[a:]
        return tuple(ans)

    @classmethod
    def inverse_inv(cls, p, q):
        # q is assumed to be standard shifted with no diagonal primes
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}

        order = len(q) * [0]
        for (i, j) in q:
            order[abs(q[(i, j)]) - 1] = (i, j, q[(i, j)] < 0)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j, signed = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]
            if not signed:
                i = i - 1
            if signed:
                for col in range(j - 1, 0, -1):
                    row = i
                    while p.get((row + 1, col), a) < a:
                        row += 1
                    if row == col:
                        i, j = col, col
                        break
                    if (row + 1, col) in p and p[(row + 1, col)] == a:
                        a = p[(row, col)]
                    else:
                        a, p[(row, col)] = p[(row, col)], a
            for row in range(i, 0, -1):
                col = j
                while p.get((row, col + 1), a) < a:
                    col += 1
                if (row, col + 1) in p and p[(row, col + 1)] == a:
                    a = p[(row, col)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def random(cls, mu):
        mu = Partition(*mu.parts[:])
        nu = mu.transpose()
        t = Tableau.from_partition(mu)

        for n in range(abs(mu), 0, -1):
            m = random.randint(1, abs(mu))
            i = 1
            while m > mu(i) and i <= len(mu):
                m -= mu(i)
                i += 1
            j = m

            while True:
                row = mu(i) - j
                col = nu(j) - i

                if row == 0 and col == 0:
                    break
                m = random.randint(1, row + col)

                if m <= row:
                    j += m
                elif row < m:
                    i += m - row
            t.mapping[(i, j)] = MarkedNumber(n)
            mu.decrement(i)
            nu.decrement(j)
            if n % 1000 == 0:
                print(n, '. . .')
        return t

    @classmethod
    def inverse_eg(cls, p, q):
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}
        order = len(q) * [0]
        for (i, j) in q:
            order[q[(i, j)] - 1] = (i, j)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]

            for row in range(i - 1, 0, -1):
                col = j
                while p.get((row, col + 1), a + 1) <= a:
                    col += 1
                if p[(row, col)] == a:
                    a = p[(row, col - 1)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def longest_eg_insertion_tableau(cls, n):
        return Tableau.from_string(';'.join([','.join([str(j) for j in range(i, n)]) for i in range(1, n)]))

    @classmethod
    def longest_inv_insertion_tableau(cls, n):
        partition = StrictPartition(*list(range(n - 1, 0, -2)))
        t = Tableau.from_partition(partition)
        for i in range(1, len(partition) + 1):
            for j in range(i, n - i + 1):
                t.mapping[(i, j)] = MarkedNumber(i + j - 1)
        return t

    @classmethod
    def longest_fpf_insertion_tableau(cls, n):
        partition = StrictPartition(*list(range(n - 2, 0, -2)))
        t = Tableau.from_partition(partition)
        for i in range(1, len(partition) + 1):
            for j in range(i, n - i):
                t.mapping[(i, j)] = MarkedNumber(i + j)
        return t

    @classmethod
    def random_sorting_network(cls, n):
        p = Tableau.longest_eg_insertion_tableau(n)
        q = cls.random(p.partition())
        return cls.inverse_eg(p, q)

    @classmethod
    def random_inv_network(cls, n):
        p = Tableau.longest_inv_insertion_tableau(n)
        q = cls.random_shifted(StrictPartition(p.partition()))
        return cls.inverse_inv(p, q)

    @classmethod
    def random_fpf_network(cls, n):
        p = Tableau.longest_fpf_insertion_tableau(n)
        q = cls.random_shifted(StrictPartition(p.partition()))
        return cls.inverse_fpf(p, q)

    def is_key_tableau(self):
        rows = self.get_rows()
        # check that shape is key diagram with weakly increasing rows
        for index, r in enumerate(rows):
            if any(r[j] < r[j + 1] for j in range(len(r) - 1)):
                print('*1')
                return False
            if any(self.get(index + 1, j + 1) is None for j in range(len(r))):
                print('*2')
                return False
        n = max([0] + [len(r) for r in rows])
        m = len(rows)
        # check that columns have distinct entries + triangle condition
        for j in range(1, n + 1):
            for i1 in range(1, m + 1):
                a = self.get(i1, j)
                if a is None:
                    continue
                a = a.number
                for i2 in range(i1 + 1, m + 1): 
                    b = self.get(i2, j)
                    if b is None:
                        continue
                    b = b.number
                    if a == b:
                        print('*3')
                        return False
                    if a > b and not ((i1, j + 1) in self and self.get(i1, j + 1).number > b):
                        print('*4')
                        return False
        return True

    def is_key_flagged(self, flag=None):
        if flag is None:
            return self.is_key_tableau() and all(self.get(i, j).number <= i for i, j in self)
        else:
            return self.is_key_tableau() and all(self.get(i, j).number <= flag[i - 1] for i, j in self)

    def is_shifted(self):
        return not any(j < i for i, j in self.mapping)

    def is_contiguous(self):
        values = {v.weight() for v in self.mapping.values()}
        return values == set(range(1, len(values) + 1))

    def is_unmarked(self):
        return not any(v.is_marked() for v in self.mapping.values())

    def is_diagonally_unmarked(self):
        return not any(self.entry(i, i).is_marked() for i in range(1, self.max_row + 1))

    def is_weakly_increasing(self):
        for i, j in self.cells():
            u = self.entry(i, j)
            v = self.entry(i, j + 1)
            w = self.entry(i + 1, j)
            if v and u == v and u.is_marked():
                return False
            if w and u == w and not u.is_marked():
                return False
            if (v and v < u) or (w and w < u):
                return False
        return True

    def is_increasing(self):
        for i, j in self.cells():
            u = self.entry(i, j)
            v = self.entry(i, j + 1)
            w = self.entry(i + 1, j)
            if (v and v <= u) or (w and w <= u):
                return False
        return True

    def is_decreasing(self):
        for i, j in self.cells():
            u = self.entry(i, j)
            v = self.entry(i, j + 1)
            w = self.entry(i + 1, j)
            if (v and v >= u) or (w and w >= u):
                return False
        return True

    def is_semistandard(self):
        return self.is_unmarked() and self.is_weakly_increasing()

    def is_standard(self):
        return self.is_unmarked() and self.is_increasing() and self.is_contiguous()

    def is_shifted_semistandard(self, require_unmarked_diagonal=False):
        if require_unmarked_diagonal and not self.is_diagonally_unmarked():
            return False
        if not self.is_weakly_increasing():
            return False
        if any((i + 1, j) in self.mapping and self.entry(i, j) == self.entry(i + 1, j) and not self.entry(i, j).is_primed() for (i, j) in self.mapping):
            return False
        if any((i, j + 1) in self.mapping and self.entry(i, j) == self.entry(i, j + 1) and self.entry(i, j).is_primed() for (i, j) in self.mapping):
            return False
        return True

    def is_shifted_standard(self, require_unmarked_diagonal=True):
        return (not require_unmarked_diagonal or self.is_diagonally_unmarked()) and self.is_increasing() and self.is_contiguous()

    @classmethod
    def get_standard(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for i, j in shape.corners():
            for t in cls.get_standard(shape - Shape({(i, j)})):
                mapping = t.mapping
                mapping[(i, j)] = MarkedNumber(len(shape))
                ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_semistandard(cls, shape):
        if type(shape) in [tuple, list]:
            shape = Partition(*shape)
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for border in shape.horizontal_border_strips():
            for t in cls.get_semistandard(shape - Shape(border)):
                if t.mapping:
                    n, mapping = t.maximum().increment(), t.mapping
                else:
                    n, mapping = MarkedNumber(1), {}
                for i, j in border:
                    mapping[(i, j)] = n
                ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_standard_shifted(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for i, j in shape.corners():
            for t in cls.get_standard_shifted(shape - Shape({(i, j)})):
                mapping = t.mapping
                mapping[(i, j)] = MarkedNumber(len(shape))
                ans.add(Tableau(mapping))
                if i != j:
                    mapping[(i, j)] = MarkedNumber(-len(shape))
                    ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_semistandard_shifted(cls, shape, n=None, diagonal_primes=False):
        if type(shape) == tuple:
            shape = StrictPartition(*shape)
        if type(shape) == Partition:
            shape = StrictPartition(shape)
        if type(shape) == StrictPartition:
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        n = len(shape) if n is None else n

        borders = {
            (a, b)
            for a in shape.horizontal_border_strips() | {Shape()}
            for b in (shape - a).vertical_border_strips(exclude_diagonal=not diagonal_primes) | {Shape()}
            if len(a) > 0 or len(b) > 0
        }

        ans = set()
        for border_h, border_v in borders:
            for k in range(n):
                for t in cls.get_semistandard_shifted(shape - border_h - border_v, k, diagonal_primes):
                    mapping = t.mapping if t.mapping else {}
                    for i, j in border_h:
                        mapping[(i, j)] = MarkedNumber(k + 1)
                    for i, j in border_v:
                        mapping[(i, j)] = MarkedNumber(-k - 1)
                    ans.add(Tableau(mapping))
        return ans

    def __le__(self, other):
        assert type(other) == Tableau
        assert set(other.mapping) == set(self.mapping)
        # return all(a <= other.mapping[x] for x, a in self.mapping.items())
        return self.row_reading_word() <= other.row_reading_word()

    def __lt__(self, other):
        assert type(other) == Tableau
        assert set(other.mapping) == set(self.mapping)
        return self.row_reading_word() < other.row_reading_word()

    def __len__(self):
        return len(self.mapping)

    def __nonzero__(self):
        return len(self.mapping) > 0

    def __eq__(self, other):
        assert type(other) == Tableau
        return self.mapping == other.mapping

    def __hash__(self):
        return hash(tuple(sorted(self.mapping.items())))

    def __repr__(self):
        width = max({len(str(v)) for v in self.mapping.values()} | {0})
        base = [['.' + (width - 1) * ' ' for i in range(self.max_column)] for j in range(self.max_row)]
        for i, j in self.mapping:
            v = str(self.mapping[(i, j)])
            base[i - 1][j - 1] = v + (width - len(v)) * ' '
        rows = [' '.join(row) for row in base]
        pad = max([len(a) for a in rows + ['']]) * ' '
        if FRENCH:
            return pad + '\n' + '\n'.join(reversed(rows)) + '\n' + pad   # French
        else:
            return pad + '\n' + '\n'.join(rows) + '\n' + pad             # English

    @classmethod
    def decreasing_part(cls, row):
        if row:
            i = 1
            while i < len(row) and row[i] < row[i - 1]:
                i += 1
            return row[:i]

    @classmethod
    def increasing_part(cls, row):
        if row:
            rev = cls.decreasing_part(list(reversed(row)))
            return tuple(reversed(rev))[1:]

    def replace_row(self, j, newrow, shifted=False):
        dictionary = {(j_, k): self.entry(j_, k) for (j_, k) in self.mapping if j != j_}
        for k_zerobased, v in enumerate(newrow):
            k = k_zerobased + 1
            assert type(v) in [int, MarkedNumber]
            v = v if type(v) == MarkedNumber else MarkedNumber(v)
            if shifted:
                k += j - 1
            dictionary[(j, k)] = v
        return Tableau(dictionary)

    def replace_column(self, j, newcol):
        dictionary = {(i, j_): self.entry(i, j_) for (i, j_) in self.mapping if j != j_}
        for i_zerobased, v in enumerate(newcol):
            i = i_zerobased + 1
            assert type(v) in [int, MarkedNumber]
            v = v if type(v) == MarkedNumber else MarkedNumber(v)
            dictionary[(i, j)] = v
        return Tableau(dictionary)

    def add_to_column(self, j, v):
        return self.replace_column(j, self.get_column(j) + (v, ))

    def add_to_row(self, j, v, shifted=False):
        return self.replace_row(j, self.get_row(j) + (v, ), shifted)

    @classmethod
    def from_shifted_growth_diagram(cls, growth, edges, corners):
        def shdiff(nu, lam):
            return next(iter(Partition.skew(nu, lam, shifted=True)))

        p, q = Tableau(), Tableau()
        n, m = len(growth) - 1, len(growth[0]) - 1 if growth else 0

        for i in range(1, n + 1):
            mu, nu = growth[i][m], growth[i - 1][m]
            for a, b in Partition.skew(mu, nu, shifted=True):
                p = p.add(a, b, i)

        for i in range(1, m + 1):
            mu, nu = growth[n][i], growth[n][i - 1]
            v = -i if edges[n][i] else i
            j = corners[n][i]
            assert mu != nu
            a, b = shdiff(mu, nu)
            q = q.add(a, b, v)

        return p, q

    @classmethod
    def bump(cls, p, column_dir, sequence):
        assert all(sequence[i + 1] > sequence[i] for i in range(len(sequence) - 1))
        if len(sequence) == 0 or p > sequence[-1]:
            newseq = sequence + (p,)
            q = None
        elif p == sequence[-1]:
            newseq = sequence
            q = None
        else:
            if p <= sequence[0]:
                i = 0
                column_dir = True
            else:
                i = [j for j in range(1, len(sequence)) if sequence[j - 1] < p <= sequence[j]][0]
            if p == sequence[i]:
                newseq = sequence
                q = sequence[i + 1]
            else:
                newseq = sequence[:i] + (p,) + sequence[i + 1:]
                q = sequence[i]
        return q, column_dir, newseq

    def remove(self, i, j):
        assert (i, j) in self
        return Tableau({a: b for a, b in self.mapping.items() if a != (i, j)})

    def num_rows(self):
        return len(self.get_rows())

    def get_rows(self):
        ans = []
        for i, j in sorted(self.mapping):
            n = self.mapping[(i, j)].number
            while i > len(ans):
                ans.append([])
            ans[i - 1].append(n)
        return ans

    def get_columns(self):
        ans = []
        for i, j in sorted(self.mapping):
            n = self.mapping[(i, j)].number
            while j > len(ans):
                ans.append([])
            ans[j - 1].append(n)
        return ans

    def row_reading_word(self):
        return tuple(
            self.mapping[key].number for key in sorted(self.mapping, key=lambda x: (-x[0], x[1]))
        )

    def column_reading_word(self):
        return tuple(
            self.mapping[key].number for key in sorted(self.mapping, key=lambda x: (x[1], -x[0]))
        )

    def shifted_reading_word(self):
        return self.shifted_crystal_word()[0]
        # a = tuple(i for i in reversed(self.column_reading_word()) if i < 0)
        # b = tuple(i for i in self.row_reading_word() if i > 0)
        # return a + b

    def descent_set(self):
        w = [abs(a) for a in reversed(self.column_reading_word()) if a < 0]
        w += [abs(a) for a in self.row_reading_word() if a > 0]
        d = {i: pos for pos, i in enumerate(w)}
        return {i for i in d if i + 1 in d and d[i + 1] < d[i]}

    def shifted_descents(self):
        w = self.shifted_reading_word()
        n = len(w)
        assert all(i in w for i in range(1, n + 1))
        d = {i: pos for pos, i in enumerate(w)}
        return {i for i in range(1, n) if d[i + 1] < d[i]}

    @classmethod
    def inverse_sagan_worley(cls, p, q):
        n = len(p)
        if n == 0:
            return ()
        if q.find(n):
            i, j = next(iter(q.find(n).mapping))
            a, p = p.entry(i, j), p.remove(i, j)
            cdir = False
        else:
            i, j = next(iter(q.find(-n).mapping))
            a, p = p.entry(i, j), p.remove(i, j)
            cdir = True
        while i > 1 or cdir:
            if cdir:
                j = j - 1
                i = max([k for (k, l) in p if l == j and p.entry(k, l) <= a])
                a, p = p.entry(i, j), p.set(i, j, a)
                cdir = (i != j)
            else:
                i = i - 1
                j = max([l for (k, l) in p if k == i and p.entry(k, l) < a])
                a, p = p.entry(i, j), p.set(i, j, a)
        return cls.inverse_sagan_worley(p, q) + (a.number,)

    @classmethod
    def inverse_rsk(cls, p, q):
        n = len(p)
        if n == 0:
            return ()
        i, j = next(iter(q.find(n).mapping))
        a = p.entry(i, j)
        p = p.remove(i, j)
        while i > 1:
            i = i - 1
            j = max([l for (k, l) in p if k == i and p.entry(k, l) < a])
            a, p = p.entry(i, j), p.set(i, j, a)
        return cls.inverse_rsk(p, q) + (a.number,)

    def rsk_insert(self, p, j=0):
        p = MarkedNumber(p) if type(p) == int else p
        if p is None:
            return (j, self)

        def rsk_bump(a, tup):
            if len(tup) == 0 or a >= tup[-1]:
                newtup = tup + (a,)
                q = None
            else:
                i = [j for j in range(len(tup)) if a < tup[j]][0]
                newtup = tup[:i] + (a,) + tup[i + 1:]
                q = tup[i]
            return q, newtup

        j += 1
        row = self.get_row(j)
        p, row = rsk_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)
        return tab.rsk_insert(p, j)

    def hecke_insert(self, p, j=0):
        if p is None:
            return (j, self)

        def hecke_bump(a, tup):
            if len(tup) == 0 or a > tup[-1]:
                newtup = tup + (a,)
                q = None
            elif a == tup[-1]:
                newtup = tup
                q = None
            else:
                i = [j for j in range(len(tup)) if a < tup[j]][0]
                if i > 0 and a == tup[i - 1]:
                    newtup = tuple(tup)
                    q = tup[i]
                else:
                    newtup = tup[:i] + (a,) + tup[i + 1:]
                    q = tup[i]
            return q, newtup

        j += 1
        row = self.get_row(j)
        p, row = hecke_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)
        if tab.is_increasing():
            return tab.hecke_insert(p, j)
        else:
            return self.hecke_insert(p, j)

    def shifted_hecke_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = self.bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = self.bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if tab.is_increasing():
            return tab.shifted_hecke_insert(p, j, column_dir, verbose=verbose)
        else:
            return self.shifted_hecke_insert(p, j, column_dir, verbose=verbose)

    def eg_tableau(cls, *x):
        u = Tableau()
        for i in x:
            _, u = u.eg_insert(i)
        return u

    def eg_path(self, *x):
        assert len(x) > 0
        u = self
        for i in x[:-1]:
            _, u = u.eg_insert(i)
        _, v = u.eg_insert(x[-1])
        (new_row, new_col) = next(iter(set(v.mapping) - set(u.mapping)))
        x = v.get(new_row, new_col)
        path = [(new_row, new_col)]
        for i in range(new_row - 1, 0, -1):
            j = 1
            while (i, j) in v.mapping:
                if v.get(i, j) == x or v.get(i, j) != u.get(i, j):
                    path.append((i, j))
                    x = v.get(i, j)
                    break
                j += 1
        return tuple(reversed(path))


    def eg_insert(self, p, j=0):
        if p is None:
            return (j, self)
        elif type(p) is int:
            p = MarkedNumber(p)

        def eg_bump(a, tup):
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if a == b:
                    b = tup[i + 1]
                    new = tup
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, new)
            return (None, tup + (a,))

        j += 1
        row = self.get_row(j)
        p, row = eg_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)

        assert tab.is_increasing()
        return tab.eg_insert(p, j)

    def sp_mixed_insert(self, p, j=1, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def sgn(x):
            return -1 if x.number < 0 else 1

        def bump(a, cdir, tup):
            for i, b in enumerate(tup):
                if abs(a) >= abs(b):
                    continue
                if not cdir and i == 0:
                    new = (MarkedNumber(abs(a) * sgn(b)),) + tup[1:]
                    assert abs(b) % 2 == 0
                    b = MarkedNumber((abs(b) - 1) * sgn(a))
                    cdir = True
                else:
                    new = tup[:i] + (MarkedNumber(abs(a) * sgn(b)),) + tup[i + 1:]
                    b = MarkedNumber(abs(b) * sgn(a))
                    cdir = abs(b) % 2 != 0
                return (b, cdir, new, i)
            return (None, cdir, tup + (a,), None)

        row, col = self.get_row(j), self.get_column(j)

        if column_dir:
            p, column_dir, col, i = bump(p, column_dir, col)
            tab = self.replace_column(j, col)

            if column_dir:
                j += 1
            else:
                j = i + 2
        else:
            p, column_dir, row, i = bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

            if column_dir:
                j = i + j + 1
            else:
                j += 1
        return tab.sp_mixed_insert(p, j, column_dir, verbose=verbose)

    def mixed_insert(self, p, j=1, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def sgn(x):
            return -1 if x.number < 0 else 1

        def bump(a, cdir, tup):
            for i, b in enumerate(tup):
                if abs(a) >= abs(b):
                    continue
                if not cdir and i == 0:
                    new = (a,) + tup[1:]
                    assert abs(b) % 2 == 0
                    b = MarkedNumber((abs(b) - 1) * sgn(b))
                    cdir = True
                else:
                    new = tup[:i] + (MarkedNumber(abs(a) * sgn(b)),) + tup[i + 1:]
                    b = MarkedNumber(abs(b) * sgn(a))
                    cdir = abs(b) % 2 != 0
                return (b, cdir, new, i)
            return (None, cdir, tup + (a,), None)

        row, col = self.get_row(j), self.get_column(j)

        if column_dir:
            p, column_dir, col, i = bump(p, column_dir, col)
            tab = self.replace_column(j, col)

            if column_dir:
                j += 1
            else:
                j = i + 2
        else:
            p, column_dir, row, i = bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

            if column_dir:
                j = i + j + 1
            else:
                j += 1
        return tab.mixed_insert(p, j, column_dir, verbose=verbose)

    def sagan_worley_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def bump(a, cdir, tup):
            for i, b in enumerate(tup):
                eq = (cdir and not a.is_primed()) or (not cdir and a.is_primed())
                if (eq and a > b) or (not eq and a >= b):
                    continue
                if a == b:
                    if not cdir and i == 0:
                        if a.is_primed():
                            b = -a
                        cdir = True
                    new = tup
                elif not cdir and i == 0:
                    if -a == b:
                        a = -a
                    new = (a,) + tup[1:]
                    cdir = True
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))
            # for i, b in enumerate(tup):
            #     if a > b:
            #         continue
            #     if not cdir and a == b:
            #         continue
            #     if not cdir and i == 0:
            #         new = (a,) + tup[1:]
            #         cdir = True
            #     else:
            #         new = tup[:i] + (a,) + tup[i + 1:]
            #     return (b, cdir, new)
            # return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if column_dir:
            p, column_dir, col = bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)
        return tab.sagan_worley_insert(p, j, column_dir, verbose=verbose)

    def count_primes(self):
        ans = 0
        for (i, j) in self.mapping:
            if self.entry(i, j).is_primed():
                ans += 1
        return ans

    def count_diagonal_primes(self):
        ans = 0
        for (i, j) in self.mapping:
            if i == j and self.entry(i, j).is_primed():
                ans += 1
        return ans

    def primed_sw_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def primed_sw_bump(j, a, cdir, tup):
            for i, b in enumerate(tup):
                eq = (cdir and not a.is_primed()) or (not cdir and a.is_primed())
                if (eq and a > b) or (not eq and a >= b):
                    continue
                if a == b:
                    if not cdir and i == 0:
                        assert a.is_primed()
                        b = -a
                        cdir = True
                    new = tup
                elif not cdir and i == 0:
                    assert not b.is_primed()
                    if -a != b:
                        a, b = MarkedNumber.swap_primes(a, b)
                    else:
                        a = -a
                    new = (a,) + tup[1:]
                    cdir = True
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = primed_sw_bump(j, p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = primed_sw_bump(j, p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if verbose:
            print(tab, '\n')
        assert tab.is_shifted_semistandard(False)
        return tab.primed_sw_insert(p, j, column_dir, verbose=verbose)

    def involution_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def involution_bump(j, a, cdir, tup):
            for i, b in enumerate(tup):
                if abs(a) > abs(b):
                    continue
                # print('j =', j, 'a =', a, 'cdir =', cdir, 'tup =', tup, 'i =', i + 1, 'b =', b)
                if abs(a) == abs(b):
                    assert not a.is_primed() or not b.is_primed()
                    b = a.increment()
                    new = tup[:i] + MarkedNumber.swap_primes(tup[i], tup[i + 1]) + tup[i + 2:]
                    cdir = cdir or (i == 0)
                elif not cdir and i == 0:
                    a, b = MarkedNumber.swap_primes(a, b)
                    new = (a,) + tup[1:]
                    cdir = True
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = involution_bump(j, p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = involution_bump(j, p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if verbose:
            print(tab, '\n')
        assert tab.is_increasing()
        return tab.involution_insert(p, j, column_dir, verbose=verbose)

    def fpf_insert(self, p, j=0, column_dir=False, verbose=False):
        if p is None:
            return (j, column_dir, self)

        def fpf_bump(a, cdir, tup):
            # inserting `a` (number) into `tup` (tuple) in row (`cdir=False`) or column direction
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if a == b:
                    b = tup[i + 1]
                    assert b == a.increment()
                    new = tup
                elif not cdir and i == 0:
                    cdir = True
                    if a.number % 2 == 0:
                        new = (a,) + tup[1:]
                    else:
                        b = a.increment().increment()
                        new = tup
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = fpf_bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = fpf_bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        assert tab.is_increasing()
        return tab.fpf_insert(p, j, column_dir, verbose=verbose)

    @cached_value(HORIZONTAL_STRIPS_CACHE)
    def _horizontal_strips(cls, mu, lam):  # noqa
        if not Partition._contains(mu, lam):
            return []

        core = [mu[i + 1] if i + 1 < len(mu) else 0 for i in range(len(mu))]
        for i in range(len(lam)):
            core[i] = max(core[i], lam[i])
        core = tuple(core)

        ans = []
        level = {core}
        while level:
            for nu in level:
                diff = {(i + 1, j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                nu = nu if nu and nu[-1] > 0 else nu[:-1]
                corners = [(i + 1, nu[i]) for i in range(len(nu)) if core[i] < nu[i]]
                ans.append((nu, diff, corners))
            level = {
                nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                for i in range(len(mu))
                for nu in level
                if nu[i] < mu[i]
            }
        return ans

    @classmethod
    def _subsets(cls, diff, corners, setvalued):
        if setvalued:
            for v in range(2**len(corners)):
                thisdiff = diff
                for i in range(len(corners)):
                    thisdiff = thisdiff if v % 2 == 0 else thisdiff | {corners[i]}
                    v = v // 2
                yield thisdiff
        else:
            yield diff

    @classmethod
    def shifted_k_flagged(cls, k, mu):  # noqa
        max_entry = max(len(mu), mu[0] if mu else 0) + k
        return {t for t in cls.get_semistandard_shifted(mu, max_entry) if t.is_shifted_k_flagged(k)}

    @classmethod
    def k_flagged(cls, k, mu):  # noqa
        max_entry = len(mu) + k
        return {t for t in cls.semistandard(max_entry, mu) if t.is_k_flagged(k)}

    @classmethod
    def semistandard(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._semistandard(max_entry, mu, nu, setvalued)

    @cached_value(SEMISTANDARD_CACHE)
    def _semistandard(cls, max_entry, mu, lam, setvalued):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam) and max_entry > 0:
            for nu, diff, corners in cls._horizontal_strips(mu, lam):
                for aug in cls._subsets(diff, corners, setvalued):
                    for tab in cls._semistandard(max_entry - 1, nu, lam, setvalued):
                        for (i, j) in aug:
                            tab = tab.add(i, j, max_entry)
                        ans.add(tab)
        return ans

    @classmethod
    def standard(cls, mu, nu=()):  # noqa
        return cls._standard(mu, nu)

    @cached_value(STANDARD_CACHE)
    def _standard(cls, mu, lam):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), mu[i]
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_partition(nu):
                    for tab in cls._standard(nu, lam):
                        ans.add(tab.add(row, col, n))
        return ans

    @classmethod
    def standard_shifted_marked(cls, mu, nu=(), diagonal_primes=False):  # noqa
        return cls._standard_shifted_marked(mu, nu, diagonal_primes)

    @cached_value(STANDARD_SHIFTED_MARKED_CACHE)
    def _standard_shifted_marked(cls, mu, lam, diagonal_primes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), (i + mu[i])
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_strict_partition(nu):
                    for tab in cls._standard_shifted_marked(nu, lam, diagonal_primes):
                        ans.add(tab.add(row, col, n))
                        if diagonal_primes or row != col:
                            ans.add(tab.add(row, col, -n))
        return ans

    @classmethod
    def even_diagonal_unprimed_shifted_rpp(cls, max_entry, mu, nu=()):  # noqa
        return cls.unprimed_shifted_rpp(max_entry, mu, nu, True)

    @classmethod
    def unprimed_shifted_rpp(cls, max_entry, mu, nu=(), even_diagonal=False):  # noqa
        ans = {t.unprime() for t in cls.semistandard_marked_rpp(max_entry, mu, nu)}
        if even_diagonal:
            ans = {t for t in ans if all(t.entry(i, i).number % 2 != 0 for i in range(1, t.max_row + 1))}
        return ans

    @classmethod
    def semistandard_marked_rpp(cls, max_entry, mu, nu=(), diagonal_nonprimes=False):  # noqa
        return cls._semistandard_marked_rpp(max_entry, mu, nu, diagonal_nonprimes)

    @cached_value(SEMISTANDARD_MARKED_RPP_CACHE)
    def _semistandard_marked_rpp(cls, max_entry, mu, lam, diagonal_nonprimes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam) and max_entry > 0:
            for nu1, diff1 in cls._shifted_rpp_horizontal_strips(mu):
                for nu2, diff2 in cls._shifted_rpp_vertical_strips(nu1):
                    if diagonal_nonprimes or not any(i == j for i, j in diff1):
                        for tab in cls._semistandard_marked_rpp(max_entry - 1, nu2, lam, diagonal_nonprimes):
                            for (i, j) in diff1:
                                tab = tab.add(i, j, max_entry)
                            for (i, j) in diff2:
                                tab = tab.add(i, j, -max_entry)
                            ans.add(tab)
        return ans

    @cached_value(SHIFTED_RPP_HORIZONTAL_STRIPS_CACHE)
    def _shifted_rpp_horizontal_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        if mu == ():
            return [(mu, set())]

        def remove_box(nu, i):
            if i < len(nu) and nu[i] > 0:
                nu = nu[:i] + (nu[i] - 1,) + nu[i + 1:]
                while nu and nu[-1] == 0:
                    nu = nu[:-1]
                if all(nu[j] > nu[j + 1] for j in range(len(nu) - 1)):
                    yield nu

        def remove_all_boxes(nu, i):
            queue = [nu]
            while queue:
                nu, queue = queue[0], queue[1:]
                yield nu
                for x in remove_box(nu, i):
                    queue.append(x)

        def skew(mu, nu):
            ans = set()
            for i, part in enumerate(mu):
                subpart = nu[i] if i < len(nu) else 0
                for j in range(subpart, part):
                    ans.add((i + 1, j + 1 + i))
            return ans

        ans = set()
        queue = [(mu, len(mu) - 1)]
        while queue:
            nu, i = queue[0]
            queue = queue[1:]
            if i >= 0:
                for nu in remove_all_boxes(nu, i):
                    ans.add(nu)
                    queue.append((nu, i - 1))

        return [(nu, skew(mu, nu)) for nu in ans]

    @cached_value(SHIFTED_RPP_VERTICAL_STRIPS_CACHE)
    def _shifted_rpp_vertical_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        if mu == ():
            return [(mu, set())]

        def remove_box(nu, i):
            for j in range(len(nu) - 1, -1, -1):
                if j + nu[j] == i + 1:
                    nu = nu[:j] + (nu[j] - 1,) + nu[j + 1:]
                    while nu and nu[-1] == 0:
                        nu = nu[:-1]
                    yield nu
                    return

        def remove_all_boxes(nu, i):
            queue = [nu]
            while queue:
                nu, queue = queue[0], queue[1:]
                yield nu
                for x in remove_box(nu, i):
                    queue.append(x)

        def skew(mu, nu):
            ans = set()
            for i, part in enumerate(mu):
                subpart = nu[i] if i < len(nu) else 0
                for j in range(subpart, part):
                    ans.add((i + 1, j + 1 + i))
            return ans

        ans = set()
        queue = [(mu, (mu[0] if mu else 0) - 1)]
        while queue:
            nu, i = queue[0]
            queue = queue[1:]
            if i >= 0:
                for nu in remove_all_boxes(nu, i):
                    ans.add(nu)
                    queue.append((nu, i - 1))

        return [(nu, skew(mu, nu)) for nu in ans]

    @classmethod
    def rpp(cls, mu, k):
        ans = set()
        for t in cls.k_flagged(k, mu):
            ans.add(Tableau({b: v.number - b[0] for b, v in t.mapping.items()}))
        return ans

    @classmethod
    def shrpp(cls, mu, k, require_even_diag=True):
        ans = set()
        for t in {t.unprime() for t in Tableau.semistandard_marked_rpp(k + 1, mu)}:
            has_all_even_diag = all(t.entry(i, i).number % 2 != 0 for i in range(1, t.max_row + 1))
            if not require_even_diag or has_all_even_diag:
                ans.add(Tableau({b: v.number - 1 for b, v in t.mapping.items()}))
        return ans

    def rpp_weight(self):
        ans = []
        bns = []
        for i, j in self:
            x = self[(i, j)].number
            assert x != 0
            if x > 0:
                while x - 1 >= len(ans):
                    ans.append(set())
                ans[x - 1].add(j)
            if x < 0:
                while -x - 1 >= len(bns):
                    bns.append(set())
                bns[-x - 1].add(i)
        n = max(len(ans), len(bns))
        while n > len(ans):
            ans.append(set())
        while n > len(bns):
            bns.append(set())
        return tuple(len(ans[i]) + len(bns[i]) for i in range(n))
