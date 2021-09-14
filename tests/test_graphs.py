from vectors import Vector
import random
import itertools
import subprocess
import time
import numpy


def write_graph(comp, part, write=True):
    if type(comp[0]) == int:
        n = sum(comp)
        des = [0]
        for c in comp:
            des += [des[-1] + c]
        edges = [(i, i + 1, 1) for i in range(1, sum(comp)) if i not in des]
    else:
        n = sum(map(len, comp))
        edges = [(istring[i], istring[i + 1], 1) for istring in comp for i in range(len(istring) - 1)]
    edges += [(istring[i], istring[i + 1], 2) for istring in part for i in range(len(istring) - 1)]

    if len({a for (a, _, _) in edges} | {b for (_, b, _) in edges}) < n:
        return False

    if not write:
        return True

    directory = "/Users/emarberg/examples/crystals/random/" + ("edges_%s_" % len(edges)) + str(time.time_ns())
    dot_filename = directory + ".dot"
    png_filename = directory + ".png"

    s = ["digraph G {"]
    for a, b, i in edges:
        s += ["    %s -> %s [ label=\"%s\" ];" % (a, b, i)]
    s += ["}"]
    s = "\n".join(s)

    with open(dot_filename, "w") as f:
        f.write(s)
    subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
    subprocess.run(["open", png_filename])
    return True


def all_comps(n, mpart=None):
    if n == 0:
        yield ()
        return
    for a in range(1, 1 + min(n, (mpart or n))):
        for c in all_comps(n - a, a):
            yield (a,) + c


def all_partitions(base, mpart=None):
    base = set(range(1, base + 1)) if type(base) == int else base
    n = len(base)
    if n == 0:
        yield ()
        return
    for a in range(1, 1 + min(n, (mpart or n))):
        for sub in itertools.combinations(base, a):
            for p in all_partitions(set(base) - set(sub), a):
                for perm in itertools.permutations(sub):
                    ans = (perm,) + p
                    if all(len(ans[i]) > len(ans[i + 1]) or min(ans[i]) < min(ans[i + 1]) for i in range(len(ans) - 1)):
                        yield ans


def rcomp(n):
    if n == 0:
        return ()
    a = random.randint(1, n)
    return (a,) + rcomp(n - a)


def rpartition(base):
    base = set(range(1, base + 1)) if type(base) == int else base
    n = len(base)
    if n == 0:
        return ()
    a = random.randint(1, n)
    subsets = list(itertools.combinations(base, a))
    sub = subsets[-1 + random.randint(1, len(subsets))]
    permutations = list(itertools.permutations(sub))
    perm = permutations[-1 + random.randint(1, len(permutations))]
    return (perm,) + rpartition(set(base) - set(sub))


def get_system(comp, part):
    if type(comp[0]) == int:
        n = sum(comp)
        comp = list(comp)
        for i in range(len(comp) - 1, -1, -1):
            comp[i] = tuple(range(sum(comp[:i]), sum(comp[:i + 1])))
    else:
        n = sum(map(len, comp))
        comp = [tuple(i - 1 for i in e) for e in comp]
    part = [tuple(i - 1 for i in e) for e in part]
    mat = []
    for istring in comp:
        i = istring[0]
        a, b = 3 * i, 3 * i + 1
        # x_a - x_b == len(istring) - 1
        row = [0 for _ in range(3 * n + 1)]
        row[a], row[b], row[-1] = 1, -1, len(istring) - 1
        mat.append(row)

        for index in range(1, len(istring)):
            i, j = istring[index - 1], istring[index]
            ai, bi, ci = 3 * i, 3 * i + 1, 3 * i + 2
            aj, bj, cj = 3 * j, 3 * j + 1, 3 * j + 2
            # x_ai - x_aj == 1
            row = [0 for _ in range(3 * n + 1)]
            row[ai], row[aj], row[-1] = 1, -1, 1
            mat.append(row)
            # x_bj - x_bi == 1
            row = [0 for _ in range(3 * n + 1)]
            row[bj], row[bi], row[-1] = 1, -1, 1
            mat.append(row)
            # x_ci == x_cj
            row = [0 for _ in range(3 * n + 1)]
            row[ci], row[cj], row[-1] = 1, -1, 0
            mat.append(row)

    for istring in part:
        i = istring[0]
        a, b = 3 * i + 1, 3 * i + 2
        # x_a - x_b == len(istring) - 1
        row = [0 for _ in range(3 * n + 1)]
        row[a], row[b], row[-1] = 1, -1, len(istring) - 1
        mat.append(row)

        for index in range(1, len(istring)):
            i, j = istring[index - 1], istring[index]
            ai, bi, ci = 3 * i + 1, 3 * i + 2, 3 * i
            aj, bj, cj = 3 * j + 1, 3 * j + 2, 3 * j
            # x_ai - x_aj == 1
            row = [0 for _ in range(3 * n + 1)]
            row[ai], row[aj], row[-1] = 1, -1, 1
            mat.append(row)
            # x_bj - x_bi == 1
            row = [0 for _ in range(3 * n + 1)]
            row[bj], row[bi], row[-1] = 1, -1, 1
            mat.append(row)
            # x_ci == x_cj
            row = [0 for _ in range(3 * n + 1)]
            row[ci], row[cj], row[-1] = 1, -1, 0
            mat.append(row)

    return mat


def test(n=4):
    i, j, k = 0, 0, 0
    pairs = [(comp, part) for comp in all_comps(n) for part in all_partitions(n)]
    q = len(pairs)
    for index in numpy.random.permutation(q):
        comp, part = pairs[index]
        i += 1
        mat = get_system(comp, part)
        if Vector.is_consistent_linear_system(mat):
            j += 1
            if write_graph(comp, part):
                k += 1
                print(comp, part)
                print()
                for row in mat:
                    print('  ', row)
                print()
                print('weyl group action:', check_weyl_action(comp, part))
                # rref = Vector.rref(mat)
                # print()
                # for row in rref:
                #     print('  ', row)
                # print()
                # input('')
        if i % 1000 == 0:
            print('considered', i, 'of', q, 'consistent', j, 'crystals', k)


def istring_reflect(x, strings):
    for s in strings:
        if x in s:
            i = [_ for _ in range(len(s)) if s[_] == x][0]
            return s[len(s) - 1 - i]


def check_weyl_action(onestrings, twostrings):
    if type(onestrings[0]) == int:
        onestrings = list(onestrings)
        for i in range(len(onestrings) - 1, -1, -1):
            onestrings[i] = tuple(range(sum(onestrings[:i]), sum(onestrings[:i + 1])))
        twostrings = [tuple(i - 1 for i in e) for e in twostrings]

    elems = {e for _ in onestrings for e in _} | {e for _ in twostrings for e in _}
    for e in elems:
        a = istring_reflect(istring_reflect(istring_reflect(e, onestrings), twostrings), onestrings)
        b = istring_reflect(istring_reflect(istring_reflect(e, twostrings), onestrings), twostrings)
        if a != b:
            return False
    return True


def test_abnormal():
    onestrings = [[1, 18, 21], [2, 22], [4, 10, 19], [6, 0, 11], [9], [12, 20], [17, 16, 24, 5], [25, 8, 13, 3], [26, 15, 14, 7, 23]]
    twostrings = [[1, 25, 26], [3, 19, 22, 0], [7, 24], [8, 15], [13, 10, 16], [14], [18, 12, 4, 17], [21, 20, 9, 2, 6], [23, 5, 11]]
    m, i = {}, 1
    for c in onestrings:
        for v in c:
            m[v], i = i, i + 1
    comp = list(map(len, onestrings))
    part = [[m[v] for v in p] for p in twostrings]
    if Vector.is_consistent_linear_system(get_system(comp, part)):
        if write_graph(comp, part, False):
            print('weyl group action (control):', check_weyl_action(onestrings, twostrings))

    # altstrings = [[1, 25, 26], [3, 19, 24, 0], [7, 22], [8, 15], [13, 10, 16], [14], [18, 12, 4, 17], [21, 20, 9, 2, 6], [23, 5, 11]]
    # part = [[m[v] for v in p] for p in altstrings]
    onestrings = [[1, 2, 4], [17, 22], [8, 14, 18], [21, 25, 27], [13], [6, 9], [12, 19, 23, 26], [3, 7, 10, 15], [5, 11, 16, 20, 24]]
    altstrings = [[1, 3, 5], [15, 18, 23, 25], [20, 22], [7, 11], [10, 14, 19], [16], [2, 6, 8, 12], [4, 9, 13, 17, 21], [24, 26, 27]]
    if Vector.is_consistent_linear_system(get_system(onestrings, altstrings)):
        if write_graph(onestrings, altstrings, False):
            print('weyl group action (experiment):', check_weyl_action(onestrings, altstrings))
            print(sorted(onestrings, key=lambda x: x[0]))
            print(sorted(altstrings, key=lambda x: x[0]))
    # heights = {1: 0}
    # while len(heights) < 27:
    #     for s in onestrings:
    #         for i in range(1, len(s)):
    #             if s[i] not in heights and s[i - 1] in heights:
    #                 heights[s[i]] = heights[s[i - 1]] + 1
    #     for s in altstrings:
    #         for i in range(1, len(s)):
    #             if s[i] not in heights and s[i - 1] in heights:
    #                 heights[s[i]] = heights[s[i - 1]] + 1
    # sort = sorted(heights, key=lambda x: heights[x])
    # m, i = {}, 1
    # for a in sort:
    #     m[a], i = i, i + 1
    # onestrings = [[m[a] for a in r] for r in onestrings]
    # altstrings = [[m[a] for a in r] for r in altstrings]
    # print(onestrings)
    # print(altstrings)
    # print()
