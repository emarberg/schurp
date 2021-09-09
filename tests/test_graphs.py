from vectors import Vector
import random
import itertools
import subprocess
import time


def write_graph(comp, part):
    des = [0]
    for c in comp:
        des += [des[-1] + c]
    edges = [(i, i + 1, 1) for i in range(1, sum(comp)) if i not in des]
    edges += [(istring[i], istring[i + 1], 2) for istring in part for i in range(len(istring) - 1)]

    if len({a for (a, _, _) in edges} | {b for (_, b, _) in edges}) < sum(comp):
        return False

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
    n = sum(comp)
    comp = list(comp)
    for i in range(len(comp) - 1, -1, -1):
        comp[i] = tuple(range(sum(comp[:i]), sum(comp[:i + 1])))
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


def test(n):
    for comp in all_comps(n):
        for part in all_partitions(n):
            mat = get_system(comp, part)
            if Vector.is_consistent_linear_system(mat):
                if write_graph(comp, part):
                    print(comp, part)
                    print()
                    for row in mat:
                        print('  ', row)
                    print()
                    # rref = Vector.rref(mat)
                    # print()
                    # for row in rref:
                    #     print('  ', row)
                    # print()
                    input('')
