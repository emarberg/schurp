from stable.vectors import Vector
from stable.partitions import Partition
import pytest


def compose(*args):
    def fn(x):
        vec = x if type(x) == Vector else Vector({x: 1})
        for op in reversed(args):
            ans = Vector(printer=vec.printer)
            for mu in vec:
                ans += vec[mu] * op(mu)
            vec = ans
        return vec
    return fn


def u(i):
    def op(mu):
        for j in range(len(mu)):
            if mu[j] == i - 1:
                nu = mu[:j] + (i,) + mu[j + 1:]
                return Vector({nu: 1})
        if i == 1:
            nu = mu + (1,)
            return Vector({nu: 1})
        return Vector()
    return op


def d(i):
    def op(mu):
        for j in reversed(range(len(mu))):
            if mu[j] == i:
                nu = (mu[:j] + (i - 1,) + mu[j + 1:]) if i > 1 else mu[:j]
                return Vector({nu: 1})
        return Vector()
    return op


def skew_symmetric_double(mu):
    shape = {(a + 1, a + b + 1) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
    shape |= {(b, a) for (a, b) in shape}
    a = 1
    while True:
        shape |= {(a, a)}
        if (a, a + 1) not in shape:
            break
        a += 1
    return shape


def skew_symmetric_undouble(shape):
    ans = []
    for a, b in shape:
        if a < b:
            while a > len(ans):
                ans += [0]
            ans[a - 1] += 1
    return tuple(ans)


def symmetric_double(mu):
    shape = {(a + 1, a + b) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
    shape |= {(b, a) for (a, b) in shape}
    return shape


def symmetric_undouble(shape):
    ans = []
    for a, b in shape:
        if a <= b:
            while a > len(ans):
                ans += [0]
            ans[a - 1] += 1
    return tuple(ans)


# def u_skew_symmetric(i):
#     i, colummwise = abs(i), i > 0
#     def op(mu):
#         shape = skew_symmetric_double(mu)
#         columns = {b for a, b in shape if a == i}
#         b = (max(columns) + 1) if columns else 1
#         if columnwise and i < b:
#             return Vector()
#         if i == 1 or (i - 1, b) in shape:
#             shape |= {(i, b), (b, i)}
#             return Vector({skew_symmetric_undouble(shape): 1})
#         else:
#             return Vector()
#     return op


# def d_skew_symmetric(i):
#     i, colummwise = abs(i), i > 0
#     def op(mu):
#         shape = skew_symmetric_double(mu)
#         columns = {b for a, b in shape if a == i and b != i}
#         b = max(columns) if columns else 0
#         if colummwise and i < b:
#             return Vector()
#         if b > 0 and (i + 1, b) not in shape:
#             shape -= {(i, b), (b, i)}
#             return Vector({skew_symmetric_undouble(shape): 1})
#         else:
#             return Vector()
#     return op


# def u_symmetric(i):
#     i, colummwise = abs(i), i > 0
#     def op(mu):
#         shape = symmetric_double(mu)
#         columns = {b for a, b in shape if a == i}
#         b = (max(columns) + 1) if columns else 1
#         if (colummwise and i < b) or (not colummwise and i > b):
#             return Vector()
#         if i == 1 or (i - 1, b) in shape:
#             shape |= {(i, b), (b, i)}
#             return Vector({symmetric_undouble(shape): 1})
#         else:
#             return Vector()
#     return op


# def d_symmetric(i):
#     i, colummwise = abs(i), i > 0
#     def op(mu):
#         shape = symmetric_double(mu)
#         columns = {b for a, b in shape if a == i}
#         b = max(columns) if columns else 0
#         if (colummwise and i < b) or (not colummwise and i > b):
#             return Vector()
#         if b > 0 and (i + 1, b) not in shape:
#             shape -= {(i, b), (b, i)}
#             return Vector({symmetric_undouble(shape): 1})
#         else:
#             return Vector()
#     return op


# def u_shifted(i, restricted_diagonal=False):
#     def op(mu):
#         if i == 0:
#             return Vector({mu: 1})

#         if i < 0:
#             i = -i - 1
#             if len(mu) < i or (0 < i < len(mu) and mu[i - 1] == mu[i] + 1) or (len(mu) == i and mu[-1] == 1):
#                 return Vector()
#             elif restricted_diagonal and i == len(mu):
#                 return Vector()
#             else:
#                 mu = mu + (0,)
#                 nu = mu[:i] + (mu[i] + 1,) + mu[i + 1:]
#                 nu = nu[:-1] if nu[-1] == 0 else nu
#                 return Vector({nu: 1})

#         if i > 0:
#             shape = {(a + 1, a + b) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
#             rows = {a for a, b in shape if i == b}
#             a = (max(rows) + 1) if rows else 0
#             if a == 0 or (a != i and (a, i - 1) not in shape):
#                 return Vector()
#             else:
#                 mu = mu + (0,)
#                 nu = mu[:a - 1] + (mu[a - 1] + 1,) + mu[a:]
#                 nu = nu[:-1] if nu[-1] == 0 else nu
#                 return Vector({nu: 1})
#     return op


# def d_shifted(i, restricted_diagonal=False):
#     def op(mu):
#         if i == 0:
#             return Vector({mu: 1})

#         if i < 0:
#             i = -i - 1
#             if len(mu) <= i or (i + 1 < len(mu) and mu[i + 1] == mu[i] - 1):
#                 return Vector()
#             else:
#                 nu = mu[:i] + (mu[i] - 1,) + mu[i + 1:]
#                 nu = nu[:-1] if nu[-1] == 0 else nu
#                 return Vector({nu: 1})

#         if i > 0:
#             shape = {(a + 1, a + b) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
#             rows = {a for a, b in shape if i == b}
#             a = (max(rows) - 1) if rows else -1
#             if a == -1 or mu[a] + a != i:
#                 return Vector()
#             elif restricted_diagonal and a + 1 == len(mu) and mu[-1] == 1:
#                 return Vector()
#             else:
#                 nu = mu[:a] + (mu[a] - 1,) + mu[i + 1:]
#                 nu = nu[:-1] if nu[-1] == 0 else nu
#                 return Vector({nu: 1})
#     return op


def u_symmetric(i):
    def op(mu):
        shape = {(a + 1, a + b) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
        a, b = 1, 1 + abs(i)
        while (a, b) in shape:
            a += 1
            b += 1
        if a > 1 and (a - 1, b) not in shape:
            return Vector()
        if a < b and (a, b - 1) not in shape:
            return Vector()
        shape |= {(a, b)}
        return Vector({symmetric_undouble(shape): 1})
    return op


def d_symmetric(i):
    def op(mu):
        shape = {(a + 1, a + b) for a in range(len(mu)) for b in range(1, mu[a] + 1)}
        cells = {(a, b) for (a, b) in shape if abs(a - b) == abs(i)}
        if len(cells) == 0:
            return Vector()
        a, b = max(cells)
        if (a + 1, b) in shape or (a, b + 1) in shape:
            return Vector()
        shape -= {(a, b)}
        return Vector({symmetric_undouble(shape): 1})
    return op


def small_schur_expressions(n):
    yield (), compose()
    for i in range(1, n + 1):
        yield (i,), compose(u(i))
        yield (-i,), compose(d(i))

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            yield (i, j), compose(u(i), u(j))
            yield (-i, -j), compose(d(i), d(j))
            yield (i, -j), compose(u(i), d(j))
            yield (-i, j), compose(d(i), u(j))

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            for k in range(1, n + 1):
                yield (i, j, k), compose(u(i), u(j), u(k))
                yield (-i, -j, -k), compose(d(i), d(j), d(k))


def calc_span(relations, tup):
    ans = set()
    add = {tup}
    while add:
        ans |= add
        new_add = set()

        for tup in add:
            for m in range(1, len(tup) + 1):
                for i in range(len(tup) - m + 1):
                    a, b, c = tup[:i], tup[i:i + m], tup[i + m:]
                    x = relations.get(b, b)
                    if x != b and len(x) <= len(b):
                        new_add.add(a + x + c)

        add = new_add - ans
    return ans


def get_dict(expected):
    dictionary = {}
    for pair in expected:
        a, b = tuple(pair)
        dictionary[a] = b
        dictionary[b] = a
    return dictionary


ZERO = (('z', 0),)


def get_span(dictionary, keys):
    span = {k: calc_span(dictionary, k) for k in keys}
    for k in span:
        if any(ZERO[0] in t for t in span[k]):
            span[k] = {ZERO}
    return span


@pytest.mark.slow
def test_schur_relations(n=4):

    def tag(args):
        if len(args) == 0:
            return 'id'
        s = ''
        for a in args:
            if a > 0:
                s += 'u_%i ' % a
            else:
                s += 'd_%i ' % -a
        return s.strip()

    square = (n + 1) * (n + 1,)
    partitions = list(Partition.subpartitions(square))
    small = {k: v for k, v in small_schur_expressions(n)}

    def equals(k1, k2):
        op1, op2 = small[k1], small[k2]
        return all(op1(mu) == op2(mu) for mu in partitions)

    expected = [
        {(), (-1, 1)},
    ] + [
        {(i, j), (j, i)} for i in range(1, n + 1) for j in range(1, n + 1) if i + 1 < j
    ] + [
         {(-i, -j), (-j, -i)} for i in range(1, n + 1) for j in range(1, n + 1) if i + 1 < j
    ] + [
        {(-i, j), (j, -i)} for i in range(1, n + 1) for j in range(1, n + 1) if i != j
    ] + [
        {(-i - 1, i + 1), (i, -i)} for i in range(1, n)
    ] + [
        {(i + 1, i, i), (i, i + 1, i)} for i in range(1, n)
    ] + [
        {(i + 1, i, i + 1), (i + 1, i + 1, i)} for i in range(1, n)
    ] + [
        {(-i, -i - 1, -i), (-i, -i, -i - 1)} for i in range(1, n)
    ] + [
        {(-i, -i - 1, -i - 1), (-i - 1, -i, -i - 1)} for i in range(1, n)
    ]
    dictionary = get_dict(expected)
    span = get_span(dictionary, small)

    def check_expected_identities(k1, k2, do_check=True):
        if {k1, k2} in expected:
            if do_check:
                try:
                    assert equals(k1, k2)
                except:
                    print()
                    print('** failure: ', tag(k1), '!=', tag(k2))
                    raise Exception
            return True
        return False

    relations = set()
    for k1 in small:
        for k2 in small:
            if k1 < k2:
                if check_expected_identities(k1, k2):
                    relations.add((k1, k2))
                    print()
                    print('* expected: ', tag(k1), '=', tag(k2))
                    continue
                if span[k1] & span[k2]:
                    continue
                if equals(k1, k2):
                    relations.add((k1, k2))
                    print()
                    print('unexpected:', tag(k1), '=', tag(k2))
                    raise Exception
    print()
    return relations


def small_q_schur_expressions(n):
    f, g = u_symmetric, d_symmetric

    yield ZERO, lambda x: Vector()
    yield (), compose()

    for i in range(0, n + 1):
        yield (('u', i),), compose(f(i))
        yield (('d', i),), compose(g(i))

    for i in range(0, n + 1):
        for j in range(0, n + 1):
            yield (('u', i), ('u', j)), compose(f(i), f(j))
            yield (('d', i), ('d', j)), compose(g(i), g(j))
            yield (('u', i), ('d', j)), compose(f(i), g(j))
            yield (('d', i), ('u', j)), compose(g(i), f(j))

    for i in range(0, n + 1):
        for j in range(0, n + 1):
            for k in range(0, n + 1):
                yield (('u', i), ('u', j), ('u', k)), compose(f(i), f(j), f(k))
                yield (('d', i), ('d', j), ('d', k)), compose(g(i), g(j), g(k))

    for i in range(0, n):
        yield (('u', i), ('u', i + 1), ('u', i), ('u', i + 1)), compose(f(i), f(i + 1), f(i), f(i + 1))
        yield (('u', i + 1), ('u', i), ('u', i + 1), ('u', i)), compose(f(i + 1), f(i), f(i + 1), f(i))
        yield (('d', i), ('d', i + 1), ('d', i), ('d', i + 1)), compose(g(i), g(i + 1), g(i), g(i + 1))
        yield (('d', i + 1), ('d', i), ('d', i + 1), ('d', i)), compose(g(i + 1), g(i), g(i + 1), g(i))


@pytest.mark.slow
def test_shifted_schur_relations(n=4):

    def tag(args):
        if len(args) == 0:
            return 'id'
        if args == ZERO:
            return '0'
        s = ''
        for a in args:
            s += '%s_{%i} ' % a
        return s.strip()

    delta = tuple(range(n + 1, 0, -1))
    partitions = list(Partition.subpartitions(delta, strict=True))
    small = {k: v for k, v in small_q_schur_expressions(n)}

    def equals(k1, k2):
        op1, op2 = small[k1], small[k2]
        return all(op1(mu) == op2(mu) for mu in partitions)

    expected = [
        {(('u', i), ('u', i + 1), ('u', i)), ZERO} for i in range(1, n)
    ] + [
        {(('d', i), ('d', i + 1), ('d', i)), ZERO} for i in range(1, n)
    ] + [
        {(('u', i), ('u', i)), ZERO} for i in range(0, n + 1)
    ] + [
        {(('u', i + 1), ('u', i), ('u', i + 1)), ZERO} for i in range(0, n)
    ] + [
        {(('d', i), ('d', i)), ZERO} for i in range(0, n + 1)
    ]+ [
        {(('d', i + 1), ('d', i), ('d', i + 1)), ZERO} for i in range(0, n)
    ] + [
        {(('u', i + 1), ('d', i)), ZERO} for i in range(0, n)
    ] + [
        {(('d', i), ('u', i + 1)), ZERO} for i in range(0, n)
    ] + [
        {(('d', i + 1), ('u', i)), ZERO} for i in range(0, n)
    ] + [
        {(('u', i), ('d', i + 1)), ZERO} for i in range(0, n)
    ] + [
        {(('u', i), ('u', j)), (('u', j), ('u', i))} for i in range(0, n + 1) for j in range(0, n + 1) if i + 1 < j
    ] + [
        {(('d', i), ('d', j)), (('d', j), ('d', i))} for i in range(0, n + 1) for j in range(0, n + 1) if i + 1 < j
    ] + [
        {(('u', i), ('d', j)), (('d', j), ('u', i))} for i in range(0, n + 1) for j in range(0, n + 1) if abs(i - j) > 1
    ]
    span = get_span(get_dict(expected), small)
    print('\n.\n')

    def check_expected_identities(k1, k2, do_check=True):
        if {k1, k2} in expected:
            if do_check:
                try:
                    assert equals(k1, k2)
                except:
                    print()
                    print('** failure: ', tag(k1), '!=', tag(k2))
                    # for mu in partitions:
                    #     op1 = small[k1]
                    #     op2 = small[k2]
                    #     if op1(mu) != op2(mu):
                    #         print('  ', mu, ':', op1(mu), '!=', op2(mu))
                    # print()
                    raise Exception
            return True
        return False

    relations = set()
    for k1 in sorted(small):
        for k2 in sorted(small):
            if k1 < k2:
                if check_expected_identities(k1, k2):
                    relations.add((k1, k2))
                    # print()
                    # print('* expected: ', tag(k1), '=', tag(k2))
                    continue
                if span[k1] & span[k2]:
                    continue

#                if all(small[k1](mu).is_zero() for mu in partitions) and k2 != ZERO:
#                     continue

                if equals(k1, k2):
                    relations.add((k1, k2))
                    print()
                    print('unexpected:', tag(k1), '=', tag(k2))
                    raise Exception
    print()
    return relations
