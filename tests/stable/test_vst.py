from stable.tableaux import Tableau
from stable.vst import combine_str, ValuedSetTableau
from stable.partitions import Partition
from stable.polynomials import Polynomial
from stable.utils import beta, jp_expansion, jp, jq_expansion, jq
import traceback


def shifted_ribbon_q_params(kappa, nu):
    boxes = Partition.shifted_shape(kappa, nu)
    boxes = sorted(boxes, key=lambda xy: (-xy[0], xy[1]))

    a = 0  # single box components
    b = 0  # multibox components
    c = 0  # free corner boxes
    for i in range(len(boxes)):
        if is_only_in_group(boxes, i):
            a += 1
        elif is_first_in_group(boxes, i):
            b += 1
        if is_trailing_q_corner_box(boxes, i):
            c += 1
    d = sum(kappa) - sum(nu) - a - 2 * b - c + 2

    return a, b, c, d


def test_one_var_jq(n=10):
    for kappa in Partition.all(n, strict=True):
        for nu in Partition.shifted_ribbon_complements(kappa):
            print('kappa =', kappa, 'nu =', nu)
            a, b, c, d = shifted_ribbon_q_params(kappa, nu)
            t = Polynomial.x(1)
            expected = (t + beta)**(d - 2) * 2**a * t**(a + b + c) * (2 * t + beta)**b
            f = jq(1, kappa, nu).polynomial()
            assert expected == f


def shifted_ribbon_p_params(kappa, nu):
    boxes = Partition.shifted_shape(kappa, nu)
    boxes = sorted(boxes, key=lambda xy: (-xy[0], xy[1]))

    a = 0  # single box components
    b = 0  # multibox components
    c = 0  # free corner boxes
    for i in range(len(boxes)):
        first_component = is_first_component(boxes, i)
        first_component_on_diagonal = boxes[0][0] == boxes[0][1]
        first_component_and_on_diagonal = first_component and first_component_on_diagonal

        if is_only_in_group(boxes, i) and not first_component_and_on_diagonal:
            a += 1
        elif is_first_in_group(boxes, i) and not first_component_and_on_diagonal:
            b += 1

        if is_trailing_p_corner_box(boxes, i):
            c += 1
    d = sum(kappa) - sum(nu) - a - 2 * b - c + 1

    return a, b, c, d


def test_one_var_jp(n=10):
    for kappa in Partition.all(n, strict=True):
        for nu in Partition.shifted_ribbon_complements(kappa):
            a, b, c, d = shifted_ribbon_p_params(kappa, nu)
            print('kappa =', kappa, 'nu =', nu)
            t = Polynomial.x(1)
            expected = (t + beta)**(d - 1) * 2**a * t**(a + b + c) * (2 * t + beta)**b
            f = jp(1, kappa, nu).polynomial()

            # print()
            # print('a =', a, 'b =', b, 'c =', c, 'd =', d)
            # Partition.print_shifted(kappa, nu)
            # print(expected, '==', f)
            assert expected == f


def test_hat_c(n=10):
    for nu in Partition.all(n, strict=True):
        for lam in Partition.shifted_ribbon_complements(nu):
            ans = Tableau.shifted_setvalued_copieri(nu, lam, diagonal_primes=True)

            kappa = Partition.trim((Partition.get(nu, 1) + 2,) + lam)
            a, b, c, d = shifted_ribbon_p_params(kappa, nu)

            print('kappa =', kappa, 'nu =', nu, 'lambda =', lam)
            try:
                u = Polynomial.monomial(1)
                f = (2 * u + beta)**b * 2**a * (u + beta)**(a + b + c - 1)
                g = 0
                for r in ans:
                    s = Partition.get(kappa, 1) - r
                    g += beta**(sum(kappa) - sum(nu) - s) * len(ans[r]) * u**(s - d)
                assert f == g
            except:
                print()
                print('a =', a, 'b =', b, 'c =', c, 'd =', d)
                Partition.print_shifted(kappa, nu)
                print()
                print(f, '==', g)
                print()
                # print(ans)
                input('\n?\n')


def test_hat_z(n=10):
    for kappa in Partition.all(n, strict=True):
        for nu in Partition.shifted_ribbon_complements(kappa):
            a, b, c, d = shifted_ribbon_p_params(kappa, nu)
            expansion = jp_expansion(jp(1, kappa, nu))
            print('kappa =', kappa, 'nu =', nu, ':', expansion)
            try:
                if sum(nu) == sum(kappa):
                    continue

                s = Partition.trim((sum(kappa) - sum(nu),))
                assert expansion[s] == 2**(a + b)
                assert expansion[()] == 0

                s = Partition.trim((sum(kappa) - sum(nu) - 1,))
                if sum(nu) + 1 == sum(kappa):
                    expected = Polynomial.zero()
                else:
                    expected = 2**(a + b) * (a + b + c - 1) + (2**(a + b - 1) * b if b != 0 else 0)
                    expected *= -beta
                assert expansion[s] == expected
            except:
                print()
                print('a =', a, 'b =', b, 'c =', c, 'd =', d)
                print()
                print(expansion, ':', expansion[s], '==', expected)
                print()
                traceback.print_exc()
                input('?')


def test_hat_y(n=10):
    for kappa in Partition.all(n, strict=True):
        for nu in Partition.shifted_ribbon_complements(kappa):
            a, b, c, d = shifted_ribbon_q_params(kappa, nu)
            expansion = jq_expansion(jq(1, kappa, nu))
            print('kappa =', kappa, 'nu =', nu, ':', expansion)
            try:
                if sum(nu) == sum(kappa):
                    continue

                s = Partition.trim((sum(kappa) - sum(nu),))
                assert expansion[s] == 2**(a + b - 1)
                assert expansion[()] == 0

                s = Partition.trim((sum(kappa) - sum(nu) - 1,))
                if sum(nu) + 1 == sum(kappa):
                    expected = Polynomial.zero()
                else:
                    expected = 2 * a + 3 * b + 2 * c - 3
                    if a + b - 2 == -1:
                        assert expected % 2 == 0
                        expected //= 2
                    elif a + b - 2 > 0:
                        expected *= 2**(a + b - 2)
                    expected *= -beta
                assert expansion[s] == expected
            except:
                print()
                print('a =', a, 'b =', b, 'c =', c, 'd =', d)
                print()
                print(expansion, ':', expansion[s], '==', expected)
                print()
                traceback.print_exc()
                input('?')


def test_hat_b(n=10):
    for nu in Partition.all(n, strict=True):
        for lam in Partition.shifted_ribbon_complements(nu):
            ans = Tableau.shifted_setvalued_copieri(nu, lam, diagonal_primes=False)

            for e in [1, 2, 3, 4]:
                kappa = Partition.trim((Partition.get(nu, 1) + e,) + lam)
                a, b, c, d = shifted_ribbon_q_params(kappa, nu)

                print('kappa =', kappa, 'nu =', nu, 'lambda =', lam)
                try:
                    u = Polynomial.monomial(1)
                    if e == 1 and Partition.get(lam, 1) == Partition.get(nu, 1) and sum(lam) > 0:
                        f = (2 * u + beta)**(b - 1) * 2**a * (u + beta)**(a + b + c)
                    elif e == 1:
                        f = (2 * u + beta)**b * 2**(a - 1) * (u + beta)**(a + b + c - 1)
                    else:
                        f = (2 * u + beta)**(b - 1) * 2**a * (u + beta)**(a + b + c - 1)

                    g = 0
                    for r in ans:
                        s = Partition.get(kappa, 1) - r
                        g += beta**(sum(kappa) - sum(nu) - s) * len(ans[r]) * u**(s - d + int(e == 1))
                    assert f == g
                except:
                    print()
                    print('a =', a, 'b =', b, 'c =', c, 'd =', d, 'e =', e)
                    Partition.print_shifted(kappa, nu)
                    print()
                    print(ans)
                    print(f, '==', g)
                    traceback.print_exc()
                    input('\n?\n')


def is_first_in_group(boxes, i):
    x, y = boxes[i]
    return i == 0 or (boxes[i - 1][0] != x and boxes[i - 1][1] != y)


def is_second_in_group(boxes, i):
    return not is_first_in_group(boxes, i) and is_first_in_group(boxes, i - 1)


def is_last_in_group(boxes, i):
    x, y = boxes[i]
    return i + 1 == len(boxes) or (boxes[i + 1][0] != x and boxes[i + 1][1] != y)


def is_only_in_group(boxes, i):
    return is_first_in_group(boxes, i) and is_last_in_group(boxes, i)


def is_trailing_q_corner_box(boxes, i):
    x, y = boxes[i]
    return ((x + 1, y - 1) in boxes and (x + 1, y) in boxes) or ((x + 1, y - 1) in boxes and (x, y - 1) in boxes)


def is_trailing_p_corner_box(boxes, i):
    x, y = boxes[i]
    return is_trailing_q_corner_box(boxes, i) or (x + 1 == y and (x + 1, y) in boxes) or (x == y)


def is_p_corner_box(boxes, i):
    x, y = boxes[i]
    return is_q_corner_box(boxes, i) or (x == y and (x - 1, y) in boxes)


def is_q_corner_box(boxes, i):
    x, y = boxes[i]
    return ((x, y - 1) in boxes and (x - 1, y) in boxes) or ((x + 1, y) in boxes and (x, y + 1) in boxes)


def is_first_component(boxes, i):
    u, v = boxes[i]
    for j in range(i - 1, -1, -1):
        a, b = boxes[j]
        if a != u and b != v:
            return False
        u, v = a, b
    return True


def is_first_component_with_multiple_boxes(boxes, i):
    if is_only_in_group(boxes, i):
        return False
    u, v = boxes[i]
    for j in range(i - 1, -1, -1):
        a, b = boxes[j]
        if a != u and b != v:
            return all(is_only_in_group(boxes, k) for k in range(j + 1))
        u, v = a, b
    return True


def is_last_component(boxes, i):
    u, v = boxes[i]
    for j in range(i + 1, len(boxes)):
        a, b = boxes[j]
        if a != u and b != v:
            return False
        u, v = a, b
    return True


def comarked_q_ribbons(nu, lam):
    boxes = sorted(Partition.shifted_shape(nu, lam), key=lambda xy: (-xy[0], xy[1]))
    no_multiple_boxes = all(is_only_in_group(boxes, i) for i in range(len(boxes)))

    def generate(boxes, i=0):
        if i == len(boxes):
            yield ()
            return

        first_component_with_multiple_boxes = is_first_component_with_multiple_boxes(boxes, i)
        first_in_group = is_first_in_group(boxes, i)
        last_in_group = is_last_in_group(boxes, i)
        only_in_group = is_only_in_group(boxes, i)
        corner_box = is_q_corner_box(boxes, i)
        last_box = i + 1 == len(boxes)

        choices = [1, 2] if corner_box else [2]
        if only_in_group:
            choices = [1, -1, 2, -2] if not last_box else [2] if no_multiple_boxes else [2, -2]
        elif first_in_group and not first_component_with_multiple_boxes:
            choices = [-1, 2, -2]
        elif last_in_group and not last_box:
            choices = [1, 2]

        x, y = boxes[i]
        for c in choices:
            for rest in generate(boxes, i + 1):
                tup = (((x, y), c),) + rest
                if no_multiple_boxes and any(c == -1 for _, c in tup[:-1]) and not any(c in [2, -2] for _, c in tup[:-1]):
                    continue
                yield tup

    ans = {}
    for tup in generate(boxes):
        tab = Tableau({box: val for box, val in tup})
        weight = tab.abs_sum() - tab.size()
        ans[weight] = ans.get(weight, []) + [tab]

    decomp = {key: len(val) for key, val in ans.items()}
    expected = {sum(key): val.substitute(0, 1) for key, val in jq_expansion(jq(1, nu, lam)).items()}
    return ans, decomp, expected


def q_ribbon_bijection(lam, r, tab):
    lamshape = Partition.shifted_shape(lam)
    kap = Partition.trim((Partition.get(lam, 1) + r,) + lam)
    boxes = sorted(Partition.shifted_shape(kap, lam), key=lambda xy: (-xy[0], xy[1]))
    boxes = {(i, j): 0 for (i, j) in boxes if (i, j) not in tab.boxes}
    for (i, j) in tab.boxes:
        if tab.get(i, j) == 1:
            h = max([x for (x, y) in lamshape if y == j - 1])
            if (h, j - 1) not in boxes:
                boxes[h, j - 1] = 0
            if (h, j) in tab.boxes and (h + 1, j - 1) in tab.boxes:
                v = boxes[h, j - 1]
                if (i - 1, j) in tab.boxes:
                    boxes[h, j - 1] = 1 if v == 0 else (-1, 1)
                if (i + 1, j) in tab.boxes:
                    boxes[h, j - 1] = -1 if v == 0 else (-1, 1)
    for (i, j) in tab.boxes:
        if (i + 1, j) in tab.boxes or (i, j - 1) in tab.boxes:
            continue


def test_comarked_q_ribbons(n=7):
    for nu, lam in Partition.shifted_ribbons(n):
        ans, decomp, expected = comarked_q_ribbons(nu, lam)
        print('    nu =', nu, 'lambda =', lam)
        if decomp != expected:
            Partition.print_shifted(nu, lam)
            print()
            print('  decomp:', decomp)
            print('expected:', expected)
            print()
            print(ans)
            input('\n')
        # assert decomp == expected


def comarked_p_ribbons(nu, lam):
    def generate(boxes, i=0):
        if i == len(boxes):
            yield ()
            return

        first_component = is_first_component(boxes, i)
        first_component_on_diagonal = boxes[0][0] == boxes[0][1]
        first_component_and_on_diagonal = first_component and first_component_on_diagonal

        first_in_group = is_first_in_group(boxes, i)
        last_in_group = is_last_in_group(boxes, i)
        only_in_group = is_only_in_group(boxes, i)
        corner_box = is_p_corner_box(boxes, i)
        last_box = i + 1 == len(boxes)

        choices = [1, 2] if corner_box else [2]
        if only_in_group and first_component_and_on_diagonal:
            choices = [1, 2] if not last_box else [2]
        elif only_in_group:
            choices = [1, -1, 2, -2] if not last_box else [2, -2]
        elif first_in_group and not first_component_and_on_diagonal:
            choices = [1, 2, -2]
        elif last_in_group and not last_box:
            choices = [1, 2]

        x, y = boxes[i]
        for c in choices:
            for rest in generate(boxes, i + 1):
                yield (((x, y), c),) + rest

    ans = {}
    boxes = sorted(Partition.shifted_shape(nu, lam), key=lambda xy: (-xy[0], xy[1]))
    for tup in generate(boxes):
        tab = Tableau({box: val for box, val in tup})
        weight = tab.abs_sum() - tab.size()
        ans[weight] = ans.get(weight, []) + [tab]

    decomp = {key: len(val) for key, val in ans.items()}
    expected = {sum(key): val.substitute(0, 1) for key, val in jp_expansion(jp(1, nu, lam)).items()}
    return ans, decomp, expected


def test_comarked_p_ribbons(n=7):
    for nu, lam in Partition.shifted_ribbons(n):
        ans, decomp, expected = comarked_p_ribbons(nu, lam)
        print('    nu =', nu, 'lambda =', lam)
        if decomp != expected:
            Partition.print_shifted(nu, lam)
            print()
            print('  decomp:', decomp)
            print('expected:', expected)
            print()
            print(ans)
            input('\n')
        # assert decomp == expected


# def test_simple():
#     test = []

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 1"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 1"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 1 1"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 1 1"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 1"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 1"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 2"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 2;. 2' 1'"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 1 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1'"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1 2'"), Tableau(". 0 0;. 0 1"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 2 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     for lhs, expected_forward, expected_middle, expected_backward in test:
#         forward = lhs.forward_transition(1)
#         middle, case = forward.middle_transition(1)
#         backward = middle.backward_transition(1)
#         image = lhs.transition(1)
#         post = image.transition(1)
#         print(combine_str(lhs, forward, middle, backward, image, post))
#         print(combine_str(lhs, expected_forward, expected_middle, expected_backward))
#         print()
#         print()
#         assert expected_forward == forward
#         assert expected_middle == middle
#         assert expected_backward == backward
#         assert lhs.is_semistandard()
#         assert forward.is_interstandard()
#         assert middle.is_interstandard()
#         assert backward.is_semistandard()
#         assert lhs == post


def print_transition(vst, i):
    f = vst.forward_transition(i)
    m, case = f.middle_transition(i)
    b = m.backward_transition(i)
    image = vst.transition(i)
    post = image.transition(i)
    print(combine_str(vst, f, m, b))
    print()
    print('reverse weight case:', case)
    assert image == b
    assert post == vst

def _test_small(k=9, dnp=True, verbose=False):
    n = 2
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            test = sorted(ValuedSetTableau.all(n, mu, nu, diagonal_primes=dnp))

            if verbose:
                _unseen = set(test)
                images = {}
                multiplicities = {}
                for vst in test:
                    image = vst.transition(1)
                    _unseen.discard(image)
                    multiplicities[image] = multiplicities.get(image, 0) + 1
                    key = image.unprime()
                    images[key] = images.get(key, []) + [vst]
                unseen = {}
                for vst in _unseen:
                    key = vst.unprime_diagonal()
                    unseen[key] = unseen.get(key, []) + [vst]

            seen = {}
            for vst in test:
                try:
                    f = vst.forward_transition(1)
                    m, case = f.middle_transition(1)
                    b = m.backward_transition(1)
                    image = vst.transition(1)
                    post = image.transition(1)

                    assert vst.is_semistandard()
                    assert f.is_interstandard()
                    assert m.is_interstandard()
                    assert b.is_semistandard()
                    assert image.is_semistandard()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)
                    assert dnp or not image.diagonal_primes()
                    assert vst.diagonal_primes() == image.diagonal_primes()
                    assert vst == post

                    if verbose:
                        print('\nverbose on, printing transitions:')
                        print_transition(vst, 1)
                        print(5 * '\n')
                        input('')
                        key = image
                        seen[key] = seen.get(key, []) + [vst]
                        assert len(seen[key]) == 1
                except:
                    print('\nvst:')
                    print_transition(vst, 1)
                    print_transition(image, 1)
                    print_transition(post, 1)
                    print_transition(post.transition(1), 1)
                    print(5 * '\n')
                    input('')

                    if verbose:
                        for preimage in seen[key]:
                            print_transition(preimage, 1)
                        print('preimages:', len(seen[key]))
                        print()
                        print('alternatives:')
                        alts = images.get(key.unprime(), [])
                        for u in alts:
                            if u not in seen[key]:
                                print_transition(u, 1)
                        print()
                        print('unseen:')
                        ukey = image.unprime_diagonal()
                        uns = unseen.get(ukey, [])
                        for u in uns:
                            print_transition(u, 1)

                    print(5 * '\n')

                    traceback.print_exc()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)
                    assert vst == post


def test_small_p(k=8, verbose=False):
    _test_small(k, False, verbose)


def test_small_q(k=8, verbose=False):
    _test_small(k, True, verbose)


def test_interstandard(k=8):
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            source = sorted(ValuedSetTableau.all(2, mu, nu, diagonal_primes=True))
            target = sorted(ValuedSetTableau.all_interstandard(mu, nu))

            for vst in source:
                try:
                    image = vst.forward_transition(1)
                    post = image.backward_transition(1)
                    assert vst.is_semistandard()
                    assert image.is_interstandard()
                    assert post.is_semistandard()
                    assert vst == post
                except:
                    print('\nswap direction:\n')
                    print(combine_str(vst, image, post))
                    traceback.print_exc()
                    input('\n?\n')

            for vst in target:
                try:
                    image = vst.backward_transition(1)
                    post = image.forward_transition(1)
                    assert vst.is_interstandard()
                    assert image.is_semistandard()
                    assert post.is_interstandard()
                    assert vst == post
                except:
                    print('\nunswap direction:\n')
                    print(combine_str(vst, image, post))
                    traceback.print_exc()
                    input('\n?\n')


def test_tau_exceptional_case(k=9):
    n = 2
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            test = sorted(ValuedSetTableau.all(n, mu, nu))

            for vst in test:
                i = [i for (i, j, _) in vst.tableau if i == j]
                if len(i) != 1:
                    continue
                i = i[0]
                #if vst.tableau.get(i, i) not in [-1, 2]:
                #    continue
                f = vst.forward_transition(1)
                m, case = f.middle_transition(1)
                b = m.backward_transition(1)
                image = vst.transition(1)

                assert f.tableau.get(i, i) == vst.tableau.get(i, i)
                if m.tableau.get(i, i) != vst.tableau.get(i, i) == 2:
                    assert m.tableau.get(i, i) == 1
                    assert (i - 1, i) not in vst.tableau.boxes or m.tableau.get(i - 1, i) == -1
                    if (i - 1, i) in vst.tableau.boxes:
                        print_transition(vst, 1)
    #                    input('\n\n\n')
                        if (i, i + 1) in vst.tableau.boxes and (i - 1, i + 1) in vst.tableau.boxes:
                            assert f.tableau.get(i, i + 1) == 2
                            assert (f.tableau.get(i - 1, i) == -1 and f.tableau.get(i - 1, i + 1) == -2) or (f.tableau.get(i - 1, i) == -2 and f.tableau.get(i - 1, i + 1) == 1 and not f.grouping.get(i, i))
                        if (i, i + 1) not in vst.tableau.boxes and (i - 1, i + 1) in vst.tableau.boxes:
                            assert (f.tableau.get(i - 1, i) == -2 and f.tableau.get(i - 1, i + 1) in [1, 2]) or (f.tableau.get(i - 1, i) == -1 and f.tableau.get(i - 1, i + 1) == -2)

def test_tau_two_box_case(k=9):
    n = 2
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            test = sorted(ValuedSetTableau.all(n, mu, nu))

            for vst in test:
                i = [i for (i, j, _) in vst.tableau if i == j]
                if len(i) != 2:
                    continue
                i = min(i)
                if vst.tableau.get(i, i) != -1 or vst.tableau.get(i + 1, i + 1) != 2:
                    continue
                f = vst.forward_transition(1)
                if f.tableau.get(i, i + 1) != 1:
                    continue
                m, case = f.middle_transition(1)
                b = m.backward_transition(1)
                image = vst.transition(1)

                print_transition(vst, 1)
                input('\n\n\n')

