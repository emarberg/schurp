from .vectors import Vector
from .polynomials import X, Y
from .partitions import Partition
from collections import defaultdict


beta = X(0)
xvar = X()
yvar = Y()


def operator_decorator(func):
    def op(vector):
        if type(vector) == tuple:
            vector = Vector({vector: 1})
        return func(vector)
    return op


def operator_shifted_add(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            row = Partition.find_shifted_inner_corner(mu, abs(index))
            if row is not None:
                ans += Vector({mu: beta * coefficient})
                continue
            row = Partition.find_shifted_outer_corner(mu, abs(index))
            if row is not None:
                nu = Partition.add(mu, row)
                ans += Vector({nu: coefficient})
        return ans
    return op


def operator_AL(n, var=yvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n + 1):
            ans += var * operator_shifted_add(i)(ans)
        return ans
    return op


def operator_AR(n, var=yvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n, -1, -1):
            ans += var * operator_shifted_add(i)(ans)
        return ans
    return op


def operator_shifted_row_remove(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            e = 0
            i = abs(index)
            row = Partition.find_shifted_inner_corner(mu, i)
            while i >= 0 and row is not None:
                i = i - 1
                mu = Partition.remove_shifted_inner_corner(mu, row)
                if mu is None:
                    break
                ans += Vector({mu: (-beta) ** e})
                e += 1
        return ans
    return op


def operator_shifted_column_remove(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            e = 0
            i = abs(index)
            row = Partition.find_shifted_inner_corner(mu, i)
            while row is not None:
                i = i + 1
                mu = Partition.remove_inner_corner(mu, row)
                ans += Vector({mu: (-beta) ** e})
                row = Partition.find_shifted_inner_corner(mu, i)
                e += 1
        return ans
    return op


def operator_R(n, var=xvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n + 1):
            ans += var * operator_shifted_row_remove(i)(ans)
        return ans
    return op


def operator_C(n, var=xvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n, 0, -1):
            ans += var * operator_shifted_column_remove(i)(ans)
        return ans
    return op


def _extract(expr, n, coefficients, differences):
    for nu, coeff in expr.items():
        coefficients[nu][n] = coeff
        if coeff != 0 and coefficients[nu].get(n - 1, 0) != 0:
            differences[nu] += [coeff - coefficients[nu][n - 1]]


def compare(max_size_of_partitions, number_of_values_of_n, op1, op2):
    for mu in Partition.all(max_size_of_partitions, strict=True):
        s = mu[0] if mu else 0
        coefficients = defaultdict(dict)
        differences = defaultdict(list)
        for n in range(s, number_of_values_of_n + s):
            lhs = op1(n)(op2(n)(mu))
            rhs = op2(n)(op1(n)(mu))
            _extract(rhs - lhs, n, coefficients, differences)

        if coefficients:
            print()
            print('mu =', mu)
            print()
            for nu in sorted(coefficients, reverse=True):
                k = max(mu[0] if mu else 0, nu[0] if nu else 0)
                if k >= number_of_values_of_n + s:
                    continue
                print('  mu =', mu, '--> nu =', nu, ': coefficient of nu in op2(op1(mu)) - op1(op2(mu)) is')
                print()
                for n in range(k, number_of_values_of_n + s):
                    c = coefficients[nu].get(n, 0)
                    print('    n =', n, ':', c)
                print()
            print()
            input('continue (press any key)')


def compare_AL_R_versus_R_AL(max_size_of_partitions, number_of_values_of_n):
    sb = 'x / (1 + %s)' % str(beta)
    for mu in Partition.all(max_size_of_partitions, strict=True):
        print()
        print('mu =', mu)
        print()
        s = mu[0] if mu else 0
        coefficients = defaultdict(dict)
        differences = defaultdict(list)
        for n in range(s, number_of_values_of_n + s):
            op_AL = operator_AL(n)  # noqa
            op_R = operator_R(n)  # noqa
            lhs = op_AL(op_R(mu))
            rhs = op_R(op_AL(mu))
            _extract(rhs - lhs, n, coefficients, differences)
        for nu in sorted(coefficients):
            k = max(mu[0] if mu else 0, nu[0] if nu else 0)
            print('  mu =', mu, '--> nu =', nu, ': coefficient of nu in R_n(AL_n(mu)) - AL_n(R_n(mu)) is')
            print()
            n = number_of_values_of_n + s - 1
            f = sum([(-beta)**i for i in range(n + 1 - k)]) * X()
            c = coefficients[nu].get(n, 0)
            e = (c - f)
            if e:
                print('    ', e, '+', sb)
            else:
                print('    ', sb)
            print()
        print()
        input('continue (press any key)')

