from permutations import Permutation
from vectors import Vector
import schubert


def raising(v, n):
    if type(v) == Permutation:
        return raising(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            ans += Vector({
                w * Permutation.s_i(i): c * i
                for i in range(1, n)
                if w(i) < w(i + 1)
            })
        return ans


def lowering(v, n):
    def a(w, i, j):
        c = w.code()
        d = (w * Permutation.t_ij(i, j)).code()
        c += (len(d) - len(c)) * (0,)
        d += (len(c) - len(d)) * (0,)
        ans = 0
        for i in range(len(c)):
            ans += abs(c[i] - d[i])
        return ans

    if type(v) == Permutation:
        return lowering(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            ans += Vector({
                w * Permutation.t_ij(i, j): c * a(w, i, j)
                for i in range(1, n)
                for j in range(i + 1, n + 1)
                if w(j) < w(i) and not any(w(j) < w(t) < w(i) for t in range(i + 1, j))
            })
        return ans


def diagonal(v, n):
    top = Permutation.longest_element(n).length()
    if type(v) == Permutation:
        return diagonal(Vector({v: 1}), n)
    elif type(v) == Vector:
        return Vector({
            w: c * (2 * w.length() - top)
            for w, c in v.items()
        })


def i_raising(v, n):
    if type(v) == Permutation:
        return i_raising(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            for i in range(1, n):
                s = Permutation.s_i(i)
                e = (i + 1) // 2
                if i == w(i) < w(i + 1) == i + 1:
                    ans += Vector({w * s: c * e})
                elif w(i) < w(i + 1):
                    ans += Vector({s * w * s: c * e})
        return ans


def i_diagonal(v, n):
    top = Permutation.longest_element(n).involution_length()
    if type(v) == Permutation:
        return i_diagonal(Vector({v: 1}), n)
    elif type(v) == Vector:
        return Vector({
            w: c * (2 * w.involution_length() - top)
            for w, c in v.items()
        })


def i_lowering(v, n):
    def a(w, i, j):
        c = w.involution_code()
        d = (w.inverse_tau_ij(i, j)).involution_code()
        c += (len(d) - len(c)) * (0,)
        d += (len(c) - len(d)) * (0,)
        ans = 0
        for i in range(len(c)):
            ans += abs(c[i] - d[i])
        return ans

    if type(v) == Permutation:
        return i_lowering(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            ans += Vector({
                w.inverse_tau_ij(i, j): c * a(w, i, j)
                for i in range(1, n)
                for j in range(i + 1, n + 1)
                if w.inverse_tau_ij(i, j).involution_length() == w.involution_length() - 1
            })
        return ans


def f_raising(v, n):
    if type(v) == Permutation:
        return f_raising(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            for i in range(1, n):
                s = Permutation.s_i(i)
                e = i
                if w(i) < w(i + 1):
                    ans += Vector({s * w * s: c * e})
        return ans


def f_diagonal(v, n):
    top = Permutation.longest_element(n).fpf_involution_length()
    if type(v) == Permutation:
        return f_diagonal(Vector({v: 1}), n)
    elif type(v) == Vector:
        return Vector({
            w: c * (2 * w.fpf_involution_length() - top)
            for w, c in v.items()
        })


def f_lowering(v, n):
    def a(w, i, j):
        t = Permutation.t_ij(i, j)
        c = w.fpf_involution_code()
        d = (t * w * t).fpf_involution_code()
        c += (len(d) - len(c)) * (0,)
        d += (len(c) - len(d)) * (0,)
        ans = 0
        for i in range(len(c)):
            ans += abs(c[i] - d[i])
        print('\n  . . . %s\n' % ans)
        return ans

    if type(v) == Permutation:
        return f_lowering(Vector({v: 1}), n)
    elif type(v) == Vector:
        ans = Vector()
        for w, c in v.items():
            ans += Vector({
                Permutation.t_ij(i, j) * w * Permutation.t_ij(i, j): c * a(w, i, j)
                for i in range(1, n)
                for j in range(i + 1, n + 1)
                if (Permutation.t_ij(i, j) * w * Permutation.t_ij(i, j)).fpf_involution_length() == w.fpf_involution_length() - 1
            })
        return ans
