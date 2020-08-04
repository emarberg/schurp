from permutations import Permutation


def test_bruhat_minimal():
    def check(w, p, r, i, j):
        v = w * Permutation.t_ij(i, j)
        return p < r and i < j and w(i) < w(j) and v(p) < v(r)

    def expand(w, p, r, i, j):
        t = Permutation.t_ij
        return w * t(i, j) * t(p, r)

    for n in range(1, 5):
        for w in Permutation.all(n):
            for (p, q) in w.outer_corners():
                for i in range(1, n + 1):
                    for j in range(1, n + 1):
                        for r in range(1, n + 1):
                            if check(w, p, r, i, j):
                                choices = {i, j, p, r}
                                print(w, p, q, r, i, j, choices)
                                assert any(
                                    check(w, i2, j2, p, r2) and
                                    expand(w, p, r, i, j) == expand(w, i2, j2, p, r2)
                                    for i2 in choices
                                    for j2 in choices
                                    for r2 in choices
                                )


def test_fpf_bruhat_minimal():
    def check(w, p, r, i, j):
        v = Permutation.t_ij(i, j) * w * Permutation.t_ij(i, j)
        return p < r and i < j and w(i) < w(j) and v(p) < v(r)

    def expand(w, p, r, i, j):
        t = Permutation.t_ij
        return t(p, r) * t(i, j) * w * t(i, j) * t(p, r)

    for n in range(1, 9):
        for w in Permutation.fpf_involutions(n):
            for (p, q) in w.outer_corners():
                if p <= q:
                    continue
                for i in range(1, n + 1):
                    for j in range(1, n + 1):
                        for r in range(1, n + 1):
                            if check(w, p, r, i, j):
                                choices = {i, j, p, r, w(i), w(j), w(p), w(r)}
                                print(w, p, q, r, i, j, choices)
                                assert any(
                                    check(w, i2, j2, p, r2) and
                                    expand(w, p, r, i, j) == expand(w, i2, j2, p, r2)
                                    for i2 in choices
                                    for j2 in choices
                                    for r2 in choices
                                )
