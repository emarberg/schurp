from .vectors import Vector
from collections import defaultdict


class Word:

    @classmethod
    def fpf_involution_wiring_diagram(cls, word, labels=None):
        n = (max(word) + 1) if word else 1
        adjust = False
        if n % 2 != 0:
            n += 1
            adjust = True

        def s(c, j):
            a, b = c[0], c[1]
            a = (a + 1) if a == j else (a - 1) if (a == j + 1) else a
            b = (b + 1) if b == j else (b - 1) if (b == j + 1) else b
            return (a, b) if a < b else (b, a)

        cyc = [(i, i + 1) for i in range(1, n, 2)]
        for a in word:
            cyc = [s(c, a) for c in cyc]

        def match(c, d):
            return c[0] < d[0] < c[1]

        bools = [any(match(cyc[i], c)for c in cyc) for i in range(len(cyc))]
        assert len(bools) == n // 2

        last = ''.join([
            '╰───╯   ' if b else '█   █   '
            for b in bools
        ])

        lines = cls._wiring_diagram_helper(word, labels, []).split('\n')
        lines = lines[:-3] + [
            (n // 2) * '│   │   ',
            (n // 2) * '│   │   ',
            last
        ] + lines[-3:]

        if adjust:
            m = 4 * n - 4
            for i in range(1, len(lines) - 6):
                lines[i] = lines[i][:m] + ('│   ' if (i - 1) % 4 else '█   ') + lines[i][m:]
            lines[-2] += '   ' + str(n)

        return '\n'.join(lines)

    @classmethod
    def involution_wiring_diagram(cls, word, labels=None):
        commutations = {}
        oneline = [i for i in range(1, max(word) + 2)] if word else []
        for index, i in enumerate(word):
            newline = oneline[:]
            newline = [j + (1 if j == i else -1 if j == i + 1 else 0) for j in newline]
            newline[i - 1], newline[i] = newline[i], newline[i - 1]
            if newline == oneline:
                commutations[index + 1] = i
                newline = oneline[:]
                newline[i - 1], newline[i] = newline[i], newline[i - 1]
            oneline = newline
        return cls._wiring_diagram_helper(word, labels, commutations)

    @classmethod
    def wiring_diagram(cls, word, labels=None):
        return cls._wiring_diagram_helper(word, labels, [])

    @classmethod
    def _wiring_diagram_helper(cls, word, labels, commutations):
        assert all(i is None or i > 0 for i in word)

        a = '\u2502'  # |
        b = '\u2588'  # box
        c = '\u2572'  # \
        d = '\u2571'  # /

        def baseline(n, i):
            return [
                n * [a + '   '],
                n * [a + '   '] + [str(i) + (': ' + str(labels[i]) if labels else '')],
                n * [a + '   '],
                n * [b + '   '],
            ]

        def switch(n, i, index):
            base = baseline(n, index)

            for jndex in commutations:
                if jndex > index:
                    j = commutations[jndex]
                    base[0][j - 1] = '    '
                    base[1][j - 1] = '    '
                    base[2][j - 1] = '    '
                    base[3][j - 1] = '    '
                    base[0][j] = '    '
                    base[1][j] = '    '
                    base[2][j] = '    '
                    base[3][j] = '    '

            if i is not None:
                if index in commutations:
                    base[0][i - 1] = a + '   '
                    base[1][i - 1] = a + '   '
                    base[2][i - 1] = a + '   '
                    base[3][i - 1] = '\u2570' + 3 * '\u2500' + '\u256f'

                    base[0][i] = a + '   '
                    base[1][i] = a + '   '
                    base[2][i] = a + '   '
                    base[3][i] = '   '

                else:
                    base[0][i - 1] = ' ' + c + ' ' + d
                    base[1][i - 1] = '  \u2573 '
                    base[2][i - 1] = ' ' + d + ' ' + c

                    base[0][i] = '    '
                    base[1][i] = '    '
                    base[2][i] = '    '
            return [''.join(bits) for bits in base]

        filtered = [i for i in word if i is not None]
        n = max(filtered) + 1 if filtered else 1
        numbers = '   '.join([str(i) for i in range(1, n + 1)])
        bullets = '   '.join(n * [b]) + '   '
        lines = ['', bullets]
        for index, i in enumerate(reversed(word)):
            lines += switch(n, i, len(word) - index)
        lines += ['', numbers]
        return '\n'.join(lines) + '\n'

    @classmethod
    def f_crystal_operator(cls, i, *args, **kwargs):
        w = ()
        for subword in args:
            w += subword

        s = [(j, a) for j, a in enumerate(w) if a in [i, i + 1]]
        ell = len(s) + 1
        while len(s) < ell:
            ell = len(s)
            for x in range(len(s) - 1):
                (j, a) = s[x]
                (k, b) = s[x + 1]
                if a > b:
                    s = s[:x] + s[x + 2:]
                    break
        s = [p for p in s if p[-1] == i]
        if len(s) == 0:
            return None
        j, _ = s[-1]
        assert w[j] == i
        w = w[:j] + (i + 1,) + w[j + 1:]

        if len(args) <= 1 and kwargs.get('unpack', True):
            return w

        ans = ()
        for a in args:
            ans += (w[:len(a)],)
            w = w[len(a):]
        return ans

    @classmethod
    def e_crystal_operator(cls, i, *args, **kwargs):
        w = ()
        for subword in args:
            w += subword

        s = [(j, a) for j, a in enumerate(w) if a in [i, i + 1]]
        ell = len(s) + 1
        while len(s) < ell:
            ell = len(s)
            for x in range(len(s) - 1):
                (j, a) = s[x]
                (k, b) = s[x + 1]
                if a > b:
                    s = s[:x] + s[x + 2:]
                    break
        s = [p for p in s if p[-1] == i + 1]
        if len(s) == 0:
            return None
        j, _ = s[0]
        assert w[j] == i + 1
        w = w[:j] + (i,) + w[j + 1:]

        if len(args) <= 1 and kwargs.get('unpack', True):
            return w

        ans = ()
        for a in args:
            ans += (w[:len(a)],)
            w = w[len(a):]
        return ans

    @classmethod
    def descent_set(cls, v):
        return set(i + 1 for i in range(len(v) - 1) if v[i] > v[i + 1])

    @classmethod
    def ascents(cls, v):
        return len([i for i in range(len(v) - 1) if v[i] <= v[i + 1]])

    @classmethod
    def peaks(cls, v):
        s = [i for i in range(1, len(v) - 1) if v[i - 1] <= v[i] >= v[i + 1]]
        return (1 + cls.peaks(v[s[0] + 1:])) if s else (1 if v else 0)

    @classmethod
    def increasing_zeta(cls, word):
        return 1 if all(word[i] < word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def decreasing_zeta(cls, word):
        return 1 if all(word[i] > word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_increasing_zeta(cls, word):
        return 1 if all(word[i] <= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_decreasing_zeta(cls, word):
        return 1 if all(word[i] >= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def right_weakly_unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def left_weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def quasisymmetrize(cls, word, zeta):
        if len(word) == 0:
            return Vector({(): 1})
        ans = defaultdict(int)
        for i in range(1, len(word) + 1):
            pre = word[:i]
            sub = word[i:]
            z = zeta(pre)
            if z != 0:
                for k, v in cls.quasisymmetrize(sub, zeta).items():
                    ans[(i,) + k] += z * v
        ans = Vector({k: v for k, v in ans.items() if v != 0})
        return ans

