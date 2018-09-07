from schubert import InvSchubert, FPFSchubert
from schubert import x as x_var, one as one_var
import subprocess
import os


class Pipedream:

    DIRECTORY = '/Users/emarberg/Dropbox/involution-words/Pipedreams/examples/'

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

    def __repr__(self):
        s = []
        for i in range(1, self.n + 1):
            s += [(self.n + 1 - i) * ['.']]
        for (i, j) in self.crossings:
            s[i - 1][j - 1] = '+'
        return '\n'.join(' '.join(row) for row in s)

    def __hash__(self):
        return hash(tuple(sorted(self.crossings)))

    def __eq__(self, other):
        assert type(other) == Pipedream
        return self.crossings == other.crossings

    def fpf_involution_ladder_moves(self):
        for (i, j) in self.crossings:
            if (i, j + 1) in self.crossings:
                continue
            x = i - 1
            while (x, j) in self.crossings and (x, j + 1) in self.crossings:
                x = x - 1
            if x == 0:
                continue
            if (x, j) in self.crossings:
                p = any((x - d, j + 1 + d) in self.crossings for d in range(x))
                q = any((x - 1 - d, j + 1 + d) in self.crossings for d in range(x - 1))
                r = any((x - d, j + 2 + d) in self.crossings for d in range(x))
                s = any((x - d, j - 1 + d) in self.crossings for d in range(x))
                t = any((x - 1 - d, j - 1 + d) in self.crossings for d in range(x - 1))
                if j > 1 and not (p or q or r or s or t):
                    yield Pipedream((self.crossings - {(i, j)}) | {(x, j - 1)})
            elif x > j + 1:
                yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

    def upper_fpf_involution_ladder_interval(self):
        level = {self}
        while level:
            new_level = set()
            for dream in level:
                yield dream
                for e in dream.fpf_involution_ladder_moves():
                    new_level.add(e)
            level = new_level

    def involution_ladder_moves(self):
        for (i, j) in self.crossings:
            if (i, j + 1) in self.crossings:
                continue
            x = i - 1
            while (x, j) in self.crossings and (x, j + 1) in self.crossings:
                x = x - 1
            if x == 0 or x < j + 1:
                continue
            if (x, j) in self.crossings:
                #   j j+1
                #   . .
                # x + . .
                #   + +
                #   + +
                #   + +
                # i +
                p = any((x - d, j + 1 + d) in self.crossings for d in range(x))
                q = any((x - 1 - d, j + 1 + d) in self.crossings for d in range(x - 1))
                r = any((x - d, j + 2 + d) in self.crossings for d in range(x))
                s = any((x - 1 - d, j + d) in self.crossings for d in range(x - 1))
                if p or q or r or s:
                    continue
            yield Pipedream((self.crossings - {(i, j)}) | {(x, j + 1)})

    def upper_involution_ladder_interval(self):
        level = {self}
        while level:
            new_level = set()
            for dream in level:
                yield dream
                for e in dream.involution_ladder_moves():
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
