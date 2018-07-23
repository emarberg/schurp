from permutations import Permutation
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

    def is_symmetric(self):
        return all((j, i) in self.crossings for (i, j) in self.crossings)

    def is_atomic(self):
        if not all((j, i) in self.crossings for (i, j) in self.crossings if i < j):
            return False
        dream = Pipedream([(i, j) for (i, j) in self.crossings if i >= j])
        z = self.permutation()
        w = dream.permutation()
        return w in z.get_atoms()

    def is_fpf_atomic(self):
        if not all((i - 1, i - 1) in self.crossings for (i, j) in self.crossings if i == j > 1):
            return False
        if not self.is_symmetric():
            return False
        dream = Pipedream([(i, j) for (i, j) in self.crossings if i > j])
        z = self.permutation()
        w = dream.permutation()
        return w in z.get_fpf_atoms()

    def atomic_part(self):
        return Pipedream({(i, j) for (i, j) in self.crossings if i >= j})

    def fpf_atomic_part(self):
        return Pipedream({(i, j) for (i, j) in self.crossings if i > j})

    def inv_monomial(self):
        if not self.is_atomic():
            return 0
        ans = one_var()
        for (i, j) in self.crossings:
            if i == j:
                ans *= x_var(i)
            elif i > j:
                ans *= x_var(i) + x_var(j)
        return ans

    def fpf_monomial(self):
        if not self.is_fpf_atomic():
            return 0
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

    def permutation(self):
        return Permutation.from_word(self.word())

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
    def test_involutions(cls, n):
        for i, perm in enumerate(Permutation.involutions(n)):
            print('. . .', i, perm)
            dreams = [Pipedream.from_word(*e) for e in perm.get_pipe_dreams()]
            test_s = sum([dream.inv_monomial() for dream in dreams])
            s = InvSchubert.get(perm)
            if s != test_s:
                return False
        return True

    @classmethod
    def test_fpf_involutions(cls, n):
        for i, perm in enumerate(Permutation.fpf_involutions(n)):
            print('. . .', i, perm)
            dreams = [Pipedream.from_word(*e) for e in perm.get_pipe_dreams()]
            test_s = sum([dream.fpf_monomial() for dream in dreams])
            s = FPFSchubert.get(perm)
            if s != test_s:
                return False
        return True

    @classmethod
    def generate_involutions(cls, n):
        for perm in Permutation.involutions(n):
            if perm(n) != n:
                Pipedream.save_involution(perm)

    @classmethod
    def save_involution(cls, perm):
        if len(perm) == 0:
            return
        dreams = [Pipedream.from_word(*e) for e in perm.get_pipe_dreams()]
        dreams = [e for e in dreams if e.is_atomic()]
        test_s = sum([dream.inv_monomial() for dream in dreams])
        s = InvSchubert.get(perm)

        filename = ''.join(map(str, perm.oneline))
        directory = cls.DIRECTORY + 'inv/' + (s != test_s) * '_' + str(len(dreams)) + '_' + filename + '/'
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
        lines += [e.atomic_part().tikz_simple() for e in dreams]
        lines += ['']
        lines += ['\\[w = %s\\]' % ''.join(map(str, perm.oneline))]
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
    def generate_fpf_involutions(cls, n):
        for perm in Permutation.fpf_involutions(n):
            if not (perm(n) == n - 1 and perm(n - 1) == n):
                Pipedream.save_fpf_involution(perm)

    @classmethod
    def save_fpf_involution(cls, perm):
        if len(perm) == 0:
            return
        dreams = [Pipedream.from_word(*e) for e in perm.get_pipe_dreams()]
        dreams = [e for e in dreams if e.is_fpf_atomic()]
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
        lines += [e.fpf_atomic_part().tikz_simple() for e in dreams]
        lines += ['']
        lines += ['\\[w = %s\\]' % ''.join(map(str, perm.oneline))]
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
