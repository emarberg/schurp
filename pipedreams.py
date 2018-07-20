import subprocess


class Pipedream:

    DIRECTORY = '/Users/emarberg/Dropbox/involution-words/Pipedreams/examples/'

    def __init__(self, crossings):
        self.crossings = {(i, j) for (i, j) in crossings}
        self.n = max([i + j for (i, j) in self.crossings]) if self.crossings else 1

    def __repr__(self):
        return '<Pipedream for %s>' % ''.join([str(s) for s in self.word()])

    def tikz(self):
        s = []
        s += ['\\begin{center}']
        s += ['\\begin{tikzpicture}']
        for i in range(1, self.n + 1):
            i_ = self.n - i
            s += ['\\node at (0.25, %s) {%s};' % (i_, i)]
            s += ['\\node at (%s, %s-0.20) {%s};' % (i, self.n, i)]
        for i in range(1, self.n + 1):
            for j in range(1, self.n + 1):
                if i + j > self.n + 1:
                    continue
                j_ = self.n - j
                if (j, i) not in self.crossings:
                    s += ['\\draw [domain=270:360] plot ({%s-0.5+0.5* cos(\\x)}, {%s+0.5+0.5*sin(\\x)});' % (i, j_)]
                if (j, i) not in self.crossings and i + j <= self.n:
                    s += ['\\draw [domain=90:180] plot ({%s+0.5+0.5*cos(\\x)}, {%s-0.5+0.5*sin(\\x)});' % (i, j_)]
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
                    word += [j + 1 - i]
        return tuple(word)

    def filename_prefix(self):
        return ''.join(map(str, self.word()))

    def filename(self):
        return self.filename_prefix() + '.tex'

    def save(self):
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
        subprocess.run(["pdflatex", self.filename()])
        subprocess.run(["mv", self.filename(), self.DIRECTORY])
        subprocess.run(["mv", self.filename_prefix() + '.pdf', self.DIRECTORY])
        subprocess.run(["rm", self.filename_prefix() + '.log'])
        subprocess.run(["rm", self.filename_prefix() + '.aux'])
