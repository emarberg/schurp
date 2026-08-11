from words import mixed_insert
from stable.tableaux import Tableau
from num2words import num2words


# from stable.partitions import *
# mus = list(Partition.all(6, strict=True))
# for n in range(1, 6):
#     for mu in mus:
#         document(n, mu)
#        print(mu, n)


def svdtab(max_entry, mu):
    ans = Tableau.setvalued_decomposition_tableaux(max_entry, mu)
    return sorted(ans, key=lambda t: (defect(t),) + content(t))


def defect(t):
    return len(t) - len(t.boxes)


def content(t):
    ans = []
    for (_, _, v) in t:
        for a in v:
            while a > len(ans):
                ans += [0]
            ans[a - 1] += 1
    return tuple(ans)




def tex(self, FRENCH=False):
    rows = []
    for i in range(1, self.max_row() + 1):
        row = []
        for j in range(1, self.max_column() + 1):
            v = self.get(i, j)
            if v is not None:
                v = (v,) if type(v) != tuple else v
                delim = '' if max(v) < 10 else ','
                v = delim.join(map(lambda x: str(x) if x > 0 else (str(-x) + "'"), v))
            row += [str(v) if v is not None else '\\none']
        rows += [' & '.join(row)]
    return '\\begin{ytableau}' + ' \\\\ '.join(reversed(rows) if FRENCH else rows) + '\\end{ytableau}'


def dmi(t):
    d = t.distribute()
    return [mixed_insert(x.reverse_row_reading_word())[0] for x in d]


def umi(t):
    e = dmi(t)
    return Tableau.union(*e)


def document(max_entry, mu):
    size = '5em'

    s = []
    s += ['\\documentclass{article}[12pt]']
    s += ['\\usepackage{fullpage}']
    s += ['\\usepackage{hyperref}']
    s += ['\\usepackage{ytableau}']
    s += ['\\ytableausetup{boxsize=%s, aligntableaux=center}' % size]
    s += ['']
    s += ['\\begin{document}']
    s += ['']
    s += ['\\noindent Shape: $\\lambda = %s$' % str(mu)]
    s += ['\\[\\]']
    s += ['Entries: nonempty subsets of $[%s]$' % str(max_entry)]
    s += ['\\[\\]']
    s += ['What is being computed: take a set-valued decomposition tableau $T \\in \\mathsf{SetDecTab}_n(\\lambda)$, form all of its distributions, take the reverse row reading words of these, apply mixed insertion, then take the union of the resulting shifted tableaux, which all have shape $\\lambda$.']
    s += ['']
    s += ['\\tableofcontents']
    s += ['']

    tabs = svdtab(max_entry, mu)
    index = -1
    for i, t in enumerate(tabs):
        if index < defect(t):
            index = defect(t)
            s += ['\\section{Defect %s}' % num2words(index)]
            s += ['']
            s += ['\\ ']
            s += ['']
        s += ['\\begin{itemize}']
        s += ['\\item[(%s)] $%s \\mapsto %s$' % (i + 1, tex(t), tex(umi(t)))]
        s += ['\\end{itemize}']
        s += ['']

    s += ['']
    s += ['\\end{document}']
    s = '\n'.join(s)

    delim = ',' if mu and max(mu) >= 10 else ''
    mustr = delim.join([str(m) for m in mu])
    fname = 'tex/svmixed_mu%s_n%s.tex' % (mustr, max_entry) 
    with open(fname, 'w') as f:
        f.write(s)

