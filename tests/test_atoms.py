from signed import SignedPermutation
import subprocess


class SignedAtomsGraph:

    DIRECTORY = '/Users/emarberg/examples/stanley-type-c/atoms/'

    def tikz(self):
        pre = [
            '\\begin{center}',
            '\\begin{tikzpicture}',
        ]

        vertices = sorted({u for u, _ in self.edges} | {v for _, v in self.edges})
        labels = {v: str(self.n) + 'node' + str(i) for i, v in enumerate(vertices)}

        nodes = [
            '\\node (%s) {%s};' % (labels[v], v.tex()) for v in labels
        ]

        edges = [
            '\\draw [->] (%s) -- (%s);' % (labels[u], labels[v]) for u, v in self.edges
        ]

        post = [
            '\\end{tikzpicture}',
            '\\end{center}',
        ]

        return '\n'.join(pre + nodes + edges + post)

    @property
    def _filename(self):
        return 'atoms_graph_b%s' % self.n

    @property
    def dot_filename(self):
        return self.DIRECTORY + 'dot/' + '%s.dot' % self._filename

    @property
    def png_filename(self):
        return self.DIRECTORY + 'png/' + '%s.png' % self._filename

    def node_label(self, i):
        if i.rank != self.n + 1:
            assert SignedPermutation.standardize(i.oneline[1:]) in self.sub_atoms
        return str(i)
        # if i.rank == self.n + 1:
        #     return str(i.inverse())

        # sh = [(a, b) for (a, b) in i.get_atom_shape(i.oneline) if a + b != 0]
        # base = ' '.join([str(a) for a in range(1, self.n + 1)])

        # def index(j):
        #     return (j - 1) * 2

        # lines = []
        # for a, b in sorted(sh, key=lambda x: x[1] - x[0]):
        #     j, k = index(a), index(b)
        #     if not (lines and all(s == ' ' for s in lines[0][j:k])):
        #         lines = [(2 * self.n) * [' ']] + lines
        #     lines[0][j] = '\u256D'
        #     lines[0][k] = '\u256E'
        #     for l in range(j + 1, k):
        #         lines[0][l] = '\u2500'
        #     for l in range(1, len(lines)):
        #         lines[l][j] = '\u2502'
        #         lines[l][k] = '\u2502'
        # base = '\n'.join([''.join(l) for l in lines] + [base])

        # return str(i.inverse()) + '\n\n' + base

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"];']
        s += ['    "%s" -> "%s";' % (self.node_label(x), self.node_label(y)) for (x, y) in self.edges]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    def __init__(self, n):
        self.n = n
        self.atoms = SignedPermutation.longest_element(n).get_atoms()
        self.sub_atoms = SignedPermutation.longest_element(n - 1).get_atoms()
        self._edges = None

    @classmethod
    def queue_stanley_decomposition(cls, n, verbose=True):
        def get_shape(oneline):
            while oneline[-1] == len(oneline):
                oneline = oneline[:-1]
            oneline = oneline[1:]
            ans = []
            while oneline:
                for i in range(len(oneline)):
                    a = oneline[i]
                    if i == 0 and a < 0:
                        ans += [(a, -a)]
                        oneline = oneline[1:]
                        break
                    if i + 1 >= len(oneline):
                        continue
                    b = oneline[i + 1]
                    if 0 < a < -b:
                        ans += [(a, -b)]
                        oneline = oneline[:i] + oneline[i + 2:]
                        break
            return ans

        def get_rs(x):
            r = 0
            s = n
            while r < s - 1 and (x(r + 1) == n or x(r + 1) == x(r) - 1):
                r += 1
            try:
                z = SignedPermutation(*x.oneline[r:]).reduce()
                a = SignedPermutation.longest_element(n - r - 1).get_atoms()
                assert z in a
            except:
                r += 1
                shape = get_shape(tuple(x(i) for i in range(r, n + 1)))
                s = x.inverse()(
                    max([a for a, b in shape if 0 < a < x(r) < b])
                )
            return r, s

        ans = []
        atoms = SignedPermutation.longest_element(n).get_atoms()
        preatoms = SignedPermutation.longest_element(n - 1).get_atoms()
        queue = [SignedPermutation(*((n + 1,) + x.oneline + (n,))) for x in preatoms]
        while queue:
            print('queue: %s' % len(queue))
            x = queue.pop(0).reduce()
            if x in atoms:
                print('\n*', x, 'is atom\n')
                ans += [x]
                continue

            n = x.rank
            r, s = get_rs(x)
            v = x * SignedPermutation.reflection_t(r, s, n)
            v_len = len(v)
            assert v_len + 1 == len(x)

            newline = v.oneline + (n + 1,)
            new_v = SignedPermutation(*newline)
            lhs = [new_v * SignedPermutation.reflection_t(r, i, n + 1) for i in range(r + 1, n + 2) if i != s]
            lhs = sorted([u.reduce() for u in lhs if v_len + 1 == len(u)])

            if lhs != sorted([u.reduce() for u in queue if u.reduce() in lhs]):
                queue += [x]
                continue

            yield (x.reduce(), v.reduce())
            for u in lhs:
                yield (u, v.reduce())

            queue = [u for u in queue if u.reduce() not in lhs]
            rhs = [v * SignedPermutation.reflection_s(r, r, n)]
            rhs += [v * SignedPermutation.reflection_t(i, r, n) for i in range(1, r)]
            rhs += [new_v * SignedPermutation.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
            rhs = [u.reduce() for u in rhs if v_len + 1 == len(u)]
            queue += rhs

            for u in rhs:
                yield(v.reduce(), u)

        assert sorted(ans) == sorted(atoms)

    @property
    def edges(self):
        if self._edges is None:
            self._edges = [
                (a, b)
                for (a, b) in self.queue_stanley_decomposition(self.n)
                if a.rank == self.n
            ]
        return self._edges

    def test(self):
        """Tests for cgraphs.tex"""
        def a_shape(w):
            oneline = [w(i) for i in range(1, w.rank + 1)]
            m = set()
            go = True
            while go and oneline:
                for i in range(len(oneline)):
                    go = False
                    if i + 1 < len(oneline) and oneline[i] > oneline[i + 1]:
                        a = oneline[i]
                        b = -oneline[i + 1]
                        assert 0 < a < b
                        m |= {(a, b), (-b, -a)}
                        oneline = oneline[:i] + oneline[i + 2:]
                        go = True
                        break
            assert all(c < 0 for c in oneline)
            m |= {(c, -c) for c in oneline}
            return m

        def q_shape(w):
            oneline = [w(i) for i in range(2, w.rank + 1)]
            m = set()
            go = True
            while go and oneline:
                for i in range(len(oneline)):
                    go = False
                    if i + 1 < len(oneline) and oneline[i] > oneline[i + 1]:
                        a = oneline[i]
                        b = -oneline[i + 1]
                        assert 0 < a < b
                        m |= {(a, b), (-b, -a)}
                        oneline = oneline[:i] + oneline[i + 2:]
                        go = True
                        break
            assert all(c < 0 for c in oneline)
            m |= {(c, -c) for c in oneline}
            return m

        def is_even(w):
            return w(1) < 0 or (w(1) - self.n - 1) % 2 == 0

        def up(w, i):
            if i + 2 <= w.rank:
                c, a, b = w(i), w(i + 1), w(i + 2)
                if a < b < c:
                    s = SignedPermutation.s_i(i, w.rank)
                    t = SignedPermutation.s_i(i + 1, w.rank)
                    return w * t * s

        def is_maximal(w):
            return all(up(w, i) is None for i in range(1, w.rank))

        def down(w, i):
            if i + 2 <= w.rank:
                b, c, a = w(i), w(i + 1), w(i + 2)
                if a < b < c:
                    s = SignedPermutation.s_i(i, w.rank)
                    t = SignedPermutation.s_i(i + 1, w.rank)
                    return w * s * t

        failures = 0

        # lengths
        for u, v in self.edges:
            try:
                if is_even(u):
                    assert len(v) == len(u) - 1
                else:
                    assert len(u) == len(v) - 1
            except:
                failures += 1

        # lemma 4.2(a)
        for u, v in self.edges:
            if not is_even(v):
                continue
            for i in range(2, self.n):
                uu = up(u, i)
                vv = up(v, i)
                try:
                    assert (vv is None) or ((uu, vv) in self.edges)
                except:
                    failures += 1
                    print('4.2(a)', u, '->', v, ', ', uu, '->', vv, ', ', i)

        # lemma 4.2(b)
        for v, w in self.edges:
            if not is_even(v):
                continue
            for i in range(2, self.n):
                vv = down(v, i)
                ww = down(w, i)
                try:
                    assert (vv is None) or ((vv, ww) in self.edges)
                except:
                    failures += 1
                    print('4.2(b)', v, '->', w, ', ', vv, '->', ww, ', ', i)

        q = {u for u, _ in self.edges} | {v for _, v in self.edges}
        q_odd = {u for u in q if not is_even(u)}
        q_even = q - q_odd
        atoms = {u for u in q_even if u(1) < 0 or any(a < u(1) < b == -a for a, b in q_shape(u))}
        q_even = q_even - atoms

        print('atoms:', len(atoms), ' even:', len(q_even), ' odd:', len(q_odd))

        # lemma 4.3
        for v in atoms:
            m = a_shape(v)
            if v(1) < 0:
                p = m - {(v(1), -v(1))}
                b = -v(1)
            else:
                a = v(1)
                b = [y for x, y in m if x == a][0]
                p = (m - {(a, b), (-b, -a)}) | {(-a, a)}
            try:
                assert not any(x < b < y for x, y in p)
                assert not any(x < c < y < d for x, y in p for c, d in p)
            except:
                failures += 1
            if is_maximal(v):
                pairs = sorted([(x, -y) for (x, y) in p if 0 < x < y])
                pairs = [x for pr in pairs for x in pr]
                fixed = sorted([x for (x, y) in p if x + y == 0])
                alpha = SignedPermutation(*([b] + fixed + pairs))
                try:
                    assert [alpha] == [u for (u, z) in self.edges if z == v]
                except:
                    print(alpha, ' v =', v, ' M\' =', p)
                    failures += 1

        # lemma 4.5(a) and theorem 4.8(a)
        for w in q_odd:
            b = w(1)
            for i in range(2, self.n):
                ww = up(w, i)
                if ww:
                    for a in range(1, b):
                        t = SignedPermutation.reflection_t(a, b, self.n)
                        if len(w) + 1 == len(t * w):
                            try:
                                assert up(t * w, i) == t * ww
                                assert len(ww) + 1 == len(t * ww)
                                assert (w, t * w) in self.edges
                            except:
                                print('w =', w, ' t =', t, ' w\' =', ww, ' i =', i, ' ? ', up(t * w, i), '!=', t * ww)
                                failures += 1

        # lemma 4.5(b) and theorem 4.8(b)
        for w in q_odd:
            b = w(1)
            for i in range(2, self.n):
                ww = down(w, i)
                if ww:
                    for c in range(b + 1, self.n + 1):
                        t = SignedPermutation.reflection_t(b, c, self.n)
                        if len(w) + 1 == len(t * w):
                            try:
                                assert down(t * w, i) == t * ww
                                assert len(ww) + 1 == len(t * ww)
                                assert (t * w, w) in self.edges
                            except:
                                print('w =', w, ' t =', t, ' w\' =', ww, ' i =', i, ' ? ', down(t * w, i), '!=', t * ww)
                                failures += 1

        # lemma 4.6
        for u in q_odd:
            p = q_shape(u)
            assert u(1) > 0
            m = p | {(-u(1), u(1))}
            try:
                assert not any(x < c < y < d for x, y in m for c, d in m)
                v = u * SignedPermutation.s_i(0, u.rank)
                assert v in atoms
                assert (u, v) in self.edges
            except:
                print(u, v, p, m)
                failures += 1

        # lemma 4.7
        for u in q_odd:
            p = q_shape(u)
            b = u(1)
            for a in range(1, b):
                if (-a, a) not in p:
                    continue
                t = SignedPermutation.reflection_t(a, b, u.rank)
                v = t * u
                if len(v) != len(u) + 1:
                    continue
                m = (p - {(-a, a)}) | {(a, b), (-b, -a)}
                try:
                    assert not any(x < c < y < d for x, y in m for c, d in m)
                    assert v in atoms
                    assert (u, v) in self.edges
                except:
                    print(u, v, a, b, p, m)
                    failures += 1

        print('failures:', failures)


def test_signed_atoms(nmax=6):
    for n in range(1, nmax):
        g = SignedAtomsGraph(n)
        g.test()
