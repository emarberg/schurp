from permutations import Permutation
from signed import SignedPermutation
from even import EvenSignedPermutation
import polynomials
from vectors import Vector
from heapq import heappush, heappop
import json
import math
import time
from pathlib import Path
import subprocess


class QPWGraph:

    def __init__(self, qpmodule, sgn=None, nbytes=8, setup=True):
        self.qpmodule = qpmodule
        self.sgn = sgn
        self.nbytes = nbytes
        self.frame = None
        self.is_cbasis_computed = False
        self.is_wgraph_computed = False

        self.wbytes = None
        self.wgraph = None
        self.wgraph_addresses = None

        assert sgn is not None or qpmodule.layer is None

        if setup:
            self.setup()
        else:
            self.is_setup = False

    def get_wgraph_size(self):
        return self._int(self.wgraph_addresses[-self.wbytes:])

    def get_wgraph_edges(self, w, include_label=False):
        assert self.is_wgraph_computed
        start = self._int(self.wgraph_addresses[w * self.wbytes:(w + 1) * self.wbytes])
        stop = self._int(self.wgraph_addresses[(w + 1) * self.wbytes:(w + 2) * self.wbytes])
        for i in range(start, stop):
            s = i * (self.sbytes + self.nbytes)
            y = self._int(self.wgraph[s:s + self.sbytes])
            if include_label:
                mu = self._int(self.wgraph[s + self.sbytes:s + self.sbytes + self.nbytes], signed=True)
                yield (y, mu)
            else:
                yield y

    def get_simple_edges(self, w):
        assert self.is_wgraph_computed
        for y in self.get_wgraph_edges(w):
            if w in set(self.get_wgraph_edges(y)):
                yield y

    def compute_wgraph(self, verbose=True):
        if not self.is_cbasis_computed:
            self.compute_cbasis()

        t0 = time.time()
        if verbose:
            print()
            print('Computing W-graph, cells, and molecules:')

        if not self.is_wgraph_computed:
            count = 0
            for w in self.qpmodule:
                for y in range(w):
                    if (self._ascents(w) | self._ascents(y)) != self._ascents(y) and self.get_cbasis_leading(y, w) != 0:
                        count += 1
                count += len({self.qpmodule.operate(w, s) for s in self.qpmodule.strict_ascents(w)})

            if verbose:
                print('* calculated size %s of edge set (%s seconds)' % (str(count), str(int(1000 * (time.time() - t0)) / 1000.0)))
                t0 = time.time()

            self.wbytes = self._bytes(count)
            self.wgraph_addresses = bytearray((self.qpmodule.size + 1) * self.wbytes)
            self.wgraph = bytearray(count * (self.sbytes + self.nbytes))

            def add_edge(start, y, mu):
                if mu != 0:
                    s = start * (self.sbytes + self.nbytes)
                    self.wgraph[s:s + self.sbytes] = y.to_bytes(self.sbytes, byteorder='big', signed=False)
                    s += self.sbytes
                    self.wgraph[s:s + self.nbytes] = mu.to_bytes(self.nbytes, byteorder='big', signed=True)
                    start += 1
                return start

            start = 0
            for w in self.qpmodule:
                self.wgraph_addresses[w * self.wbytes:(w + 1) * self.wbytes] = start.to_bytes(self.wbytes, byteorder='big', signed=False)
                for y in range(w):
                    if (self._ascents(w) | self._ascents(y)) != self._ascents(y):
                        start = add_edge(start, y, self.get_cbasis_leading(y, w))
                for y in sorted({self.qpmodule.operate(w, s) for s in self.qpmodule.strict_ascents(w)}):
                    start = add_edge(start, y, 1)
            self.wgraph_addresses[-self.wbytes:] = start.to_bytes(self.wbytes, byteorder='big', signed=False)

            if verbose:
                print('* wrote edges (%s seconds)' % str(int(1000 * (time.time() - t0)) / 1000.0))
                t0 = time.time()

        self.is_wgraph_computed = True
        self.cells = self._compute_cells(self.get_wgraph_edges)

        if verbose:
            print('* calculated cells (%s seconds)' % str(int(1000 * (time.time() - t0)) / 1000.0))
            t0 = time.time()

        self.molecules = self._compute_cells(self.get_simple_edges)

        if verbose:
            print('* calculated molecules (%s seconds)' % str(int(1000 * (time.time() - t0)) / 1000.0))
            print()

    def slow_compute_wgraph(self):
        if not self.is_cbasis_computed:
            self.compute_cbasis()

        wgraph = []
        for x in self.qpmodule:
            asc_x = self._weak_ascents(x) | self._strict_ascents(x)
            for y in self.qpmodule:
                asc_y = self._weak_ascents(y) | self._strict_ascents(y)
                if (asc_x | asc_y) != asc_y:
                    mu = self.get_cbasis_leading(x, y) + self.get_cbasis_leading(y, x)
                    if mu != 0:
                        wgraph.append((x, y, mu))
        wgraph.sort()

        edges = {}
        for x, y, mu in wgraph:
            edges[x] = edges.get(x, []) + [(y, mu)]

        return edges

    def print_wgraph(self):
        if not self.is_wgraph_computed:
            self.compute_wgraph()

        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [fontname="courier"];']

        for x in self.qpmodule:
            s += ['    "%s";' % str(self.permutation(x))]

        s += ['    "%s" -> "%s";' % (str(self.permutation(x)), str(self.permutation(y))) for x in self.qpmodule for y in self.get_wgraph_edges(x)]
        s += ['}']
        s = '\n'.join(s)

        directory = self.qpmodule.get_directory()
        Path(directory).mkdir(parents=True, exist_ok=True)

        dotfile = directory + ('wgraph.dot' if self.sgn is None else 'wgraph.unsigned.dot' if self.sgn else 'wgraph.signed.dot')
        with open(dotfile, 'w') as f:
            f.write(s)

        pngfile = directory + ('wgraph.png' if self.sgn is None else 'wgraph.unsigned.png' if self.sgn else 'wgraph.signed.png')
        subprocess.run(["dot", "-Tpng", dotfile, "-o", pngfile])

    def get_molecules_as_permutations(self):
        assert self.is_wgraph_computed
        return {tuple(sorted(self.permutation(i) for i in c)) for c in self.molecules}

    def get_cells_as_permutations(self):
        assert self.is_wgraph_computed
        return {tuple(sorted(self.permutation(i) for i in c)) for c in self.cells}

    def _compute_cells(self, edges):
        cells = []

        def strong_connect(v, stack, visited, lowlinks, onstack, index):
            visited[v] = index
            lowlink[v] = index
            index += 1
            stack.append(v)
            onstack[v] = True
            for w in edges(v):
                if visited[w] is None:
                    stack, visited, lowlinks, onstack, index = strong_connect(w, stack, visited, lowlink, onstack, index)
                    lowlink[v] = min(lowlink[v], lowlink[w])
                elif onstack[w]:
                    lowlink[v] = min(lowlink[v], visited[w])
            if lowlink[v] == visited[v]:
                cell = set()
                while True:
                    w = stack.pop()
                    cell.add(w)
                    onstack[w] = False
                    if w == v:
                        break
                cells.append(cell)
            return stack, visited, lowlinks, onstack, index

        visited = [None for _ in self.qpmodule]
        lowlink = [None for _ in self.qpmodule]
        onstack = [False for _ in self.qpmodule]
        index = 0
        stack = []

        for v in self.qpmodule:
            if visited[v] is None:
                stack, visited, lowlinks, onstack, index = strong_connect(v, stack, visited, lowlink, onstack, index)

        return cells

    def setup(self, verbose=True):
        t0 = time.time()
        if verbose:
            print()
            print(self)
            print()
            print('* computing heights', end='') # noqa

        self.heights = []
        for n in self.qpmodule:
            nht = self.qpmodule.height(n)
            if nht >= len(self.heights):
                self.heights.append(1)
            else:
                self.heights[nht] += 1

        self.cumheights = [0]
        for h in self.heights[:-1]:
            self.cumheights.append(self.cumheights[-1] + h)

        t1 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t1 - t0)))
            print('* computing size', end='') # noqa

        self.size = 0
        for i in range(len(self.heights)):
            for j in range(i + 1, len(self.heights)):
                self.size += self._space(i, j) * self.heights[i] * self.heights[j]
        assert self.size < 0.8e+10

        t2 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t2 - t1)))
            print('* computing suboffsets', end='') # noqa

        self.hbytes = self._bytes(max(self.heights))
        self.suboffsets = bytearray(self.hbytes * len(self.qpmodule))

        offset = 0
        start = 0
        for n in self.qpmodule:
            if offset > 0 and self.qpmodule.height(n) > self.qpmodule.height(n - 1):
                offset = 0
            v = self.heights[self.qpmodule.height(n)] - 1 - offset
            self.suboffsets[start:start + self.hbytes] = v.to_bytes(self.hbytes, byteorder='big')
            start += self.hbytes
            offset += 1

        t3 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t3 - t2)))
            print('* computing height offsets', end='') # noqa

        self.offsets = {}
        total_offsets = {}
        for j in range(len(self.heights)):
            self.offsets[j] = (j + 1) * [0]
            for i in range(j - 2, -1, -1):
                self.offsets[j][i] = self.offsets[j][i + 1] + self.heights[i + 1] * self._space(i + 1, j)
            total_offsets[j] = self.offsets[j][0] + self.heights[0] * self._space(0, j)

        t4 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t4 - t3)))
            print('* computing addresses', end='') # noqa

        self.abytes = self._bytes(self.size)
        self.addresses = bytearray(self.abytes * len(self.qpmodule))

        a = 0
        start = 0
        for j in self.qpmodule:
            self.addresses[start:start + self.abytes] = a.to_bytes(self.abytes, byteorder='big')
            start += self.abytes
            a += total_offsets[self.qpmodule.height(j)]

        t5 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t5 - t4)))
            print('* computing intervals', end='') # noqa

        self.sbytes = self._bytes(self.qpmodule.size)
        self.even_intervals, self.even_invleft, self.even_invright = self._compute_intervals(0)
        self.odd_intervals, self.odd_invleft, self.odd_invright = self._compute_intervals(1)

        self.gbytes = max(1, math.ceil(self.qpmodule.rank / 8.0))

        def scents(fn):
            return self._tobytearray([self._bitmap(set(fn(i))) for i in self.qpmodule], self.gbytes)

        if not self.sgn:
            self.weak_ascents = scents(self.qpmodule.weak_ascents)
            self.strict_ascents = scents(self.qpmodule.strict_ascents)
            self.descents = scents(lambda i: set(self.qpmodule.weak_descents(i)) | set(self.qpmodule.strict_descents(i)))
        else:
            self.weak_ascents = scents(self.qpmodule.weak_descents)
            self.strict_ascents = scents(self.qpmodule.strict_ascents)
            self.descents = scents(lambda i: set(self.qpmodule.weak_ascents(i)) | set(self.qpmodule.strict_descents(i)))

        t6 = time.time()
        if verbose:
            print(' %s milliseconds' % int(1000 * (t6 - t5)))
            print('* finished')
            print()

        self.is_setup = True

    def __eq__(self, other):
        assert type(other) == type(self)
        return self.sgn == other.sgn and \
            self.nbytes == other.nbytes and \
            self.heights == other.heights and \
            self.cumheights == other.cumheights and \
            self.size == other.size and \
            self.hbytes == other.hbytes and \
            self.offsets == other.offsets and \
            self.abytes == other.abytes and \
            self.sbytes == other.sbytes and \
            self.gbytes == other.gbytes and \
            self.is_setup == other.is_setup and \
            self.suboffsets == other.suboffsets and \
            self.addresses == other.addresses and \
            self.weak_ascents == other.weak_ascents and \
            self.strict_ascents == other.strict_ascents and \
            self.descents == other.descents and \
            self.even_intervals == other.even_intervals and \
            self.even_invleft == other.even_invleft and \
            self.even_invright == other.even_invright and \
            self.odd_intervals == other.odd_intervals and \
            self.odd_invleft == other.odd_invleft and \
            self.odd_invright == other.odd_invright and \
            self.is_cbasis_computed == other.is_cbasis_computed and \
            self.is_wgraph_computed == other.is_wgraph_computed and \
            self.frame == other.frame and \
            self.wgraph == other.wgraph and \
            self.wgraph_addresses == other.wgraph_addresses and \
            self.wbytes == other.wbytes

    def get_directory(self):
        if not self.sgn:
            return self.qpmodule.get_directory() + 'unsigned/'
        else:
            return self.qpmodule.get_directory() + 'signed/'

    @classmethod
    def read(cls, directory, sgn=None, nbytes=8):
        qpmodule = QPModule.read(directory)
        directory += 'unsigned/' if not sgn else 'signed/'
        wgraph = QPWGraph(qpmodule, sgn, nbytes, setup=False)
        wgraph._read()
        return wgraph

    def _read(self):
        metafile = self.get_directory() + 'wgraph.meta'
        with open(metafile, 'r') as file:
            dictionary = json.loads(file.read())
        assert self.sgn == dictionary['sgn']
        assert self.nbytes == dictionary['nbytes']
        self.heights = dictionary['heights']
        self.cumheights = dictionary['cumheights']
        self.size = dictionary['size']
        self.hbytes = dictionary['hbytes']
        self.offsets = {int(k): v for (k, v) in dictionary['offsets'].items()}
        self.abytes = dictionary['abytes']
        self.sbytes = dictionary['sbytes']
        self.gbytes = dictionary['gbytes']
        self.wbytes = dictionary['wbytes']
        self.is_setup = dictionary['is_setup']
        self.is_cbasis_computed = dictionary['is_cbasis_computed']
        self.is_wgraph_computed = dictionary['is_wgraph_computed']

        if self.is_setup:
            aux = self.get_directory() + 'aux/'
            self.suboffsets = self._read_bytes(aux + 'suboffsets.b')
            self.addresses = self._read_bytes(aux + 'addresses.b')
            self.weak_ascents = self._read_bytes(aux + 'weak_ascents.b')
            self.strict_ascents = self._read_bytes(aux + 'strict_ascents.b')
            self.descents = self._read_bytes(aux + 'descents.b')

            self.even_intervals = []
            self.even_invleft = []
            self.even_invright = []
            self.odd_intervals = []
            self.odd_invleft = []
            self.odd_invright = []
            for i in range(self.qpmodule.rank):
                self.even_intervals.append(self._read_bytes(aux + 'even_intervals.%s.b' % i))
                self.even_invleft.append(self._read_bytes(aux + 'even_invleft.%s.b' % i))
                self.even_invright.append(self._read_bytes(aux + 'even_invright.%s.b' % i))
                self.odd_intervals.append(self._read_bytes(aux + 'odd_intervals.%s.b' % i))
                self.odd_invleft.append(self._read_bytes(aux + 'odd_invleft.%s.b' % i))
                self.odd_invright.append(self._read_bytes(aux + 'odd_invright.%s.b' % i))

        if self.is_cbasis_computed:
            bytefile = self.get_directory() + 'wgraph.cbasis.b'
            with open(bytefile, 'rb') as file:
                self.frame = bytearray(file.read())

        if self.is_wgraph_computed:
            bytefile = self.get_directory() + 'wgraph.edges.b'
            with open(bytefile, 'rb') as file:
                self.wgraph = bytearray(file.read())

            bytefile = self.get_directory() + 'wgraph.addresses.b'
            with open(bytefile, 'rb') as file:
                self.wgraph_addresses = bytearray(file.read())

            self.cells = self._compute_cells(self.get_wgraph_edges)
            self.molecules = self._compute_cells(self.get_simple_edges)

    @classmethod
    def _read_bytes(cls, filename):
        with open(filename, 'rb') as file:
            return bytearray(file.read())

    @classmethod
    def _write_bytes(cls, filename, array):
        with open(filename, 'wb') as file:
            file.write(array)

    def write(self):
        self.qpmodule.write()
        directory = self.get_directory()
        Path(directory).mkdir(parents=True, exist_ok=True)

        metafile = directory + 'wgraph.meta'
        with open(metafile, 'w') as file:
            file.write(json.dumps(self.metadata()))

        if not self.is_setup:
            return

        aux = directory + 'aux/'
        Path(aux).mkdir(parents=True, exist_ok=True)

        self._write_bytes(aux + 'suboffsets.b', self.suboffsets)
        self._write_bytes(aux + 'addresses.b', self.addresses)
        self._write_bytes(aux + 'weak_ascents.b', self.weak_ascents)
        self._write_bytes(aux + 'strict_ascents.b', self.strict_ascents)
        self._write_bytes(aux + 'descents.b', self.descents)

        for i in range(self.qpmodule.rank):
            self._write_bytes(aux + 'even_intervals.%s.b' % i, self.even_intervals[i])
            self._write_bytes(aux + 'even_invleft.%s.b' % i, self.even_invleft[i])
            self._write_bytes(aux + 'even_invright.%s.b' % i, self.even_invright[i])
            self._write_bytes(aux + 'odd_intervals.%s.b' % i, self.odd_intervals[i])
            self._write_bytes(aux + 'odd_invleft.%s.b' % i, self.odd_invleft[i])
            self._write_bytes(aux + 'odd_invright.%s.b' % i, self.odd_invright[i])

        if not self.is_cbasis_computed:
            return

        self._write_bytes(directory + 'wgraph.cbasis.b', self.frame)

        if not self.is_wgraph_computed:
            return

        self._write_bytes(directory + 'wgraph.edges.b', self.wgraph)
        self._write_bytes(directory + 'wgraph.addresses.b', self.wgraph_addresses)

    def metadata(self):
        assert self.is_setup
        return {
            'sgn': self.sgn,
            'nbytes': self.nbytes,
            'heights': self.heights,
            'cumheights': self.cumheights,
            'size': self.size,
            'hbytes': self.hbytes,
            'offsets': self.offsets,
            'abytes': self.abytes,
            'sbytes': self.sbytes,
            'gbytes': self.gbytes,
            'wbytes': self.wbytes,
            'is_setup': self.is_setup,
            'is_cbasis_computed': self.is_cbasis_computed,
            'is_wgraph_computed': self.is_wgraph_computed,
        }

    def _ascents(self, i):
        return self._weak_ascents(i) | self._strict_ascents(i)

    def _weak_ascents(self, i):
        return self._int(self.weak_ascents[i * self.gbytes:(i + 1) * self.gbytes])

    def _strict_ascents(self, i):
        return self._int(self.strict_ascents[i * self.gbytes:(i + 1) * self.gbytes])

    def _descents(self, i):
        return self._int(self.descents[i * self.gbytes:(i + 1) * self.gbytes])

    @classmethod
    def _bitmap(cls, s):
        ans = 0
        for i in sorted(s):
            ans += 2**i
        return ans

    @classmethod
    def _bytes(cls, n):
        assert n >= 0
        ans = 1
        n = n // 256
        while n > 0:
            ans, n = ans + 1, n // 256
        return ans

    def _compute_intervals(self, parity):
        def is_descent(n, i):
            if self.qpmodule.is_strict_descent(n, i):
                return True
            if not self.sgn:
                return self.qpmodule.is_weak_descent(n, i)
            else:
                return self.qpmodule.is_weak_ascent(n, i)

        intervals = [[] for _ in range(self.qpmodule.rank)]
        for n in self.qpmodule:
            for i in range(self.qpmodule.rank):
                if is_descent(n, i) and self._height(n) % 2 == parity:
                    intervals[i].append(n)

        invert_left = [self.qpmodule.size * [0] for _ in range(self.qpmodule.rank)]
        catchup = [[] for _ in range(self.qpmodule.rank)]
        locations = self.qpmodule.rank * [0]
        for n in self.qpmodule:
            for i in range(self.qpmodule.rank):
                catchup[i].append(n)
                if is_descent(n, i) and self._height(n) % 2 == parity:
                    for m in catchup[i]:
                        invert_left[i][m] = locations[i]
                    catchup[i] = []
                    locations[i] += 1

        invert_right = [self.qpmodule.size * [0] for _ in range(self.qpmodule.rank)]
        catchup = [[] for _ in range(self.qpmodule.rank)]
        locations = [len(i) for i in intervals]
        for n in range(len(self.qpmodule) - 1, -1, -1):
            for i in range(self.qpmodule.rank):
                catchup[i].append(n)
                if is_descent(n, i) and self._height(n) % 2 == parity:
                    for m in catchup[i]:
                        invert_right[i][m] = locations[i]
                    catchup[i] = []
                    locations[i] -= 1

        intervals = [self._tobytearray(i, self.sbytes) for i in intervals]
        invert_left = [self._tobytearray(i, self.sbytes) for i in invert_left]
        invert_right = [self._tobytearray(i, self.sbytes) for i in invert_right]
        return intervals, invert_left, invert_right

    @classmethod
    def _tobytearray(cls, iterable, nbytes, signed=False):
        array = bytearray(len(iterable) * nbytes)
        start = 0
        for a in iterable:
            array[start:start + nbytes] = a.to_bytes(nbytes, byteorder='big', signed=signed)
            start += nbytes
        return array

    def _even_invleft(self, n, s):
        start = n * self.sbytes
        return self._int(self.even_invleft[s][start:start + self.sbytes])

    def _even_invright(self, n, s):
        start = n * self.sbytes
        return self._int(self.even_invright[s][start:start + self.sbytes])

    def _odd_invleft(self, n, s):
        start = n * self.sbytes
        return self._int(self.odd_invleft[s][start:start + self.sbytes])

    def _odd_invright(self, n, s):
        start = n * self.sbytes
        return self._int(self.odd_invright[s][start:start + self.sbytes])

    def permutation(self, n):
        return self.qpmodule.permutation(n, True)

    def _height(self, n):
        return self.qpmodule.height(n)

    def _space(self, hi, hj):
        return (1 + (hj - hi - 1) // 2) * self.nbytes

    @classmethod
    def _int(cls, step, signed=False):
        return int.from_bytes(step, byteorder='big', signed=signed)

    def __repr__(self):
        a = str((self.qpmodule.family, self.qpmodule.rank, self.qpmodule.layer, self.sgn))
        b = len(self.qpmodule)
        s = ['QPWGraph for %s (%s elements)' % (a, b)]
        return '\n'.join(s)

    def _address_cbasis(self, i, j, hi, hj):
        start_i = i * self.hbytes
        start_j = j * self.abytes
        space = self._space(hi, hj)
        start = self._int(self.addresses[start_j:start_j + self.abytes]) + \
            self.offsets[hj][hi] + \
            self._int(self.suboffsets[start_i:start_i + self.hbytes]) * space
        return start, space

    def address_cbasis(self, i, j):
        hi = self._height(i)
        hj = self._height(j)
        assert hi < hj
        return self._address_cbasis(i, j, hi, hj)

    def get_cbasis_leading(self, i, j):
        if self._height(i) >= self._height(j):
            return 0
        if (self._height(j) - self._height(i)) % 2 == 0:
            return 0
        start, space = self.address_cbasis(i, j)
        stop = start + space
        return self._int(self.frame[stop - self.nbytes:stop], signed=True)

    def get_cbasis(self, i, j):
        if i == j:
            return (1).to_bytes(self.nbytes, byteorder='big', signed=True)

        hi = self._height(i)
        hj = self._height(j)

        if hi >= hj:
            return (0).to_bytes(self.nbytes, byteorder='big', signed=True)

        start, space = self._address_cbasis(i, j, hi, hj)
        return self.frame[start:start + space]

    def get_cbasis_polynomial(self, i, j):
        if i == j:
            return polynomials.q(0)
        c = self.get_cbasis(i, j)
        ans = 0 * polynomials.q(0)
        start = 0
        for e in range(1 + (j - i - 1) // 2):
            ans += polynomials.q(e) * self._int(c[start:start + self.nbytes], signed=True)
            start += self.nbytes
        return ans

    def set_cbasis(self, i, j, v, set_bytes=True):
        start, size = self.address_cbasis(i, j)
        if set_bytes:
            self.frame[start:start + size] = v
        else:
            self.frame[start:start + size] = v.to_bytes(size, byteorder='big', signed=True)

    def _safe_add(self, start, f, shift=0, mu=1):
        if mu == 0:
            return
        astart = start + shift * self.nbytes
        fstart = 0
        while fstart < len(f):
            a = self._int(self.frame[astart:astart + self.nbytes], signed=True)
            b = self._int(f[fstart:fstart + self.nbytes], signed=True)
            v = (a + mu * b).to_bytes(self.nbytes, byteorder='big', signed=True)
            self.frame[astart:astart + self.nbytes] = v
            astart += self.nbytes
            fstart += self.nbytes

    def _safe_set(self, start, i, j):
        f = self.get_cbasis(i, j)
        self.frame[start:start + len(f)] = f

    @classmethod
    def _firstbit(cls, x):
        if x == 0:
            return
        ans = 0
        while x % 2 == 0:
            x = x >> 1
            ans += 1
        return ans

    def _get_interval(self, i, sj, s, hj):
        if hj % 2 == 0:
            a, b = self._even_invleft(i, s), self._even_invright(sj - 1, s)
            interval = self.even_intervals[s]
        else:
            a, b = self._odd_invleft(i, s), self._odd_invright(sj - 1, s)
            interval = self.odd_intervals[s]
        for x in range(a, b):
            yield self._int(interval[x * self.sbytes:(x + 1) * self.sbytes])

    def _slowcompute(self, verbose=True):
        assert self.is_setup
        assert not self.is_cbasis_computed

        t0 = time.time()
        self.frame = bytearray(self.size)

        if verbose:
            print('Computing canonical basis:')
        progress = 0

        for j in self.qpmodule:
            hj = self._height(j)
            wdes = set(self.qpmodule.weak_descents(j))
            sdes = set(self.qpmodule.strict_descents(j))
            des = wdes | sdes

            for i in range(j - 1, -1, -1):
                hi = self._height(i)
                if hi == hj:
                    continue

                start, _ = self.address_cbasis(i, j)

                if des & set(self.qpmodule.weak_ascents(i)):
                    continue

                if verbose:
                    newprogress = int(100 * start / self.size)
                    if newprogress > progress:
                        print('*', newprogress, 'percent done (%s seconds elapsed)' % str(int(1000 * (time.time() - t0)) / 1000.0))
                        progress = newprogress

                t = set(self.qpmodule.strict_ascents(i)) & des
                if t:
                    x = self.qpmodule.operate(i, next(iter(t)))
                    self._safe_add(start, self.get_cbasis(x, j))
                    continue

                s = next(iter(sdes))
                si = self.qpmodule.operate(i, s)
                sj = self.qpmodule.operate(j, s)

                self._safe_add(start, self.get_cbasis(si, sj))
                self._safe_add(start, self.get_cbasis(i, sj), shift=1)
                for x in range(i, sj):
                    if s in self.qpmodule.weak_ascents(x) or s in self.qpmodule.strict_ascents(x):
                        continue
                    hx = self._height(x)
                    if (hj - hx) % 2 != 0:
                        continue
                    mu = self.get_cbasis_leading(x, sj)
                    self._safe_add(start, self.get_cbasis(i, x), (hj - hx) // 2, -mu)

        if verbose:
            print('Done computing (%s seconds elapsed)' % str(int(1000 * (time.time() - t0)) / 1000.0))
        self.is_cbasis_computed = True

    def compute_cbasis(self, verbose=True):
        assert self.is_setup
        assert not self.is_cbasis_computed

        t0 = time.time()
        self.frame = bytearray(self.size)

        if verbose:
            print('Computing canonical basis:')
        progress = 0

        start = 0
        nextstart = 0
        for j in self.qpmodule:
            hj = self._height(j)
            s = self.qpmodule.get_strict_descent(j)
            for i in range(self.cumheights[hj] - 1, -1, -1):
                hi = self._height(i)
                start = nextstart
                nextstart += self._space(hi, hj)

                if verbose:
                    newprogress = int(100 * start / self.size)
                    if newprogress > progress:
                        print('*', newprogress, 'percent done (%s seconds elapsed)' % str(int(1000 * (time.time() - t0)) / 1000.0))
                        progress = newprogress

                if self._descents(j) & self._weak_ascents(i) != 0:
                    continue

                x = self._firstbit(self._descents(j) & self._strict_ascents(i))
                if x is not None:
                    x = self.qpmodule.operate(i, x)
                    self._safe_set(start, x, j)
                    continue

                si = self.qpmodule.operate(i, s)
                sj = self.qpmodule.operate(j, s)

                self._safe_set(start, si, sj)
                self._safe_add(start, self.get_cbasis(i, sj), shift=1)
                for x in self._get_interval(i, sj, s, hj):
                    hx = self._height(x)
                    mu = self.get_cbasis_leading(x, sj)
                    self._safe_add(start, self.get_cbasis(i, x), (hj - hx) // 2, -mu)

        if verbose:
            print('Done computing (%s seconds elapsed)' % str(int(1000 * (time.time() - t0)) / 1000.0))
        self.is_cbasis_computed = True


class QPModuleElement:

    def __init__(self, qpmodule, n):
        self.qpmodule = qpmodule
        self.n = n

    def __eq__(self, other):
        assert type(self) == type(other)
        assert self.qpmodule == other.qpmodule
        return self.n == other.n

    def __hash__(self):
        return hash((self.n, self.qpmodule))

    def __repr__(self):
        return str(self.qpmodule.permutation(self.n))

    def __mul__(self, i):
        result = self.qpmodule.operate(self.n, i)
        if result == self.n:
            if i in self.qpmodule.weak_ascents(self.n):
                return Vector({self: polynomials.q(1)})
            else:
                return Vector({self: -polynomials.q(-1)})
        new = QPModuleElement(self.qpmodule, result)
        if result < self.n:
            return Vector({self: polynomials.q(1) - polynomials.q(-1), new: 1})
        else:
            return Vector({new: 1})


class QPModule:

    DIRECTORY = '/Users/emarberg/Desktop/examples/models/'

    HECKE_A = 'RIGHT_HECKE_A'
    HECKE_BC = 'RIGHT_HECKE_BC'
    HECKE_D = 'RIGHT_HECKE_D'

    TWO_SIDED_HECKE_A = 'TWO_SIDED_HECKE_A'
    TWO_SIDED_HECKE_BC = 'TWO_SIDED_HECKE_BC'
    TWO_SIDED_HECKE_D = 'TWO_SIDED_HECKE_D'

    GELFAND_A = 'GELFAND_A'
    GELFAND_BC = 'GELFAND_BC'
    GELFAND_D = 'GELFAND_D'

    def __init__(self, family, rank, layer, size, height_bytes, frame=None, printer=None):
        self.family = family
        self.rank = rank
        self.layer = layer
        self.size = size
        self.stepsize = max(1, math.ceil(math.log2((size + 2) / 8.0)))
        self.height_bytes = height_bytes
        self.frame = frame
        self.printer = printer
        self.permutation_cache = {}

    def expand(self, n):
        if self.printer:
            return self.printer(self.reduced_word(n))
        else:
            word = self.reduced_word(n)
            return len(word), word

    def length(self, n):
        return self.height(n)

    def height(self, n):
        return self[n][0]

    def __eq__(self, other):
        if self.layer != other.layer:
            return False
        if self.family != other.family or self.rank != other.rank or self.size != other.size:
            return False
        if self.stepsize != other.stepsize or self.height_bytes != other.height_bytes:
            return False
        return self.frame == other.frame

    def __hash__(self):
        return hash((self.family, self.rank, self.layer))

    def __len__(self):
        return self.size

    def __repr__(self):
        k = 24
        s = [self.family + str(self.rank) + ' : LAYER ' + str(self.layer)]
        s += ['']
        if self.frame:
            for i in range(min(self.size, k)):
                ht, w = self.expand(i)
                des = str(list(self.weak_descents(i)))
                asc = str(list(self.weak_ascents(i)))
                s += ['[ height ' + str(ht) + ' : ' + str(w) + ' : weak descents ' + des + ' : weak ascents ' + asc + ' ]']
            if self.size > k:
                s += ['', '(... and %s more)' % (self.size - k)]
        return '\n' + '\n'.join(s) + '\n'

    def reduced_word(self, n):
        for i in self.strict_descents(n):
            return self.reduced_word(self.operate(n, i)) + (i,)
        return ()

    def __iter__(self):
        return iter(range(self.size))

    def __getitem__(self, n):
        start = n * (self.height_bytes + self.rank * self.stepsize)
        ans = [self._int(self.frame[start:start + self.height_bytes])]
        start += self.height_bytes
        return ans + [
            self._int(self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize])
            for i in range(self.rank)
        ]

    @classmethod
    def _is_weak_ascent(cls, frame):
        return all(b == 0xFF for b in frame)

    @classmethod
    def _is_weak_descent(cls, frame):
        return all(b == 0xFF for b in frame[:-1]) and frame[-1] == 0xFE

    def is_strict_descent(self, n, i):
        return i in set(self.strict_descents(n))

    def is_strict_ascent(self, n, i):
        return i in set(self.strict_ascents(n))

    def is_weak_descent(self, n, i):
        return i in set(self.weak_descents(n))

    def is_weak_ascent(self, n, i):
        return i in set(self.weak_ascents(n))

    @classmethod
    def _int(cls, step, signed=False):
        return int.from_bytes(step, byteorder='big', signed=signed)

    def element(self, n):
        assert 0 <= n < self.size
        return Vector({QPModuleElement(self, n): 1})

    def operate(self, n, *args):
        for i in args:
            assert i in range(self.rank)
            start = n * (self.height_bytes + self.rank * self.stepsize) + self.height_bytes
            step = self.frame[start + i * self.stepsize:start + (i + 1) * self.stepsize]
            if self._is_weak_descent(step) or self._is_weak_ascent(step):
                pass
            else:
                n = self._int(step)
        return n

    def _steps(self, n):
        start = n * (self.height_bytes + self.rank * self.stepsize) + self.height_bytes
        for _ in range(self.rank):
            yield self.frame[start:start + self.stepsize]
            start += self.stepsize

    def strict_ascents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_descent(step):
                if self._int(step) > n:
                    yield i

    def get_strict_descent(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_descent(step):
                if self._int(step) < n:
                    return i

    def strict_descents(self, n):
        for i, step in enumerate(self._steps(n)):
            if not self._is_weak_ascent(step) and not self._is_weak_descent(step):
                if self._int(step) < n:
                    yield i

    def weak_ascents(self, n):
        for i, step in enumerate(self._steps(n)):
            if self._is_weak_ascent(step):
                yield i

    def weak_descents(self, n):
        for i, step in enumerate(self._steps(n)):
            if self._is_weak_descent(step):
                yield i

    def get_directory(self):
        return self.directory(self.family, self.rank, self.layer)

    @classmethod
    def directory(cls, family, rank, layer=None):
        file = cls.DIRECTORY + family + str(rank)
        if layer is not None:
            file += '_BLOCK' + str(layer)
        return file + '/'

    def metadata(self):
        return {
            'family': self.family,
            'rank': self.rank,
            'layer': self.layer,
            'size': self.size,
            'stepsize': self.stepsize,
            'height_bytes': self.height_bytes,
        }

    def write(self):
        directory = self.get_directory()
        Path(directory).mkdir(parents=True, exist_ok=True)
        bytefile = directory + 'module.b'
        with open(bytefile, 'wb') as file:
            file.write(self.frame)
        metafile = directory + 'module.meta'
        with open(metafile, 'w') as file:
            file.write(json.dumps(self.metadata()))

    @classmethod
    def read_two_sided_hecke_a(cls, n, layer=None):
        return cls.read(cls.directory(cls.TWO_SIDED_HECKE_A, 2 * n))

    @classmethod
    def read_two_sided_hecke_bc(cls, n, layer=None):
        return cls.read(cls.directory(cls.TWO_SIDED_HECKE_BC, 2 * n))

    @classmethod
    def read_two_sided_hecke_d(cls, n, layer=None):
        return cls.read(cls.directory(cls.TWO_SIDED_HECKE_D, 2 * n))

    @classmethod
    def read_hecke_a(cls, n, layer=None):
        return cls.read(cls.directory(cls.HECKE_A, n))

    @classmethod
    def read_hecke_bc(cls, n, layer=None):
        return cls.read(cls.directory(cls.HECKE_BC, n))

    @classmethod
    def read_hecke_d(cls, n, layer=None):
        return cls.read(cls.directory(cls.HECKE_D, n))

    @classmethod
    def read_gelfand_a(cls, n, k):
        return cls.read(cls.directory(cls.GELFAND_A, n, k))

    @classmethod
    def read_gelfand_bc(cls, n, k):
        return cls.read(cls.directory(cls.GELFAND_BC, n, k))

    @classmethod
    def read_gelfand_d(cls, n, k):
        return cls.read(cls.directory(cls.GELFAND_D, n, k))

    @classmethod
    def read(cls, directory):
        metafile = directory + 'module.meta'
        with open(metafile, 'r') as file:
            dictionary = json.loads(file.read())

        family = dictionary['family']
        rank = dictionary['rank']
        layer = dictionary['layer']
        size = dictionary['size']
        stepsize = dictionary['stepsize']
        height_bytes = int(dictionary['height_bytes'])
        module = cls(family, rank, layer, size, height_bytes)
        assert module.stepsize == stepsize

        bytefile = directory + 'module.b'
        with open(bytefile, 'rb') as file:
            module.frame = bytearray(file.read())
        return module

    @classmethod
    def create(cls, rank, size, stepsize, height_bytes, minima, operate, verbose=True):
        assert size * (rank * stepsize + height_bytes) < 0.8e+10

        t0 = time.time()
        if verbose:
            print()
            print('RANK', rank, 'SIZE', size)
            print()
            print('* creating bytearray of size %s bytes' % (size * rank * stepsize))

        frame = bytearray(size * (rank * stepsize + height_bytes))
        level = [(wht, time.time(), w) for (wht, w) in minima]
        origins = {w: [] for (_, w) in minima}
        position, start = 0, 0

        while level:
            wht, _, w = heappop(level)
            assert position < size
            frame[start:start + height_bytes] = wht.to_bytes(height_bytes, byteorder='big')
            start += height_bytes
            descents = {}
            for (i, n, o) in origins[w]:
                frame[o:o + stepsize] = position.to_bytes(stepsize, byteorder='big')
                descents[i] = n
            for i in range(rank):
                if i in descents:
                    frame[start:start + stepsize] = descents[i].to_bytes(stepsize, byteorder='big')
                else:
                    ht, y = operate(wht, w, i)
                    if y == (w, True):
                        # i is weak ascent
                        frame[start:start + stepsize] = bytes(stepsize * [0xFF])
                    elif y == (w, False):
                        # i is weak descent
                        frame[start:start + stepsize] = bytes((stepsize - 1) * [0xFF] + [0xFE])
                    else:
                        if y not in origins:
                            heappush(level, (ht, time.time(), y))
                        origins[y] = origins.get(y, []) + [(i, position, start)]
                start += stepsize
            position += 1
            del origins[w]

        t1 = time.time()
        if verbose:
            print('* finished, time elapsed: %s milliseconds' % int(1000 * (t1 - t0)))
            print()
        return frame

    @classmethod
    def create_hecke_classical(cls, family, n, size, height_bytes, coxeter_class):
        assert n >= 0
        rank = 2 * n if family in [cls.TWO_SIDED_HECKE_A, cls.TWO_SIDED_HECKE_BC, cls.TWO_SIDED_HECKE_D] else n
        module = cls(family, rank, None, size, height_bytes)
        stepsize = module.stepsize

        def multiply(ht, w, j):
            i = j % n
            if coxeter_class == Permutation:
                s = coxeter_class.s_i(i + 1)
            else:
                s = coxeter_class.s_i(i, w.rank)
            ws = w * s if i == j else s * w
            ht = ws.length()
            return ht, ws

        def printer(word):
            ht, v = 0, coxeter_class.identity(n)
            for i in word:
                ht, v = multiply(ht, v, i)
            return ht, v

        minima = [(0, coxeter_class.identity(n))]
        module.frame = cls.create(rank, size, stepsize, height_bytes, minima, multiply)
        module.printer = printer
        return module

    @classmethod
    def create_hecke_a(cls, n, layer=None):
        assert n >= 0
        size = math.factorial(n + 1)
        return cls.create_hecke_classical(cls.HECKE_A, n, size, 1, Permutation)

    @classmethod
    def create_hecke_bc(cls, n, layer=None):
        assert n >= 2
        size = math.factorial(n) * 2**n
        return cls.create_hecke_classical(cls.HECKE_BC, n, size, 1, SignedPermutation)

    @classmethod
    def create_hecke_d(cls, n, layer=None):
        assert n >= 2
        size = math.factorial(n) * 2**(n - 1)
        return cls.create_hecke_classical(cls.HECKE_D, n, size, 1, EvenSignedPermutation)

    @classmethod
    def create_two_sided_hecke_a(cls, n, layer=None):
        assert n >= 0
        size = math.factorial(n + 1)
        return cls.create_hecke_classical(cls.TWO_SIDED_HECKE_A, n, size, 1, Permutation)

    @classmethod
    def create_two_sided_hecke_bc(cls, n, layer=None):
        assert n >= 2
        size = math.factorial(n) * 2**n
        return cls.create_hecke_classical(cls.TWO_SIDED_HECKE_BC, n, size, 1, SignedPermutation)

    @classmethod
    def create_two_sided_hecke_d(cls, n, layer=None):
        assert n >= 2
        size = math.factorial(n) * 2**(n - 1)
        return cls.create_hecke_classical(cls.TWO_SIDED_HECKE_D, n, size, 1, EvenSignedPermutation)

    @classmethod
    def create_gelfand_a(cls, n, k):
        assert 0 <= 2 * k <= n + 1

        size = math.factorial(n + 1) // math.factorial(k) // 2**k // math.factorial(n + 1 - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_A, n, k, size, height_bytes)
        stepsize = module.stepsize

        minima = [(0, gelfand_a_start(n, k))]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, gelfand_a_conjugate)
        module.printer = gelfand_a_printer(n, k)
        return module

    @classmethod
    def create_gelfand_bc(cls, n, k):
        assert 0 <= 2 * k <= n

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_BC, n, k, size, height_bytes)
        stepsize = module.stepsize

        minima = [(0, gelfand_bc_start(n, k))]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, gelfand_bc_conjugate)
        module.printer = gelfand_bc_printer(n, k)
        return module

    @classmethod
    def create_gelfand_d(cls, n, k):
        assert n >= 2
        assert 0 <= 2 * k <= n

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        assert size % 2 == 0
        size = size // 2

        height_bytes = 1
        module = cls(cls.GELFAND_D, n, k, size, height_bytes)
        stepsize = module.stepsize

        minima = [(0, gelfand_d_start(n, k))]
        module.frame = cls.create(n, size, stepsize, height_bytes, minima, gelfand_d_conjugate)
        module.printer = gelfand_d_printer(n, k)
        return module

    @classmethod
    def generator(cls, family, n, k=None, sgn=None):
        if family == cls.HECKE_A or family == cls.TWO_SIDED_HECKE_A:
            return Permutation()
        elif family == cls.HECKE_BC:
            return SignedPermutation.identity(n)
        elif family == cls.HECKE_D:
            return EvenSignedPermutation.identity(n)
        elif family == cls.TWO_SIDED_HECKE_BC:
            return SignedPermutation.identity(n // 2)
        elif family == cls.TWO_SIDED_HECKE_D:
            return EvenSignedPermutation.identity(n // 2)
        elif family == cls.GELFAND_A:
            return Permutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [n + 2 + i for i in range(n + 1 - 2 * k)] +
                [2 * k + 1 + i for i in range(n + 1 - 2 * k)]
            )) if not sgn else Permutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [i for i in range(2 * n + 2 - 2 * k, 2 * k, -1)]
            ))
        elif family == cls.GELFAND_BC:
            return SignedPermutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [n + 1 + i for i in range(n - 2 * k)] +
                [2 * k + 1 + i for i in range(n - 2 * k)]
            )) if not sgn else SignedPermutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [i for i in range(2 * n - 2 * k, 2 * k, -1)]
            ))
        elif family == cls.GELFAND_D:
            return EvenSignedPermutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [n + 1 + i for i in range(n - 2 * k)] +
                [2 * k + 1 + i for i in range(n - 2 * k)]
            )) if not sgn else EvenSignedPermutation(*(
                [1 + i + (-1)**i for i in range(2 * k)] +
                [i for i in range(2 * n - 2 * k, 2 * k, -1)]
            ))
        else:
            raise Exception

    def permutation(self, n, sgn=None):
        if n not in self.permutation_cache:
            w = self.generator(self.family, self.rank, self.layer, sgn)
            for i in self.reduced_word(n):
                if self.family == self.HECKE_A:
                    w *= Permutation.s_i(i + 1)
                elif self.family == self.HECKE_BC:
                    w *= SignedPermutation.s_i(i, w.rank)
                elif self.family == self.HECKE_D:
                    w *= EvenSignedPermutation.s_i(i, w.rank)
                elif self.family == self.TWO_SIDED_HECKE_A:
                    m = self.rank // 2
                    s = Permutation.s_i((i % m) + 1)
                    w = w * s if i % m == i else s * w
                elif self.family == self.TWO_SIDED_HECKE_BC:
                    m = self.rank // 2
                    s = SignedPermutation.s_i(i % m, w.rank)
                    w = w * s if i % m == i else s * w
                elif self.family == self.TWO_SIDED_HECKE_D:
                    m = self.rank // 2
                    s = EvenSignedPermutation.s_i(i % m, w.rank)
                    w = w * s if i % m == i else s * w
                elif self.family == self.GELFAND_A:
                    s = Permutation.s_i(i + 1)
                    w = s * w * s
                elif self.family == self.GELFAND_BC:
                    s = SignedPermutation.s_i(i, w.rank)
                    w = s * w * s
                elif self.family == self.GELFAND_D:
                    s = EvenSignedPermutation.s_i(i, w.rank)
                    w = s * w * s
                else:
                    raise Exception
            if self.family == self.GELFAND_D and sgn and self.layer % 2 != self.rank % 2:
                w = w.star()
            if self.family == self.HECKE_A:
                self.permutation_cache[n] = tuple(w(i) for i in range(1, self.rank + 2))
            elif self.family == self.TWO_SIDED_HECKE_A:
                self.permutation_cache[n] = tuple(w(i) for i in range(1, self.rank // 2 + 2))
            else:
                self.permutation_cache[n] = tuple(w.oneline)
        return self.permutation_cache[n]

    @classmethod
    def slow_create_gelfand_a(cls, n, k):
        assert 0 <= 2 * k <= n + 1

        m = 2 * n + 2 - 2 * k
        a = [Permutation.s_i(i) for i in range(n + 2, m)]
        s = {i: Permutation.s_i(i + 1) for i in range(n + 1)}
        w = cls.generator(cls.GELFAND_A, n, k)

        size = math.factorial(n + 1) // math.factorial(k) // 2**k // math.factorial(n + 1 - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_A, n, k, size, height_bytes)
        return cls.create_gelfand_classical(module, a, s, w)

    @classmethod
    def slow_create_gelfand_bc(cls, n, k):
        assert 0 <= 2 * k <= n

        m = 2 * n - 2 * k
        a = [SignedPermutation.s_i(i, m) for i in range(n + 1, m)]
        s = {i: SignedPermutation.s_i(i, m) for i in range(n)}
        w = cls.generator(cls.GELFAND_BC, n, k)

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        height_bytes = 1
        module = cls(cls.GELFAND_BC, n, k, size, height_bytes)
        return cls.create_gelfand_classical(module, a, s, w)

    @classmethod
    def slow_create_gelfand_d(cls, n, k):
        assert n >= 2
        assert 0 <= 2 * k <= n

        m = 2 * n - 2 * k
        a = [EvenSignedPermutation.s_i(i, m) for i in range(n + 1, m)]
        s = {i: EvenSignedPermutation.s_i(i, m) for i in range(n)}
        w = cls.generator(cls.GELFAND_D, n, k)

        size = math.factorial(n) * 2**(n - 2 * k) // math.factorial(k) // math.factorial(n - 2 * k)
        assert size % 2 == 0
        size = size // 2

        height_bytes = 1
        module = cls(cls.GELFAND_D, n, k, size, height_bytes)
        return cls.create_gelfand_classical(module, a, s, w)

    @classmethod
    def create_gelfand_classical(cls, module, a, s, w):
        def conjugate(ht, u, i):
            if u * s[i] * u in a:
                return ht, (u, True)
            v = s[i] * u * s[i]
            if v == u:
                return ht, (u, False)
            return ht + 1, v

        def printer(word):
            ht, v = 0, w
            for i in word:
                ht, v = conjugate(ht, v, i)
            return ht, v

        minima = [(0, w)]
        module.frame = cls.create(
            module.rank,
            module.size,
            module.stepsize,
            module.height_bytes,
            minima,
            conjugate
        )
        module.printer = printer
        return module


def gelfand_a_conjugate(ht, w, i):
    i += 1
    if w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if w(i) == i + 1 and w(i + 1) == i:
        return ht, (w, False)
    s = Permutation.s_i(i)
    return ht + 1, s * w * s


def gelfand_a_start(n, k):
    return Permutation(*[1 + i + (-1)**i for i in range(2 * k)])


def gelfand_a_printer(n, k):
    w = gelfand_a_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_a_conjugate(ht, v, i)
        return ht, v

    return printer


def gelfand_bc_conjugate(ht, w, i):
    s = SignedPermutation.s_i(i, w.rank)
    if i == 0 and w(1) in [-1, 1]:
        return ht + 1, w * s
    if i > 0 and abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
        return ht, (w, False)
    if i > 0 and w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if i > 0 and w(i) == -i and w(i + 1) == -i - 1:
        return ht, (w, True)
    return ht + 1, s * w * s


def gelfand_bc_start(n, k):
    return SignedPermutation(*(
        [1 + i + (-1)**i for i in range(2 * k)] +
        [i for i in range(2 * k + 1, n + 1)]
    ))


def gelfand_bc_printer(n, k):
    w = gelfand_bc_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_bc_conjugate(ht, v, i)
        return ht, v

    return printer


def gelfand_d_conjugate(ht, w, i):
    if i == 0:
        s = SignedPermutation.ds_i(-1, w.rank)
        t = SignedPermutation.ds_i(1, w.rank)
        if abs(w(1)) != 1 and abs(w(2)) != 2:
            if w * s == s * w:
                return ht, (w, False)
            return ht + 1, s * w * s
        if (w(1) == 1 and w(2) == 2) or (w(1) == -1 and w(2) == -2):
            return ht + 1, s * w * t
        if (w(1) == 1 and w(2) == -2) or (w(1) == -1 and w(2) == 2):
            return ht, (w, True)
        if (abs(w(1)) == 1 and abs(w(2)) != 2) or (abs(w(1)) != 1 and abs(w(2)) == 2):
            return ht + 1, (s * w * s).dstar()
        raise Exception

    if abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
        return ht, (w, False)
    if w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if w(i) == -i and w(i + 1) == -i - 1:
        return ht, (w, True)

    s = SignedPermutation.ds_i(i, w.rank)
    return ht + 1, s * w * s


def gelfand_d_start(n, k):
    return SignedPermutation(*(
        [1 + i + (-1)**i for i in range(2 * k)] +
        [i for i in range(2 * k + 1, n + 1)]
    ))


def gelfand_d_printer(n, k):
    w = gelfand_d_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_d_conjugate(ht, v, i)
        return ht, v

    return printer
