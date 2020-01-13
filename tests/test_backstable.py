from permutations import Permutation
from signed import SignedPermutation
from polynomials import X


def backstable_compatible_sequences(seq, degree_cutoff, i_max=None):
    if len(seq) == 0:
        yield ()
    else:
        seq, z = seq[:-1], seq[-1]
        i_max = max(z, 0) if i_max is None or max(z, 0) < i_max else i_max
        for i in range(i_max, degree_cutoff, -1):
            j_max = (i - 1) if (seq and seq[-1] < z) else i
            for ans in backstable_compatible_sequences(seq, degree_cutoff - (i - 1), j_max):
                yield ans + (i,)


def principal_specialization_summand(a_seq, degree_cutoff):
    ans = 0
    for i_seq in backstable_compatible_sequences(a_seq, degree_cutoff):
        ans += X(0)**(sum(i_seq) - len(i_seq))
    return ans


def principal_specialization(w, degree_cutoff):
    ans = 0
    for a_seq in w.get_reduced_words():
        ans += principal_specialization_summand(a_seq, degree_cutoff)
    return ans


def type_c_principal_specialization(w, degree_cutoff):
    ans = 0
    for a_seq in w.get_signed_reduced_words():
        ans += principal_specialization_summand(a_seq, degree_cutoff)
    return ans * 2**w.ell_zero()


def type_c_comaj_formula_summand(a_seq, degree_cutoff):
    def generate(expon, i):
        if i == 0:
            yield expon
        else:
            for x in [expon, expon - a_seq[i - 1]]:
                while True:
                    x -= i
                    if x < degree_cutoff:
                        break
                    for e in generate(x, i - 1):
                        yield e
    #
    ans = 0
    m = len(a_seq)
    comaj = sum([i + 1 for i in range(m - 1) if a_seq[i] < a_seq[i + 1]]) + sum(a_seq)
    for e in generate(comaj, m):
        ans += X(0)**e
    return ans


def comaj_formula_summand(a_seq, degree_cutoff):
    def generate(expon, i):
        if i == 0:
            yield expon
        else:
            while True:
                expon -= i
                if expon < degree_cutoff:
                    return
                for e in generate(expon, i - 1):
                    yield e
    #
    ans = 0
    m = len(a_seq)
    comaj = sum([i + 1 for i in range(m - 1) if a_seq[i] < a_seq[i + 1]])
    for e in generate(comaj + sum(a_seq), m):
        ans += X(0)**e
    return ans


def comaj_formula(w, degree_cutoff):
    ans = 0
    for a_seq in w.get_reduced_words():
        ans += comaj_formula_summand(a_seq, degree_cutoff)
    return ans


def type_c_comaj_formula(w, degree_cutoff):
    ans = 0
    for a_seq in w.get_reduced_words():
        ans += type_c_comaj_formula_summand(a_seq, degree_cutoff)
    return ans


def factorial(n):
    return 1 if n <= 0 else factorial(n - 1) * n


def test_principal_specialization_summands(n=4, degree_cutoff=0):
    total = factorial(n)
    for i, w in enumerate(Permutation.all(n)):
        print(total - i, 'w =', w)
        for a_seq in w.get_reduced_words():
            r = comaj_formula_summand(a_seq, degree_cutoff)
            s = principal_specialization_summand(a_seq, degree_cutoff)
            print('  ', a_seq, '->', r, '=?=', s)
            if r != s:
                return
        print()
    assert False


def test_principal_specialization(n=4, degree_cutoff=-6):
    total = factorial(n)
    for i, w in enumerate(Permutation.all(n)):
        r = comaj_formula(w, degree_cutoff)
        s = principal_specialization(w, degree_cutoff)
        print(total - i, 'w =', w, '->', s == r)
        assert r == s


def test_type_c_principal_specialization(n=3, degree_cutoff=-10):
    total = factorial(n) * 2**n
    for i, w in enumerate(SignedPermutation.all(n)):
        r = type_c_comaj_formula(w, degree_cutoff)
        s = type_c_principal_specialization(w, degree_cutoff)
        print(total - i, 'w =', w, '=', w.get_reduced_word(), '->', r, '=?=', s)
        assert r == s
