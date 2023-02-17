from tableaux import Tableau
from permutations import Permutation
from words import involution_insert, fpf_insert


def tableau_from_words(*args):
    mapping = {}
    for i in range(len(args)):
        for j in range(len(args[i])):
            mapping[i + 1, j + 1] = args[i][j]
    return Tableau(mapping)


def shifted_tableau_from_words(*args):
    mapping = {}
    for i in range(len(args)):
        for j in range(len(args[i])):
            mapping[i + 1, i + j + 1] = args[i][j]
    return Tableau(mapping)


def lowest_weight_fpf_insertion(e):
    p, q = fpf_insert(*e)
    mu = q.partition().tuple()
    lowest_tab = Tableau.shifted_lowest_weight(mu)
    lowest_fac = Tableau.inverse_fpf_insertion(p, lowest_tab)
    return shifted_tableau_from_words(*reversed(lowest_fac))


def lowest_weight_inv_insertion(e):
    p, q = involution_insert(*e)
    mu = q.partition().tuple()
    lowest_tab = Tableau.shifted_lowest_weight(mu)
    lowest_fac = Tableau.inverse_involution_insertion(p, lowest_tab)
    return shifted_tableau_from_words(*reversed(lowest_fac))


def highest_weight_fpf_insertion(e):
    m = max(e + (0,)) + 1
    m += int(m % 2 != 0)
    p, q = fpf_insert(*[m - i for i in e])
    mu = q.partition().tuple()
    highest_tab = Tableau.shifted_highest_weight(mu)
    highest_fac = Tableau.inverse_fpf_insertion(p, highest_tab)
    tab = tableau_from_words(*[tuple(reversed(a)) for a in highest_fac])
    tab = Tableau({(j, i): m - abs(tab[i, j]) for (i, j) in tab})
    return tab


def highest_weight_inv_insertion(e):
    m = max(e + (0,)) + 1
    p, q = involution_insert(*[m - i for i in e])
    mu = q.partition().tuple()
    highest_tab = Tableau.shifted_highest_weight(mu)
    highest_fac = Tableau.inverse_involution_insertion(p, highest_tab)
    tab = tableau_from_words(*[tuple(reversed(a)) for a in highest_fac])
    tab = Tableau({(j, i): m - abs(tab[i, j]) for (i, j) in tab})
    return tab
    ## alternate version
    # p, q = involution_insert(*e)
    # mu = q.partition().tuple()
    # highest_tab = Tableau.shifted_highest_weight(mu)
    # highest_fac = Tableau.inverse_involution_insertion(p, highest_tab)
    # tab = tableau_from_words(*[tuple(reversed(a)) for a in highest_fac])
    # return tab


def print_lowest_weight_fpf_insertion(e):
    function = lowest_weight_fpf_insertion
    name = 'lowest-weight fpf-insertion'
    eg_compare = False
    print_insertion(e, function, name, eg_compare)


def print_lowest_weight_inv_insertion(e):
    function = lowest_weight_inv_insertion
    name = 'lowest-weight inv-insertion'
    eg_compare = False
    print_insertion(e, function, name, eg_compare)


def print_highest_weight_fpf_insertion(e):
    function = highest_weight_fpf_insertion
    name = 'highest-weight fpf-insertion'
    eg_compare = True
    print_insertion(e, function, name, eg_compare)


def print_highest_weight_inv_insertion(e):
    function = highest_weight_inv_insertion
    name = 'highest-weight inv-insertion'
    eg_compare = True
    print_insertion(e, function, name, eg_compare)


def print_insertion(e, insertion_function, insertion_name, eg_compare=False):
    if len(e) == 0:
        print('word =', e, 'is empty')
        return

    e_strs = [str(i) for i in e]
    e_strs[-1] += '\t\t(number inserted into previous)'

    highest = [insertion_function(e[:i + 1]) for i in range(len(e))]

    if eg_compare:
        aug = [Tableau()] + highest
        egstep = [aug[j + 1] == aug[j].eg_insert(e[j])[1] for j in range(len(highest))]
    
        b_strs = [' ' if b else '*' for b in egstep]
        b_strs[-1] += '\t\t(steps NOT given by EG insertion)'

        strs = [str(highest[j]) + '\n' + e_strs[j] + '\n\n' + b_strs[j] for j in range(len(e))]
    else:
        strs = [str(highest[j]) + '\n' + e_strs[j] for j in range(len(e))]

    print(combine_tableau_strings(5, *strs))
    print()
    print(insertion_name + ' for word =', e)
    print()


def pad_lines(s):
    s = s.split('\n')
    m = max([len(a) for a in s] + [0])
    s = [a + (m - len(a)) * ' ' for a in s]
    return '\n'.join(s)


def combine_tableau_strings(k, *args):
    jn = k * ' '
    s = [pad_lines(str(a)) for a in args]
    s = [str(a).split('\n') for a in s]
    m = max([len(a) for a in s + ['']])
    s = [(m - len(a)) * [len(a[0]) * ' '] + a for a in s]
    s = [jn.join(a) for a in zip(*s)]
    return '\n'.join(s)


def test_inv(n=5):
    w = Permutation.longest_element(n)
    words = list(w.get_involution_words())
    for e in words:
        lowest = [lowest_weight_inv_insertion(e[:j + 1]) for j in range(len(e))]
        assert all(lowest[j] == involution_insert(*e[:j + 1])[0] for j in range(len(e)))
        highest = [highest_weight_inv_insertion(e[:j + 1]) for j in range(len(e))]
        assert all(t.is_increasing() for t in highest)


def test_fpf(n=6):
    assert n % 2 == 0
    w = Permutation.longest_element(n)
    words = list(w.get_fpf_involution_words())
    for e in words:
        lowest = [lowest_weight_fpf_insertion(e[:j + 1]) for j in range(len(e))]
        assert all(lowest[j] == fpf_insert(*e[:j + 1])[0] for j in range(len(e)))
        highest = [highest_weight_fpf_insertion(e[:j + 1]) for j in range(len(e))]
        assert all(t.is_increasing() for t in highest)
        print_lowest_weight_fpf_insertion(e)
        print_highest_weight_fpf_insertion(e)

