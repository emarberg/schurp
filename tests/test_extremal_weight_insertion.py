from tableaux import Tableau
from permutations import Permutation
from words import involution_insert


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


def lowest_weight_inv_insertion(e):
    p, q = involution_insert(*e)
    mu = q.partition().tuple()
    lowest_tab = Tableau.shifted_lowest_weight(mu)
    lowest_fac = Tableau.inverse_involution_insertion(p, lowest_tab)
    return shifted_tableau_from_words(*reversed(lowest_fac))


def highest_weight_inv_insertion(e):
    p, q = involution_insert(*e)
    mu = q.partition().tuple()
    highest_tab = Tableau.shifted_highest_weight(mu)
    highest_fac = Tableau.inverse_involution_insertion(p, highest_tab)
    return tableau_from_words(*[tuple(reversed(a)) for a in highest_fac])


def combine_tableau_strings(*args):
    s = [str(a).split('\n') for a in args]
    m = max([len(a) for a in s + ['']])
    s = [(m - len(a)) * [len(a[0]) * ' '] + a for a in s]
    s = ['   '.join(a) for a in zip(*s)]
    return '\n'.join(s)


def test_inv(n):
    #for w in sorted(Permutation.involutions(n), key=lambda w: w.involution_length()):
    print('n =', n)
    for w in [Permutation.longest_element(n)]:
        words = list(w.get_involution_words())
        for i, e in enumerate(words):
            if i % (len(words) // 100 + 1) == 0:
                print('  ', 100 * i / len(words), '%') 
            #p, q = involution_insert(*e)
            
            # expected_lowest_tab = lowest_weight_inv_insertion(e)
            # lowest = [lowest_weight_inv_insertion(e[:i + 1]) for i in range(len(e))]

            # expected_highest_tab = highest_weight_inv_insertion(e)
            highest = [highest_weight_inv_insertion(e[:i + 1]) for i in range(len(e))]
            
            try:
                assert all(t.is_decreasing() for t in highest)
                # assert expected_lowest_tab == p
                # assert all(t.is_increasing() for t in lowest)
            except:
                print('word =', e)
                print(combine_tableau_strings(*highest))
                # print(combine_tableau_strings(*lowest))
                print()
                print()
                assert False
