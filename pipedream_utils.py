from permutations import Permutation
from tableaux import Tableau
from tests.test_pipedreams import sigma_fpf as _sigma_fpf, sigma_k
from words import Word
import words as WORDS


def pipedreams(w):
	# w should be a Permutation
	# example: 
	#   w = Permutation.from_word(3, 2, 3, 4)
	#   pipedreams(w)
	return list(w.get_pipe_dreams())


def flagged_tab(k, mu):
	# example: flagged_tab(k=2, mu=(5, 4, 4, 2))
	return list(Tableau.k_flagged(k, mu))


def rpp(k, mu):
	# TODO
	pass


def biword_star(D, n=None):
	n = max([0] + [i + j for (i, j) in D]) if n is None else n
	cols = sorted([(i, n - (i + j - 1)) for (i, j) in D])
	row_one = [c[0] for c in cols]
	row_two = [c[1] for c in cols]
	return row_one, row_two


def convert_biword_to_list_of_words(biword):
	row_one, row_two = biword
	words = []
	for index in range(len(row_one)):
		i = row_one[index]
		a = row_two[index]
		while i - 1 >= len(words):
			words += [[]]
		words[i - 1] += [a]
	return [Word(*w) for w in words]


def eg_insert(biword):
	words = convert_biword_to_list_of_words(biword)
	return WORDS.eg_insert(*words)


def P_EG(biword):
	return eg_insert(biword)[0]


def Q_EG(biword):
	return eg_insert(biword)[1]


def sigma(k, mu):
	return sigma_k(mu, k)


def biword_fpf(D):
	words = D.column_reading_words()
	row_one, row_two = [], []
	for i in range(len(words)):
		for a in words[i]:
			row_one += [i + 1]
			row_two += [a]
	return row_one, row_two


def fpf_pipedreams(w):
	# w should be a Permutation that is a fixed-point-free involution
	# example:
	#
	# w = Permutation.from_fpf_involution_word(4, 6, 5)
	# fpf_pipedreams(w)
	return list(w.get_fpf_involution_pipe_dreams())


def fpf_insert(biword):
	words = convert_biword_to_list_of_words(biword)
	return WORDS.fpf_insert(*words)


def sigma_fpf(k, mu):
	return _sigma_fpf(mu, k)


def flagged_tab_fpf(k, mu):
	# TODO 
	pass


def rpp_fpf(k, mu):
	return list(Tableau.even_diagonal_unprimed_shifted_rpp(k + 1, mu))


