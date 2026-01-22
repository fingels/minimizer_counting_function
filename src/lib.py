import cmath
from src.utils import *
from typing import Callable

# TODO : refactor with integers instead of letters ?
# TODO : option to return list of values instead of just one
# TODO : add tests for VCF
# WARNING : VCF does not have upper and lower bounds implemented

class LexMinimizerCountingFunction(object):
    '''
    Counting function for lexicographical minimizers
    '''
    __slots__ = ('autocorrelation_matrix', 'alphabet', 'minimizer', 'length',
                 'antemer_max_prefix_size', 'postmer_max_size', 'prefix_letters_vectors', 'a_max', 'number_of_greater_letters',
                 'number_of_greater_words','postmer_small_values_alphabet_zero', 'postmer_small_values_alphabet_not_zero',
                 'postmer_prefix_letters_vectors','postmer_a_max')

    def __init__(self, string, alphabet={'A', 'T', 'C', 'G'}, number_of_greater_letters_dic=None):

        self.length: int = len(string)
        self.number_of_greater_words: int = number_of_greater_words(string)
        self.minimizer: str = string + ' '  # space char added for postmers when prefix is =m
        self.alphabet: set[str] = alphabet

        self.number_of_greater_letters: dict[str, int]
        if number_of_greater_letters_dic is None:
            self.number_of_greater_letters = number_of_greater_letters(alphabet)
        else:
            self.number_of_greater_letters = number_of_greater_letters_dic
        self.antemer_max_prefix_size: int = self.length
        self.postmer_max_size: int = cmath.inf

        self.prefix_letters_vectors: list[dict[str, int]] = [dict() for _ in range(self.length + 1)]
        self.a_max: list[str] = [' '] * (self.length + 1)
        self.autocorrelation_matrix: list[list[str]] = []

        for i in range(self.length):
            self.autocorrelation_matrix.append(["="] + [None] * i)

            self.a_max[i + 1] = set(self.minimizer[0])

            self.prefix_letters_vectors[i + 1] = {a: cmath.inf for a in self.alphabet}

        self.prefix_letters_vectors[1] = {a: 0 for a in self.alphabet}
        self.prefix_letters_vectors[1][self.minimizer[0]] = 2

        for j in range(1, self.length):
            flag = False
            symb = None

            for i in range(j, self.length):
                if not flag:
                    cand = string[j:i + 1]
                    ref = string[:i - j + 1]
                    if cand == ref:
                        symb = "="
                        self.prefix_letters_vectors[i + 1][self.minimizer[i - j + 1]] = min(
                            self.prefix_letters_vectors[i + 1][self.minimizer[i - j + 1]], j + 1)
                        self.a_max[i + 1].add(self.minimizer[i - j + 1])
                    elif cand < ref:
                        symb = "<"
                        flag = True
                        self.antemer_max_prefix_size = min(self.antemer_max_prefix_size, i + 1)
                        self.postmer_max_size = min(self.postmer_max_size, j + 1)

                    else:
                        symb = '>'
                        flag = True

                self.autocorrelation_matrix[i][j] = symb

            for a in alphabet:
                if self.prefix_letters_vectors[j + 1][a] == cmath.inf:
                    self.prefix_letters_vectors[j + 1][a] = (j + 2) * (a == string[0])

        for i in range(self.length):
            self.a_max[i + 1] = max(self.a_max[i + 1])

        self.postmer_small_values_alphabet_zero : list[set[str]] = [set() for _ in range(self.length+1)]
        self.postmer_small_values_alphabet_not_zero : list[set[str]] = [set() for _ in range(self.length+1)]

        self.postmer_prefix_letters_vectors: list[dict[str, Callable[[int],int]]] = [dict() for _ in range(self.length + 1)]
        self.postmer_a_max: list[Callable[[int],str]] = [' '] * (self.length + 1)

        for i in range(1,self.length+1):
            a_max_dic = {0: ' '}
            for a in self.alphabet:
                if self.prefix_letters_vectors[i][a]==0:
                    self.postmer_small_values_alphabet_zero[i].add(a)
                    self.postmer_prefix_letters_vectors[i][a] = (lambda x:0)
                else:
                    self.postmer_small_values_alphabet_not_zero[i].add(a)
                    self.postmer_prefix_letters_vectors[i][a] = (lambda x, t=self.prefix_letters_vectors[i][a], m =self.length: t * (x >= m + t - 1))
                    a_max_dic[self.length + self.prefix_letters_vectors[i][a]-1] = a

            self.postmer_small_values_alphabet_zero[i].discard(self.minimizer[i])
            self.postmer_small_values_alphabet_not_zero[i].discard(self.minimizer[i])

            beta_values = sorted(a_max_dic.keys())
            intervals = []
            for j in range(len(beta_values)):
                intervals.append(max([a_max_dic[b] for b in beta_values[:j+1]]))

            self.postmer_a_max[i]= (lambda x,inter=tuple(intervals), bval = tuple(beta_values): ''.join([inter[j] * (bval[j] <= x < bval[j+1]) for j in range(len(bval)-1)])+ inter[-1] * (x >= bval[-1]))

    def antemer_lower_bound(self, alpha):
        array = [0] * (alpha + 1)

        array[0] = 1

        for j in range(1, alpha + 1):

            array[j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]

            for i in range(1, self.antemer_max_prefix_size):

                if j - (i + 1) >= 0:
                    array[j] += min(self.number_of_greater_letters[self.minimizer[i]],
                                    self.number_of_greater_letters[self.a_max[i]]) * array[
                                    j - (i + 1)]

        return array

    def antemer_upper_bound(self, alpha):
        array = [0] * (alpha + 1)

        array[0] = 1

        for j in range(1, alpha + 1):

            array[j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]

            for i in range(1, self.antemer_max_prefix_size):

                if j - (i + 1) >= 0:
                    array[j] += min(self.number_of_greater_letters[self.minimizer[i]],
                                    self.number_of_greater_letters[self.a_max[i]]) * array[
                                    j - (i + 1)]

                if self.a_max[i] > self.minimizer[i]:
                    if j - self.prefix_letters_vectors[i][self.a_max[i]] >= 0:
                        array[j] += array[j - self.prefix_letters_vectors[i][self.a_max[i]] + 1] - \
                                    self.number_of_greater_letters[self.minimizer[0]] * array[
                                        j - self.prefix_letters_vectors[i][self.a_max[i]]]
                    elif j - self.prefix_letters_vectors[i][self.a_max[i]] + 1 >= 0:
                        array[j] += array[j - self.prefix_letters_vectors[i][self.a_max[i]] + 1]
        return array

    def antemer(self, alpha):

        array_prefix = []
        array = [0] * (alpha + 1)

        for i in range(self.antemer_max_prefix_size):
            array_prefix.append([0] * (alpha + 1))

        for j in range(alpha + 1):
            for i in range(self.antemer_max_prefix_size):

                if i > j:
                    array_prefix[i][j] = 0
                elif j == 0:
                    array_prefix[i][j] = 1
                elif i == 0:
                    array_prefix[i][j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]
                elif i == j:
                    prod = 1
                    for l in range(i):
                        prod *= (self.autocorrelation_matrix[i - 1][l] == ">") or (
                                self.autocorrelation_matrix[i - 1][l] == '=' and
                                self.autocorrelation_matrix[self.length - 1][i - l] == '<')
                    array_prefix[i][j] = prod
                else:

                    array_prefix[i][j] = min(self.number_of_greater_letters[self.minimizer[i]],
                                             self.number_of_greater_letters[self.a_max[i]]) * array[
                                             j - (i + 1)]

                    if self.a_max[i] > self.minimizer[i]:
                        for new_prefix in range(i - self.prefix_letters_vectors[i][self.a_max[i]] + 2,
                                                self.antemer_max_prefix_size):
                            array_prefix[i][j] += array_prefix[new_prefix][
                                j - self.prefix_letters_vectors[i][self.a_max[i]] + 1]

            array[j] = sum([array_prefix[l][j] for l in range(self.antemer_max_prefix_size)])

        return array

    def postmer_lower_bound(self, beta):

        array = [0] * (beta + 1)
        array_prefix = [0] * (beta + 1)

        for j in range(self.length):
            array[j] = len(self.alphabet) ** j

        array[self.length] = self.number_of_greater_words
        array_prefix[self.length] = 1

        for j in range(self.length + 1, beta + 1):

            array[j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]

            for i in range(1, self.length + 1):

                a = self.postmer_a_max[i](j)

                array[j] +=  min(self.number_of_greater_letters[self.minimizer[i]], self.number_of_greater_letters[a]) * array[j - (i + 1)]
                if i == self.length:
                    array_prefix[j] += min(self.number_of_greater_letters[self.minimizer[i]], self.number_of_greater_letters[a]) * array[j - (self.length + 1)]

        return array_prefix

    def postmer_upper_bound(self, beta):
        array = [0] * (beta + 1)
        array_prefix = [0] * (beta + 1)

        for j in range(self.length):
            array[j] = len(self.alphabet) ** j

        array[self.length] = self.number_of_greater_words
        array_prefix[self.length] = 1

        for j in range(self.length + 1, beta + 1):

            array[j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]

            for i in range(1, self.length + 1):

                a = self.postmer_a_max[i](j)

                array[j] += min(self.number_of_greater_letters[self.minimizer[i]],
                                             self.number_of_greater_letters[a]) * array[j - (i + 1)]

                if a > self.minimizer[i]:
                    array[j] += array[j - self.postmer_prefix_letters_vectors[i][a](j) + 1] - self.number_of_greater_letters[self.minimizer[0]] * array[j - self.postmer_prefix_letters_vectors[i][a](j)]

                if i == self.length:
                    array_prefix[j] += min(self.number_of_greater_letters[self.minimizer[i]],
                                             self.number_of_greater_letters[self.postmer_a_max[i](j)]) * array[j - (self.length + 1)]
                    if a > self.minimizer[i]:
                        array_prefix[j] += array[j - self.postmer_prefix_letters_vectors[i][a](j) + 1] - self.number_of_greater_letters[self.minimizer[0]] * array[j - self.postmer_prefix_letters_vectors[i][a](j)]

        return array_prefix

    def postmer(self, beta):
        array_prefix = []
        array = [0] * (beta + 1)
        for i in range(self.length + 1):
            array_prefix.append([0] * (beta + 1))

        for j in range(beta + 1):
            for i in range(self.length + 1):
                if i > j:
                    array_prefix[i][j] = 0
                elif j == 0:
                    array_prefix[i][j] = 1
                elif j < self.length:
                    if i == 0:
                        array_prefix[i][j] = (len(self.alphabet) - 1) * array[j - 1]
                    elif i == j:
                        array_prefix[i][j] = 1
                    else:
                        array_prefix[i][j] = len(self.postmer_small_values_alphabet_zero[i]) * array[j - (i + 1)]

                        for a in self.postmer_small_values_alphabet_not_zero[i]:
                            for new_prefix in range(i - self.prefix_letters_vectors[i][a] + 2,
                                                    j - self.prefix_letters_vectors[i][a] + 2):
                                array_prefix[i][j] += array_prefix[new_prefix][
                                    j - self.prefix_letters_vectors[i][a] + 1]

                elif j == self.length:
                    if i == self.length:
                        array_prefix[i][j] = 1
                    else:
                        array_prefix[i][j] = self.number_of_greater_letters[self.minimizer[i]] * len(self.alphabet) ** (
                                    j - (i + 1))
                else:
                    if i == 0:
                        array_prefix[i][j] = self.number_of_greater_letters[self.minimizer[0]] * array[j - 1]
                    else:

                        a = self.postmer_a_max[i](j)

                        array_prefix[i][j] = min(self.number_of_greater_letters[self.minimizer[i]],
                            self.number_of_greater_letters[a]) * array[j - (i + 1)]

                        if a > self.minimizer[i]:
                            for new_prefix in range(i-self.postmer_prefix_letters_vectors[i][a](j)+2,self.length+1):
                                array_prefix[i][j] += array_prefix[new_prefix][j - self.postmer_prefix_letters_vectors[i][a](j) +1]

            array[j] = sum([array_prefix[l][j] for l in range(self.length + 1)])

        return array_prefix[self.length]

    def kmer(self, k):

        assert k >= self.length, 'k must be larger or equal to the length of the minimizer'

        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer(k - self.length)
        postmer_array = self.postmer(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return sum

    def kmer_upper_bound(self, k):
        assert k >= self.length, 'k must be larger or equal to the length of the minimizer'
        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer_upper_bound(k - self.length)
        postmer_array = self.postmer_upper_bound(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return min((beta_max + 1) * 4 ** (k - self.length), sum)

    def kmer_lower_bound(self, k):
        assert k >= self.length, 'k must be larger or equal to the length of the minimizer'
        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer_lower_bound(k - self.length)
        postmer_array = self.postmer_lower_bound(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return max(1, sum)

class VigeminCountingFunction(object):
    '''
    Counting functions for vigemins
    '''
    __slots__ = ('autocorrelation_matrix', 'alphabet', 'minimizer', 'length',
                 'antemer_max_prefix_size', 'postmer_max_size', 'prefix_letters_vectors','postmer_high_values_alphabet','postmer_small_values_alphabet_zero', 'postmer_small_values_alphabet_not_zero','key','suffix_key_convolution','antemer_alphabet_zero','antemer_alphabet_not_zero','alphabet_i')

    def __init__(self, string,key, alphabet={'A', 'T', 'C', 'G'}):

        assert len(key)==len(string)
        self.key: str = key + ' '
        self.length: int = len(string)
        self.minimizer: str = string + ' '  # space char added for postmers when prefix is =m
        self.alphabet: set[str] = alphabet

        self.antemer_max_prefix_size: int = self.length
        self.postmer_max_size: int = cmath.inf

        self.prefix_letters_vectors: list[dict[str, int]] = [dict() for _ in range(self.length + 1)]
        self.autocorrelation_matrix: list[list[str]] = []

        self.alphabet_i : list[set[str]] = [0]*(self.length+1)

        for i in range(self.length):
            self.alphabet_i[i]=set()
            for a in self.alphabet:
                if xor_char(a,self.key[i]) > xor_char(self.minimizer[i],self.key[i]):
                    self.alphabet_i[i].add(a)

        self.alphabet_i[self.length]=self.alphabet

        for i in range(self.length):
            self.autocorrelation_matrix.append(["="] + [None] * i)

            self.prefix_letters_vectors[i + 1] = {a: cmath.inf for a in self.alphabet}

        self.prefix_letters_vectors[1] = {a: 0 for a in self.alphabet}
        self.prefix_letters_vectors[1][self.minimizer[0]] = 2

        for j in range(1, self.length):
            flag = False
            symb = None

            for i in range(j, self.length):
                if not flag:
                    cand = xor_word(string[j:i + 1],key[:i-j+1])
                    ref = xor_word(string[:i - j + 1],key[:i-j+1])

                    if cand == ref:
                        symb = "="
                        self.prefix_letters_vectors[i + 1][self.minimizer[i - j + 1]] = min(
                            self.prefix_letters_vectors[i + 1][self.minimizer[i - j + 1]], j + 1)
                    elif cand < ref:
                        symb = "<"
                        flag = True
                        self.antemer_max_prefix_size = min(self.antemer_max_prefix_size, i + 1)
                        self.postmer_max_size = min(self.postmer_max_size, j + 1)

                    else:
                        symb = '>'
                        flag = True

                self.autocorrelation_matrix[i][j] = symb

            for a in alphabet:
                if self.prefix_letters_vectors[j + 1][a] == cmath.inf:
                    self.prefix_letters_vectors[j + 1][a] = (j + 2) * (a == string[0])

        self.suffix_key_convolution: list[bool] = [False]*self.antemer_max_prefix_size

        for i in range(self.antemer_max_prefix_size):
            if xor_word(self.minimizer[:self.length-i-1]+' ',self.key[i+1:]) > xor_word(self.minimizer[i+1:],self.key[i+1:]):
                self.suffix_key_convolution[i] = True

        self.antemer_alphabet_zero: list[set[str]] = [set() for _ in range(self.length+1)]
        self.antemer_alphabet_not_zero: list[set[str]] = [set() for _ in range(self.length+1)]

        for i in range(1, self.length):

            correct_alphabet = self.alphabet_i[i] & (self.alphabet_i[0] | {self.minimizer[0]})
            for j in range(1,i):
                if self.autocorrelation_matrix[i-1][j] == '=':
                    correct_alphabet = correct_alphabet & (self.alphabet_i[i-j] | {self.minimizer[i-j]})

            for a in correct_alphabet:
                if self.prefix_letters_vectors[i][a]==0:
                    self.antemer_alphabet_zero[i].add(a)
                else:
                    self.antemer_alphabet_not_zero[i].add(a)

        self.postmer_small_values_alphabet_zero : list[set[str]] = [set() for _ in range(self.length+1)]
        self.postmer_small_values_alphabet_not_zero : list[set[str]] = [set() for _ in range(self.length+1)]
        self.postmer_high_values_alphabet: list[dict[int, set[str]]] = [dict() for _ in range(self.length+1)]

        for i in range(1,self.length+1):
            for a in self.alphabet:
                if self.prefix_letters_vectors[i][a]==0:
                    self.postmer_small_values_alphabet_zero[i].add(a)
                else:
                    self.postmer_small_values_alphabet_not_zero[i].add(a)

            self.postmer_small_values_alphabet_zero[i].discard(self.minimizer[i])
            self.postmer_small_values_alphabet_not_zero[i].discard(self.minimizer[i])

            self.postmer_high_values_alphabet[i][i] = (self.alphabet_i[0] | {self.minimizer[0]})
            for k in range(1, i):
                if self.autocorrelation_matrix[i - 1][k] == '=':
                    self.postmer_high_values_alphabet[i][k] = (self.alphabet_i[i - k] | {self.minimizer[i - k]})

    def antemer(self, alpha):

        array_prefix = []
        array = [0] * (alpha + 1)

        for i in range(self.antemer_max_prefix_size):
            array_prefix.append([0] * (alpha + 1))

        for j in range(alpha + 1):
            for i in range(self.antemer_max_prefix_size):

                if i > j:
                    array_prefix[i][j] = 0
                elif j == 0:
                    array_prefix[i][j] = 1
                elif i == 0:
                    array_prefix[i][j]= len(self.alphabet_i[0]) * array[j-1]
                elif i == j:
                    prod = 1
                    for l in range(i):
                        prod *= (self.autocorrelation_matrix[i - 1][l] == ">") or (
                                self.autocorrelation_matrix[i - 1][l] == '=' and
                                self.suffix_key_convolution[i-1-l])

                    array_prefix[i][j] = prod
                else:
                    array_prefix[i][j] = len(self.antemer_alphabet_zero[i]) * array[j - (i + 1)]

                    for a in self.antemer_alphabet_not_zero[i]:
                        for new_prefix in range(i - self.prefix_letters_vectors[i][a] + 2, self.antemer_max_prefix_size):
                            array_prefix[i][j] += array_prefix[new_prefix][
                                j - self.prefix_letters_vectors[i][a] + 1]

            array[j] = sum([array_prefix[l][j] for l in range(self.antemer_max_prefix_size)])

        return array

    def postmer(self, beta):
        array_prefix = []
        array = [0] * (beta + 1)
        for i in range(self.length + 1):
            array_prefix.append([0] * (beta + 1))

        for j in range(beta + 1):
            for i in range(self.length + 1):
                if i > j:
                    array_prefix[i][j] = 0
                elif j == 0:
                    array_prefix[i][j] = 1
                elif j < self.length:
                    if i == 0:
                        array_prefix[i][j] = (len(self.alphabet) - 1) * array[j - 1]
                    elif i == j:
                        array_prefix[i][j] = 1
                    else:
                        array_prefix[i][j] = len(self.postmer_small_values_alphabet_zero[i]) * array[j - (i + 1)]

                        for a in self.postmer_small_values_alphabet_not_zero[i]:
                            for new_prefix in range(i - self.prefix_letters_vectors[i][a] + 2,
                                                    j - self.prefix_letters_vectors[i][a] + 2):
                                array_prefix[i][j] += array_prefix[new_prefix][
                                    j - self.prefix_letters_vectors[i][a] + 1]

                elif j == self.length:
                    if i == self.length:
                        array_prefix[i][j] = 1
                    else:
                        array_prefix[i][j] = len(self.alphabet_i[i]) * len(self.alphabet) ** (
                                    j - (i + 1))
                else:

                    if i == 0:
                        array_prefix[i][j] = len(self.alphabet_i[0]) * array[j - 1]
                    else:

                        correct_alphabet = self.alphabet_i[i]
                        for key in self.postmer_high_values_alphabet[i].keys():
                            if j>=self.length+key:
                                correct_alphabet = correct_alphabet & self.postmer_high_values_alphabet[i][key]

                        correct_alphabet_zero=set()
                        correct_alphabet_not_zero = set()

                        for a in correct_alphabet:
                            if j>=self.length-1+self.prefix_letters_vectors[i][a] and self.prefix_letters_vectors[i][a] != 0:
                                correct_alphabet_not_zero.add(a)
                            else:
                                correct_alphabet_zero.add(a)

                        array_prefix[i][j] = len(correct_alphabet_zero) * array[j - (i + 1)]

                        for a in correct_alphabet_not_zero:
                            for new_prefix in range(i-self.prefix_letters_vectors[i][a]+2,self.length+1):
                                array_prefix[i][j] += array_prefix[new_prefix][j - self.prefix_letters_vectors[i][a] +1]

            array[j] = sum([array_prefix[l][j] for l in range(self.length + 1)])

        return array_prefix[self.length]

    def kmer(self, k):

        assert k >= self.length, 'k must be larger or equal to the length of the minimizer'

        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer(k - self.length)
        postmer_array = self.postmer(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return sum

