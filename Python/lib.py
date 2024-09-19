import cmath
from utils import *


# TODO : refactor with integers instead of letters ?
# TODO : preprocess more the postmers

class PreProcess(object):
    __slots__ = 'autocorrelation_matrix', 'alphabet', 'minimizer', 'length', 'antemer_max_prefix_size', 'postmer_max_size', 'prefix_letters_vectors', 'a_max', 'number_of_greater_letters', 'number_of_greater_words'

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
        self.postmer_max_size: int | float = cmath.inf

        self.prefix_letters_vectors: list[dict[str, int]] = [dict()] * (self.length + 1)
        self.a_max: list[str] = [''] * (self.length + 1)
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
                tilde_minj = {}

                a_max = ''
                for a in self.alphabet:
                    tilde_minj[a] = self.prefix_letters_vectors[i][a] * (
                                j >= self.length + self.prefix_letters_vectors[i][a] - 1)
                    if tilde_minj[a] != 0:
                        a_max = max(a_max, a)

                alphabet_partition_A = []
                alphabet_partition_B = []
                for s in [a for a in self.alphabet if a > self.minimizer[i] and a >= a_max]:
                    if tilde_minj[s] == 0:
                        alphabet_partition_A.append(s)
                    else:
                        alphabet_partition_B.append(s)

                array[j] += len(alphabet_partition_A) * array[j - (i + 1)]

                if i == self.length:
                    array_prefix[j] += len(alphabet_partition_A) * array[j - (self.length + 1)]

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
                tilde_minj = {}

                a_max = ''
                for a in self.alphabet:
                    tilde_minj[a] = self.prefix_letters_vectors[i][a] * (
                                j >= self.length + self.prefix_letters_vectors[i][a] - 1)
                    if tilde_minj[a] != 0:
                        a_max = max(a_max, a)

                alphabet_partition_A = []
                alphabet_partition_B = []
                for s in [a for a in self.alphabet if a > self.minimizer[i] and a >= a_max]:
                    if tilde_minj[s] == 0:
                        alphabet_partition_A.append(s)
                    else:
                        alphabet_partition_B.append(s)

                array[j] += len(alphabet_partition_A) * array[j - (i + 1)]

                for a in alphabet_partition_B:
                    array[j] += array[j - tilde_minj[a] + 1] - self.number_of_greater_letters[self.minimizer[0]] * \
                                array[
                                    j - tilde_minj[a]]

                if i == self.length:
                    array_prefix[j] += len(alphabet_partition_A) * array[j - (self.length + 1)]
                    for a in alphabet_partition_B:
                        array_prefix[j] += array[j - tilde_minj[a] + 1] - self.number_of_greater_letters[
                            self.minimizer[0]] * array[
                                               j - tilde_minj[a]]

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

                        alphabet_zero = 0
                        alphabet_not_zero = []

                        for s in self.alphabet - {self.minimizer[i]}:
                            if self.prefix_letters_vectors[i][s] == 0:
                                alphabet_zero += 1
                            else:
                                alphabet_not_zero.append(s)

                        array_prefix[i][j] = alphabet_zero * array[j - (i + 1)]

                        for a in alphabet_not_zero:
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

                        tilde_minj = {}

                        a_max = ''
                        for a in self.alphabet:
                            tilde_minj[a] = self.prefix_letters_vectors[i][a] * (
                                        j >= self.length + self.prefix_letters_vectors[i][a] - 1)
                            if tilde_minj[a] != 0:
                                a_max = max(a_max, a)

                        alphabet_partition_A = []
                        alphabet_partition_B = []
                        for s in [a for a in self.alphabet if a > self.minimizer[i] and a >= a_max]:
                            if tilde_minj[s] == 0:
                                alphabet_partition_A.append(s)
                            else:
                                alphabet_partition_B.append(s)

                        assert len(alphabet_partition_B) <= 1

                        array_prefix[i][j] = len(alphabet_partition_A) * array[j - (i + 1)]
                        for a in alphabet_partition_B:
                            for new_prefix in range(i - tilde_minj[a] + 2, self.length + 1):
                                array_prefix[i][j] += array_prefix[new_prefix][j - tilde_minj[a] + 1]

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
        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer_upper_bound(k - self.length)
        postmer_array = self.postmer_upper_bound(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return min((beta_max + 1) * 4 ** (k - self.length), sum)

    def kmer_lower_bound(self, k):
        beta_max = min(self.postmer_max_size - 2, k - self.length)

        antemer_array = self.antemer_lower_bound(k - self.length)
        postmer_array = self.postmer_lower_bound(beta_max + self.length)

        sum = 0

        for beta in range(0, beta_max + 1):
            sum += antemer_array[k - self.length - beta] * postmer_array[beta + self.length]

        return max(1, sum)
