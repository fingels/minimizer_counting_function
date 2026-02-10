from src.lib import LexMinimizerCountingFunction
from src.utils import *
import unittest
import os
import cmath

class TestRegularKmerEnumeration(unittest.TestCase):

    def setUp(self):
        self.data_dir = 'tmp/test/'
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        self.k = 10
        self.m = 5
        self.greater_letters_dic = number_of_greater_letters()
        self.minimizers = {}

    def tearDown(self):
        for file in os.listdir(self.data_dir):
            os.remove(self.data_dir+file)
        os.removedirs(self.data_dir)


    def test_enumerate_kmer_and_minimizers(self):
        for i in range(4**self.k):
            kmer = int_to_kmer(i, self.k)
            min_index = find_minimizer(kmer,self.m)
            found_minimizer = kmer[min_index:min_index+self.m]

            if found_minimizer not in self.minimizers.keys():
                self.minimizers[found_minimizer]=0

            self.minimizers[found_minimizer]+=1

        assert len(self.minimizers.keys())==4**self.m

    def test_found_minimizers_count(self):
        for minimizer in self.minimizers.keys():
            obj = LexMinimizerCountingFunction(minimizer, number_of_greater_letters_dic=self.greater_letters_dic)

            assert obj.kmer(self.k) == self.minimizers[minimizer]

            assert obj.kmer_upper_bound(self.k) >= obj.kmer(self.k)
            assert obj.kmer_lower_bound(self.k) <= obj.kmer(self.k)

class TestRegularKmerUnique(unittest.TestCase):

    def setUp(self):
        self.data_dir = 'tmp/test/'
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        self.greater_letters_dic = number_of_greater_letters()
        self.minimizer = 'ACACAA'
        self.alpha_value = 10
        self.beta_value = 10

        self.max_k_value = 31

        self.m = len(self.minimizer)

        self.obj = LexMinimizerCountingFunction(self.minimizer, number_of_greater_letters_dic=self.greater_letters_dic)

        if self.obj.postmer_max_size != cmath.inf:
            self.beta_value = self.obj.postmer_max_size - 1

        self.antemers = self.obj.antemer(self.alpha_value)
        self.postmers = self.obj.postmer(self.beta_value + self.m)

    def tearDown(self):
        for file in os.listdir(self.data_dir):
            os.remove(self.data_dir+file)
        os.removedirs(self.data_dir)

    def test_multiple_values(self):

        all_values = []

        for k in range(self.m,self.max_k_value+1):
            all_values.append(self.obj.kmer(k))

        all_values_bis = self.obj.kmer(self.max_k_value,return_all_values=True)

        assert tuple(all_values_bis)==tuple(all_values)

    def test_antemers(self):
        for alpha in range(self.alpha_value):
            candidates = []
            for i in range(4**alpha):
                kmer = int_to_kmer(i,alpha)
                candidates.append(kmer+self.minimizer)

            acceptables = []

            for kmer in candidates:
                min_index = find_minimizer(kmer,self.m)
                found_minimizer = kmer[min_index:min_index+self.m]
                if found_minimizer==self.minimizer and min_index==alpha:
                    acceptables.append(kmer)

            assert self.antemers[alpha]==len(acceptables)

    def test_postmers(self):
        for beta in range(self.m,self.beta_value+self.m):
            candidates = []
            for i in range(4**(beta-self.m)):
                kmer = int_to_kmer(i,beta-self.m)
                candidates.append(self.minimizer+kmer)

            acceptables = []

            for kmer in candidates:
                min_index = find_minimizer(kmer,self.m)
                found_minimizer = kmer[min_index:min_index+self.m]
                if found_minimizer>=self.minimizer:
                    acceptables.append(kmer)

            assert self.postmers[beta]==len(acceptables)