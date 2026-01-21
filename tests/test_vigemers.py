from src.lib import VigeminCountingFunction
from src.utils import *
import unittest
import os
import random
import cmath

random.seed(1234)

class TestVigemerEnumeration(unittest.TestCase):

    def setUp(self):
        self.data_dir = 'tmp/test/'
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        self.k = 10
        self.m = 5
        self.key = generate_random_k_mers(self.m)
        self.vigemins = {}

    def tearDown(self):
        for file in os.listdir(self.data_dir):
            os.remove(self.data_dir+file)
        os.removedirs(self.data_dir)

    def test_enumerate_kmer_and_vigemins(self):

        for i in range(4**self.k):
            kmer = int_to_kmer(i, self.k)
            min_index = find_vigemin(kmer,self.key,self.m)
            found_vigemin = kmer[min_index:min_index+self.m]

            if found_vigemin not in self.vigemins.keys():
                self.vigemins[found_vigemin]=0

            self.vigemins[found_vigemin]+=1

    def test_found_vigemin_count(self):
        for vigemin in self.vigemins.keys():
            assert VigeminCountingFunction(vigemin,self.key).kmer(self.k) == self.vigemins[vigemin]

class TestVigemerUnique(unittest.TestCase):

    def setUp(self):
        self.data_dir = 'tmp/test/'
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        self.minimizer = 'ACACAA'
        self.key = 'CTCAAT'
        self.alpha_value = 10
        self.beta_value = 10

        self.m = len(self.minimizer)

        self.obj = VigeminCountingFunction(self.minimizer, self.key)

        if self.obj.postmer_max_size != cmath.inf:
            self.beta_value = self.obj.postmer_max_size - 1

        self.antemers = self.obj.antemer(self.alpha_value)
        self.postmers = self.obj.postmer(self.beta_value + self.m)

    def tearDown(self):
        for file in os.listdir(self.data_dir):
            os.remove(self.data_dir+file)
        os.removedirs(self.data_dir)

    def test_antemers(self):
        for alpha in range(self.alpha_value):
            candidates = []
            for i in range(4**alpha):
                kmer = int_to_kmer(i,alpha)
                candidates.append(kmer+self.minimizer)

            acceptables = []

            for kmer in candidates:
                min_index = find_vigemin(kmer,self.key,self.m)
                found_vigemin = kmer[min_index:min_index+self.m]
                if found_vigemin==self.minimizer and min_index==alpha:
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
                min_index = find_vigemin(kmer,self.key,self.m)
                found_vigemin = kmer[min_index:min_index+self.m]
                if xor_word(found_vigemin,self.key)>=xor_word(self.minimizer,self.key):
                    acceptables.append(kmer)

            assert self.postmers[beta]==len(acceptables)
