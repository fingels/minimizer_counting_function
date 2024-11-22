import sys
import pickle
import time
from src.lib import *

minimizer = 'ACACAA'
k=10


m=len(minimizer)
greater_letters_dic = number_of_greater_letters()

obj = MinimizerCountingFunction(minimizer,number_of_greater_letters_dic=greater_letters_dic)

relmat = obj.autocorrelation_matrix
minj = obj.prefix_letters_vectors
a_max = obj.a_max
prefix_max_size = obj.antemer_max_prefix_size
postmer_max = obj.postmer_max_size

print(postmer_max)

print(obj.antemer_upper_bound(k))
print(obj.antemer(k))
print(obj.antemer_lower_bound(k))
print('\n')
print(obj.postmer_upper_bound(k+m)[m:])
print(obj.postmer(k+m)[m:])
print(obj.postmer_lower_bound(k+m)[m:])
print('\n')

for k in range(11):
    print(k+m,obj.kmer_lower_bound(k+m),obj.kmer(k+m),obj.kmer_upper_bound(k+m))



