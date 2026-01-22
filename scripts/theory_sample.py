from src.lib import *
import time
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#####################################

k = 61
m = 21
n_sample = 100000

#####################################

greater_letters_dic = number_of_greater_letters()

minimizer_dic = {}

numbers = np.linspace(0,4**m-1, n_sample, dtype=int)

count = []

t = time.time()

j=0
for i in numbers:
    sys.stdout.write("\rProgression: %i/%i" % (j+1,len(numbers)))
    minimizer = int_to_kmer(i,m)

    obj = LexMinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)

    N = obj.kmer(k)

    minimizer_dic[minimizer]= N

    count.append(N)

    j+=1

t = time.time()-t

print('\nDone in %s seconds' % t)

with open('../Data/sample_k='+str(k)+'_m='+str(m)+'_N='+str(n_sample)+'.p', 'wb') as f:
    pickle.dump(minimizer_dic, f)

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30


#####################################

minimizers_list = sorted(minimizer_dic.keys())


fig, ax = plt.subplots(figsize=(12, 8))

plt.axhline(y=4**(k-m),c='C3',label='_nolegend_')
ax.plot(numbers,count,lw=0.5)

n = len(numbers) //4
x_ticks_refs = [0,n,2*n,3*n,len(numbers)-1]

x_ticks_positions= [numbers[i] for i in x_ticks_refs]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_refs]

plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

plt.text(numbers[-1], 4**(k-m+1), r'$4^{k-m}$', ha='right', va='center',fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

plt.show()
fig.tight_layout()

fig.savefig('../Figures_theory/sample_k=' + str(k) + '_m=' + str(m) +'_N='+str(n_sample)+ '.pdf')

#####################################t