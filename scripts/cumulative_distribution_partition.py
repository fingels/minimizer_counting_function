from src.lib import *
import time
import sys
import matplotlib.pyplot as plt
from matplotlib import rc

def complement(string):
    dic = {"A":"T", 'C':'G', 'G':'C', 'T':'A'}
    new_string = ''
    for char in string:
        new_string+=dic[char]
    return new_string

#####################################

k = 21
m = 5

#####################################

greater_letters_dic = number_of_greater_letters()

minimizer_dic = {}

bound_up = {}
bound_low = {}

equal_bounds =0
tight_bounds = 0

somme = 0

t = time.time()
for i in range(4**m):
    sys.stdout.write("\rProgression: %i/%i" % (i+1,4**m))
    minimizer = int_to_kmer(i,m)

    obj = MinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)

    N = obj.kmer(k)

    somme+=N

    minimizer_dic[minimizer]= N

assert somme == 4**k

t = time.time()-t

print('\nDone in %s seconds' % t)

#
# with open('../Data/enumeration_k='+str(k)+'_m='+str(m)+'.p', 'wb') as f:
#     pickle.dump(minimizer_dic, f)

minimizers_list = sorted(minimizer_dic.keys())

x = list(range(len(minimizers_list)))

prev = 0
cumulative_count = []
compl_sum = []

for minimizer in minimizers_list:

    temp = prev + minimizer_dic[minimizer]

    cumulative_count.append(temp)

    prev = temp

    compl = complement(minimizer)

    compl_sum.append(minimizer_dic[minimizer]+minimizer_dic[compl])


#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

# plt.axhline(y=4**(k-m),c='C3',label='_nolegend_')
# ax.plot(x,cumulative_count,lw=0.5)
ax.plot(x,compl_sum,lw=0.5)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

# plt.text(len(minimizers_list)-1, 4**(k-m+0.5), r'$4^{k-m}$', ha='right', va='center',fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

plt.show()
fig.tight_layout()

# fig.savefig('../Figures_theory/enumeration_k=' + str(k) + '_m=' + str(m) + '.pdf')
