import time
import sys
from src.lib import *
import matplotlib.pyplot as plt
from matplotlib import rc
compl = {'A': 'T', 'C':'G', 'G':'C', 'T':'A'}

def reverse_complement(kmer):
    return ''.join([compl[a] for a in kmer][::-1])

def canonical(kmer):
    return min(kmer,reverse_complement(kmer))

###########################################################################

k = 12
m = 4

###########################################################################

minimizer_dic = {}
c_minimizer_dic = {}

deja_vu = set()

t = time.time()
for i in range(4**k):
    sys.stdout.write("\rProgression: %i/%i" % (i+1,4**k))

    # regular kmer

    kmer = int_to_kmer(i,k)

    i_min = find_minimizer(kmer,m)
    minimizer = kmer[i_min:i_min+m]

    if minimizer not in minimizer_dic.keys():
        minimizer_dic[minimizer]=0

    minimizer_dic[minimizer]+=1

    # canonical kmer

    can = canonical(kmer)

    if can not in deja_vu:

        deja_vu.add(can)

        ci_min = find_minimizer(can,m)
        c_minimizer = can[ci_min:ci_min+m]

        if c_minimizer not in c_minimizer_dic.keys():
            c_minimizer_dic[c_minimizer]=0

        c_minimizer_dic[c_minimizer]+=1

t = time.time()-t

print('\nDone in %s seconds' % t)

minimizers_list = sorted(minimizer_dic.keys())
x = list(range(len(minimizers_list)))

count = []
c_count = []

for minimizer in minimizers_list:

    count.append(minimizer_dic[minimizer])

    if minimizer in c_minimizer_dic.keys():
        c_count.append(c_minimizer_dic[minimizer])
    else:
        c_count.append(1/4)

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

plt.axhline(y=4**(-1),c='C3',ls='--',alpha=0.5,label='_nolegend_')
ax.plot(x,count,lw=0.5)
ax.plot(x,c_count,lw=0.5)

leg = ax.legend(['Regular', 'Canonical'], fontsize=fontsize, loc='upper right')

for line in leg.get_lines():
    line.set_linewidth(4.0)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

plt.text(len(minimizers_list)-1, 4**(-1+0.25), 'Empty', ha='right', va='center',fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

plt.show()
fig.tight_layout()
fig.savefig('../Figures_theory/canonical_brute_force_k=' + str(k) + '_m=' + str(m) + '.pdf')