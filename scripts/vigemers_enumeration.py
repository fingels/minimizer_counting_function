from src.utils import *
from src.lib import VigeminCountingFunction as VCF
import sys
import time
import matplotlib.pyplot as plt
from matplotlib import rc

#####################################################

k = 31
m = 10
#####################################################

keys=[]

keys.append('A' * m)  # lexicographical order
keys.append('A' + 'T' * (m - 1))  # anti-lexicographical order
keys.append(('AT' * m)[:m])  # reverse lexicographical order>

for char in ['C', 'G', 'T']:
    keys.append(char + generate_random_k_mers(m - 1))

minimizer_dic = {}
somme = {}
for key in keys:
    minimizer_dic[key] = {}
    somme[key] =  0

minimizers_list = []

t = time.time()
for i in range(4**m):
    sys.stdout.write("\rProgression: %f%%" % ((i+1)/4**m*100))
    minimizer = int_to_kmer(i,m)

    minimizers_list.append(minimizer)

    for key in keys:

        obj = VCF(minimizer, key)

        N = obj.kmer(k)

        somme[key]+=N

        minimizer_dic[key][minimizer]= N

for key in keys:
    assert somme[key] == 4**k

t = time.time()-t
print("\nDone in %s seconds" % t)

minimizers_list = sorted(minimizers_list)

x = list(range(len(minimizers_list)))

count = {}

for key in keys:
    count[key] = []

    for minimizer in minimizers_list:

        if minimizer_dic[key][minimizer]==0:
            count[key].append(4**(-1))
        else:
            count[key].append(minimizer_dic[key][minimizer])

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(20, 5),sharey=True)

for ax in [ax1, ax2]:
    ax.axhline(y=4**(k-m),c='C3',label='_nolegend_')
    ax.axhline(y=4**(-1),c='C3',label='_nolegend_',ls='--',alpha=0.3)

ax1.text(1, 4**(-1), r'Empty', ha='left', va='bottom',fontsize=fontsize)
ax1.text(len(minimizers_list)-1, 4**(k-m), r'$4^{k-m}$', ha='right', va='bottom',fontsize=fontsize)
# ax2.text(1, 4**(k-m), r'$4^{k-m}$', ha='right', va='bottom',fontsize=fontsize)

i=0
for key in keys[:3]:
    ax1.plot(x,count[key],lw=0.5,color='C'+str(i),label=r'$\gamma=\texttt{'+key+'}$')
    i+=1

for key in keys[3:]:
    ax2.plot(x, count[key], lw=0.5,color='C'+str(i),label=r'$\gamma=\texttt{'+key+'}$')
    i+=1

for ax in [ax1, ax2]:
    n = len(minimizers_list) //4
    x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
    x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
    ax.set_xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

ax1.yaxis.get_major_locator().set_params(integer=True)
ax1.tick_params(axis='y', labelsize=fontsize)
ax1.yaxis.get_offset_text().set_fontsize(fontsize)
ax1.set_yscale("log",base=4)

leg = fig.legend(fontsize=fontsize, ncol=3,loc='upper center', bbox_to_anchor=(0.5, 0))

for line in leg.get_lines():
    line.set_linewidth(4.0)

plt.show()

fig.savefig('../Figures_theory/vigemin_enumeration_k=' + str(k) + '_m=' + str(m) + '.pdf',bbox_inches="tight")

####################################

fig, ax = plt.subplots(figsize=(12, 8))

plt.axhline(y=4**(k-m),c='C3',label='_nolegend_')
plt.text(len(minimizers_list)-1, 4**(k-m+0.5), r'$4^{k-m}$', ha='right', va='center',fontsize=fontsize)

ax.axhline(y=4**(-1),c='C3',label='_nolegend_',ls='--',alpha=0.3)
ax.text(1, 4**(-1), r'Empty', ha='left', va='bottom',fontsize=fontsize)

i=0
for key in keys:
    ax.plot(x,sorted(count[key],reverse=True),lw=0.5,color='C'+str(i),label=r'$\gamma=\texttt{'+key+'}$')
    i+=1

# n = len(minimizers_list) //4
# x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
# x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
# plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)
plt.xticks(fontsize=fontsize)

leg = fig.legend(fontsize=fontsize, ncol=2,loc='upper center', bbox_to_anchor=(0.5, 0))

for line in leg.get_lines():
    line.set_linewidth(4.0)

plt.show()
fig.savefig('../Figures_theory/vigemin_sorted_enumeration_k=' + str(k) + '_m=' + str(m) + '.pdf',bbox_inches="tight")
