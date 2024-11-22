import csv
import matplotlib.pyplot as plt
from matplotlib import rc
import math
from src.lib import *
import sys

# filename = "minimizers_results_rna_fusion.csv"
# filename = "minimizers_results_chrY_k31.csv"
# filename = "minimizers_results_ecoli_k21.csv"
filename = "minimizers_results61.csv"

dico = {}

sum = 0

with open('../Data/'+filename, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')

    for row in spamreader:

        k,_,minimizer,count = row

        if count != 'count':
            dico[minimizer]=int(count)
            sum+=int(count)

k = int(k)
logsum = math.log(sum,4)

minimizers_list = sorted(dico.keys())

greater_letters_dic = number_of_greater_letters()

x = list(range(len(minimizers_list)))

counts = []
count_theo = []

freq_empirical = []
freq_theo = []

#to tune the figures, parse only the beginning minimizers instead of all of them
# N = 10000
N = len(minimizers_list)

i = 0

for minimizer in minimizers_list[:N]:
    sys.stdout.write("\rProgression: %i/%i" % (i + 1, len(minimizers_list)))

    counts.append(dico[minimizer])
    freq_empirical.append(math.log(dico[minimizer],4)-logsum)

    obj = MinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)
    th = obj.kmer(k)
    count_theo.append(th)
    freq_theo.append(math.log(th,4)-k)

    i+=1

###########################################################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

###########################################################################

fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(x[:N],counts,lw=0.5)
ax.plot(x[:N],count_theo,lw=0.5)

leg = ax.legend(['Empirical', 'Theory'],fontsize=fontsize,loc='center left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

plt.show()
fig.tight_layout()
# fig.savefig('../Images/imbalance_'+filename+'.pdf')
fig.savefig('../Images/theory_vs_'+filename+'.pdf')

###########################################################################

fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(x[:N],freq_empirical,alpha=0.75,lw=0.75)
ax.plot(x[:N],freq_theo,alpha=0.5,lw=0.5)

leg = ax.legend(['Empirical', 'Theory'],fontsize=fontsize,loc='center left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)

plt.show()
fig.tight_layout()
fig.savefig('../Images/frequences_'+filename+'.pdf')