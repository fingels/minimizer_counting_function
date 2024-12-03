from src.lib import *
import matplotlib.pyplot as plt
from matplotlib import rc
import cmath
from scipy.stats import linregress

minimizer_list = ['AAAAAA','ACACAA','ACACAC','CAAAAA','GAAAAA','TAAAAA']

m=max([len(mini) for mini in minimizer_list])
greater_letters_dic = number_of_greater_letters()

obj_list = {mini:MinimizerCountingFunction(mini,number_of_greater_letters_dic=greater_letters_dic) for mini in minimizer_list}

x = []
y = {mini: [] for mini in minimizer_list}
logy = {mini:[] for mini in minimizer_list}

kmax= 100

for k in range(m,m+kmax):
    x.append(k)

    for mini in minimizer_list:

        val = obj_list[mini].kmer(k)

        y[mini].append(val)
        logy[mini].append(cmath.log(val,4))

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

legend = []

for mini in minimizer_list:
    ax.plot(x,y[mini])

    legend.append(r'\texttt{'+mini+'}')


leg = ax.legend(legend,fontsize=fontsize,loc='upper left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

ax.xaxis.get_major_locator().set_params(integer=True)
plt.xticks(fontsize=fontsize)

plt.show()
fig.tight_layout()

fig.savefig('../Figures_theory/k_values_examples.pdf')

#####################################

# Regression

fig, ax = plt.subplots(figsize=(12, 8))

legend = []

i=0

for mini in minimizer_list:

    col = 'C'+str(i)
    i+=1

    ax.plot(x,logy[mini],c=col)
    slope, intercept, r, _, _ = linregress(x, logy[mini])
    reg = [intercept+slope*i for i in x]
    ax.plot(x,reg,c=col,ls='--',label='_nolegend_')

    legend.append(r'\texttt{'+mini+'}')


leg = ax.legend(legend,fontsize=fontsize,loc='upper left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.xaxis.get_major_locator().set_params(integer=True)
plt.xticks(fontsize=fontsize)


plt.show()
fig.tight_layout()
fig.savefig('../Figures_theory/k_values_examples_with_regression.pdf')