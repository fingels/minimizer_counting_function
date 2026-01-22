from src.lib import *
import matplotlib.pyplot as plt
from matplotlib import rc
import cmath
from scipy.stats import linregress
import sys

#####################################

kmax = 100
m = 6

#####################################

greater_letters_dic = number_of_greater_letters()

minimizers_list = []


slopes = []
intercepts = []
rvals = []
errs = []

xx = list(range(m,m+kmax+1))

for i in range(4**m):
    sys.stdout.write("\rProgression: %i/%i" % (i+1,4**m))
    minimizer = int_to_kmer(i,m)

    minimizers_list.append(minimizer)

    obj = LexMinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)

    values =[]

    for k in xx:
       values.append(cmath.log(obj.kmer(k),4))

    slope, intercept, r, _, err = linregress(xx, values)

    errs.append(err)
    slopes.append(slope)
    intercepts.append(intercept)
    rvals.append(r**2)

minimizers_list.sort()
x = list(range(len(minimizers_list)))

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

plt.axhline(y=1,c='gray',label='_nolegend_',lw=0.5)

ax.plot(x,slopes)
ax.plot(x,rvals)

leg = ax.legend([r'$\widehat{A_w}$', r'$R^2$'], fontsize=fontsize, loc='lower left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

# ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)

plt.show()
fig.tight_layout()

fig.savefig('../Figures_theory/regression_slopes_enumeration.pdf')

#####################################

