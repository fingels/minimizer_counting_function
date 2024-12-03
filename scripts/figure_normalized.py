import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import math

#####################################

# km_couples = [(21,2),(21,3),(21,4)]
# km_couples = [(21,3),(31,3),(61,3),(101,3)]

km_couples = [(61,8),(61,10),(101,8),(101,10)]
# km_couples = [(21,10),(31,10),(61,10),(101,10)]

#####################################

plot_dic = {}

for k,m in km_couples:

    with open('../Data/enumeration_k='+str(k)+'_m='+str(m)+'.p', "rb") as f:
        minimizer_dic = pickle.load(f)

    minimizers_list = sorted(minimizer_dic.keys())

    normalization_factor = (k-m) + math.log(k-m+1,4)

    x_coordinates = [i / (4 ** m -1) for i in range(len(minimizers_list))]
    y_coordinates = []

    for minimizer in minimizers_list:
        y_coordinates.append(math.log(minimizer_dic[minimizer],4)/normalization_factor)

    plot_dic[(k,m)] = {'x' : x_coordinates, 'y': y_coordinates, 'hline': (k-m)/normalization_factor}

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

legend = []

for k,m in plot_dic.keys():

    ax.plot(plot_dic[(k,m)]['x'],plot_dic[(k,m)]['y'],lw=0.5,alpha=0.5)

    legend.append(r'$k='+str(k)+', m='+str(m)+'$')

    # plt.axhline(plot_dic[(k,m)]['hline'])

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

leg = ax.legend(legend, fontsize=fontsize, loc='lower left')

for line in leg.get_lines():
    line.set_linewidth(4.0)

plt.show()
fig.tight_layout()
fig.savefig('../Figures_theory/normalized.pdf')

#####################################