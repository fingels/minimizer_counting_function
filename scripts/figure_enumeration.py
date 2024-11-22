import pickle
import matplotlib.pyplot as plt
from matplotlib import rc

#####################################

k = 21
m = 8
compute_bounds = False

#####################################

with open('../Data/enumeration_k='+str(k)+'_m='+str(m)+'.p', "rb") as f:
    minimizer_dic = pickle.load(f)

if compute_bounds:
    with open('../Data/bounds_k=' + str(k) + '_m=' + str(m) + '.p', "rb") as f:
        bounds = pickle.load(f)
        bound_up = bounds['up']
        bound_low = bounds['low']

minimizers_list = sorted(minimizer_dic.keys())

x = list(range(len(minimizers_list)))

count = []

up = []
low = []

for minimizer in minimizers_list:

    count.append(minimizer_dic[minimizer])

    if compute_bounds:
        up.append(bound_up[minimizer])
        low.append(bound_low[minimizer])

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

plt.axhline(y=4**(k-m),c='C3',label='_nolegend_')
ax.plot(x,count,lw=0.5)

if compute_bounds:
    ax.plot(x,up,lw=0.5,ls=(0,(3,3)))
    ax.plot(x,low,lw=0.5,ls=(3,(3,3)))

    leg = ax.legend(['Value', 'Upper bound', 'Lower bound'],fontsize=fontsize,loc='lower left')

    for line in leg.get_lines():
        line.set_linewidth(4.0)

n = len(minimizers_list) //4
x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

plt.text(len(minimizers_list)-1, 4**(k-m+0.5), r'$4^{k-m}$', ha='right', va='center',fontsize=fontsize)

ax.yaxis.get_major_locator().set_params(integer=True)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.set_yscale("log",base=4)

plt.show()
fig.tight_layout()

if compute_bounds:
    fig.savefig('../Images/bounds_k=' + str(k) + '_m=' + str(m) + '.pdf')
else:
    fig.savefig('../Images/enumeration_k=' + str(k) + '_m=' + str(m) + '.pdf')

#####################################