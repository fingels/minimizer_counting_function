import pickle
import matplotlib.pyplot as plt
from matplotlib import rc

#####################################

k = 31
m = 10
compute_bounds = True

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

up_err = []
low_err = []

for minimizer in minimizers_list:

    count.append(minimizer_dic[minimizer])

    if compute_bounds:
        up.append(bound_up[minimizer])
        low.append(bound_low[minimizer])

        up_err.append((bound_up[minimizer]-minimizer_dic[minimizer])/minimizer_dic[minimizer])
        low_err.append((minimizer_dic[minimizer]-bound_low[minimizer])/minimizer_dic[minimizer])

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

    leg = ax.legend([r'$\pi_k(w)$', r'$\pi^+_k(w)$', r'$\pi^-_k(w)$'],fontsize=fontsize,loc='lower left')

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
    fig.savefig('../Figures_theory/bounds_k=' + str(k) + '_m=' + str(m) + '.pdf')
else:
    fig.savefig('../Figures_theory/enumeration_k=' + str(k) + '_m=' + str(m) + '.pdf')

#####################################t

if compute_bounds:

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(x,up_err,lw=0.5,c='C1')
    ax.plot(x,low_err,lw=0.5,c='C2')

    leg = ax.legend([r'$\eta^+(w)$', r'$\eta^-(w)$'],fontsize=fontsize,loc='upper right')

    for line in leg.get_lines():
        line.set_linewidth(4.0)

    n = len(minimizers_list) //4
    x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list)-1]
    x_ticks_labels = [r'\texttt{'+minimizers_list[i][:3]+'}$\cdots$' for i in x_ticks_positions]
    plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

    ax.yaxis.get_major_locator().set_params(integer=True)
    plt.yticks(fontsize=fontsize)
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    ax.set_yscale("log",base=10)

    plt.show()
    fig.tight_layout()
    fig.savefig('../Figures_theory/error_rate_bounds_k=' + str(k) + '_m=' + str(m) + '.pdf')