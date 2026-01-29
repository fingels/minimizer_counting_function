from src.lib import VigeminCountingFunction
from src.utils import *
import matplotlib.pyplot as plt
from matplotlib import rc
import math
import numpy as np
import time
import sys

def regress_and_predict(data_x,data_y,val_to_pred):

    cov_mat = np.cov(data_x, data_y)

    cov = cov_mat[0][1]
    vx = cov_mat[0][0]
    vy = cov_mat[1][1]
    mx = np.mean(data_x)
    my = np.mean(data_y)

    a = cov / vx
    b = my - a * mx
    #
    # r = cov / np.sqrt(vx * vy)

    return a*val_to_pred+b

#####################################

k_test = 100
k_min = 15
k_max = 25
m = 10
n_test = 10**3

#####################################

data_x = list(range(k_min,k_max))

keys=[]
for char in ['A','C', 'G', 'T']:
    keys.append(char + generate_random_k_mers(m - 1))

minimizer_dic = {}
for key in keys:
    minimizer_dic[key] = {}

minimizers_list = []

t = time.time()

for i in range(n_test):
    sys.stdout.write("\rProgression: %i / %i" % (i+1, n_test))

    minimizer = generate_random_k_mers(m)

    minimizers_list.append(minimizer)

    for key in keys:

        obj = VigeminCountingFunction(minimizer, key)

        data_y = []

        for k in data_x:

            val = obj.kmer(k)

            if val == 0:
                val = -1
            else:
                val = math.log(val, 4)

            data_y.append(val)

        val_pred = regress_and_predict(data_x,data_y,k_test)

        val_theo = obj.kmer(k_test)

        if val_theo == 0:
            val_theo = -1
        else:
            val_theo = math.log(val_theo, 4)

        minimizer_dic[key][minimizer] = (val_theo,val_pred)

t = time.time() - t
print("\nDone in %s seconds" % t)

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

minimizers_list = sorted(minimizers_list)

pred = {}
theo = {}

mint= float('inf')
maxt = 0

for key in keys:
    pred[key] = []
    theo[key] = []

    for minimizer in minimizers_list:

        val_theo, val_pred = minimizer_dic[key][minimizer]

        pred[key].append(val_pred)
        theo[key].append(val_theo)

        if val_theo>maxt:
            maxt = val_theo
        if val_theo<mint:
            mint = val_theo

x = np.linspace(mint,maxt,100)

#####################################

fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(x,x,color='black',alpha=0.5,ls="dashed")

i=0
for key in keys:
    ax.scatter(theo[key],pred[key],s=5,label=r"$\gamma= \texttt {" + key + "}$",color='C'+str(i))
    i+=1

ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlabel(r'$\log\pi_k^\gamma(w)$',fontsize=fontsize)
ax.set_ylabel(r'$\widehat{\log\pi_k^\gamma(w)}$',fontsize=fontsize)

leg = fig.legend(fontsize=fontsize, ncol=2,loc='upper center', bbox_to_anchor=(0.5, 0))

plt.show()
fig.savefig('../Figures_theory/vigemin_approximation_k=' + str(k_test) + '_m=' + str(m) + '.pdf',bbox_inches="tight")

