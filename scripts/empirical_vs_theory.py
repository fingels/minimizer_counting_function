import csv
import matplotlib.pyplot as plt
from matplotlib import rc
import math
from src.lib import *
import sys

###########################################################################

# name = "fusion61"
# name = "chrY_k31"
# name = "ecoli_k21"
name = "Hg_chr1"

imbalance_only= False

###########################################################################

filename = 'minimizers_results_'+name+'.csv'

dico = {}

sum = {}

with open('../Data/'+filename, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')

    for row in spamreader:

        k,m,minimizer,count = row

        if count != 'count':

            m= int(m)

            if m not in dico.keys():
                dico[m]={}
                sum[m]=0

            dico[m][minimizer]=int(count)
            sum[m]+=int(count)


k = int(k)
greater_letters_dic = number_of_greater_letters()

counts = {}
count_theo = {}

freq_empirical = {}
freq_theo = {}

x= {}
minimizers_list={}

for m in dico.keys():
    logsum = math.log(sum[m],4)

    minimizers_list[m] = sorted(dico[m].keys())

    x[m] = list(range(len(minimizers_list[m])))

    counts[m] = []
    count_theo[m] = []

    freq_empirical[m] = []
    freq_theo[m] = []

    #to tune the figures, parse only the beginning minimizers instead of all of them
    # N = 1000
    N = len(minimizers_list[m])

    i = 0

    print('Processing '+filename+' with the following values : k=%i, m=%i' % (k,m))

    for minimizer in minimizers_list[m][:N]:
        sys.stdout.write("\rProgression: %i/%i" % (i + 1, len(minimizers_list[m])))

        counts[m].append(dico[m][minimizer])

        if not imbalance_only:
            freq_empirical[m].append(math.log(dico[m][minimizer],4)-logsum)

            obj = MinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)
            th = obj.kmer(k)
            count_theo[m].append(th)
            freq_theo[m].append(math.log(th,4)-k)

        i+=1

    print('\n')

###########################################################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

###########################################################################

print('Building figures...')

for m in dico.keys():

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(x[m][:N],counts[m],lw=0.5)

    if not imbalance_only:
        ax.plot(x[m][:N],count_theo[m],lw=0.5)

        leg = ax.legend(['Empirical', 'Theory'],fontsize=fontsize,loc='center left')

        for line in leg.get_lines():
            line.set_linewidth(4.0)

    n = len(minimizers_list[m]) //4
    x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list[m])-1]
    x_ticks_labels = [r'\texttt{'+minimizers_list[m][i][:3]+'}$\cdots$' for i in x_ticks_positions]
    plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

    if imbalance_only:

        from matplotlib.ticker import ScalarFormatter

        class ScalarFormatterForceFormat(ScalarFormatter):
            def _set_format(self):  # Override function that finds format to use.
                self.format = "%1.1f"  # Give format here

        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        ax.yaxis.set_major_formatter(yfmt)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 2))
    else:
        ax.set_yscale("log", base=4)

    # ax.yaxis.get_major_locator().set_params(integer=True)
    plt.yticks(fontsize=fontsize)
    ax.yaxis.get_offset_text().set_fontsize(fontsize)

    # plt.show()
    fig.tight_layout()

    if imbalance_only:
        fig.savefig('../Figures_empirical/imbalance_'+name+'_k='+str(k)+'_m='+str(m)+'.pdf')
    else:
        fig.savefig('../Figures_empirical/theory_vs_'+name+'_k='+str(k)+'_m='+str(m)+'.pdf')

    ###########################################################################

    if not imbalance_only:

        fig, ax = plt.subplots(figsize=(12, 8))

        ax.plot(x[m][:N],freq_empirical[m],alpha=0.75,lw=0.75)
        ax.plot(x[m][:N],freq_theo[m],alpha=0.5,lw=0.5)

        leg = ax.legend(['Empirical', 'Theory'],fontsize=fontsize,loc='lower left')

        for line in leg.get_lines():
            line.set_linewidth(4.0)

        n = len(minimizers_list[m]) //4
        x_ticks_positions = [0,n,2*n,3*n,len(minimizers_list[m])-1]
        x_ticks_labels = [r'\texttt{'+minimizers_list[m][i][:3]+'}$\cdots$' for i in x_ticks_positions]
        plt.xticks(x_ticks_positions,x_ticks_labels,fontsize=fontsize)

        ax.yaxis.get_major_locator().set_params(integer=True)
        plt.yticks(fontsize=fontsize)
        ax.yaxis.get_offset_text().set_fontsize(fontsize)

        # plt.show()
        fig.tight_layout()
        fig.savefig('../Figures_empirical/frequences_'+name+'_k='+str(k)+'_m='+str(m)+'.pdf')

print('... done.')