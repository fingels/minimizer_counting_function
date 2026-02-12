from src.lib import *
import sys
import time
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc

###########################################################################

k = 31
m = 6
N_keys = 8
seq_size = 10**5

###########################################################################

keys = generate_random_k_mers(m,N_keys)

assert len(set(keys))==N_keys

###########################################################################

# TODO: SI ON VEUT IMPORTER UNE SEQUENCE DEPUIS UN FASTA IL FAUT CHANGER CES LIGNES

t = time.time()
print('Generating the sequence...')

S = np.random.choice(a=['A','C','G','T'],size=seq_size)
S = ''.join(S)

t = time.time()-t
print('Done in %f s.' % t)

###########################################################################

t=time.time()
print('\nParsing the sequence to compute the %i partitions...' % N_keys)

vigemin_dict = {key : {"pos": set(), "minimizers": dict()} for key in keys}

for i in range(seq_size-k+1):
    sys.stdout.write("\rProgression: %i/%i " % (i + 1, seq_size-k+1))

    kmer = S[i:i+k]

    ### COMPUTING THE VIGEMINS

    for key in keys:
        idx_vigemin = find_vigemin(kmer,key)
        vigemin = kmer[idx_vigemin:idx_vigemin+m]

        vigemin_dict[key]["pos"].add(i+idx_vigemin)

        if vigemin not in vigemin_dict[key]['minimizers'].keys():
            vigemin_dict[key]['minimizers'][vigemin] = 0

        vigemin_dict[key]['minimizers'][vigemin] += 1

t = time.time()-t
print('\nDone in %f s., hence %f s/key/kmer.' % (t,t/(seq_size-k+1)/N_keys))

for key in keys:
    vigemin_dict[key]["buckets"] = sorted(vigemin_dict[key]["minimizers"].values(),reverse=True)
    vigemin_dict[key]["density"] = len(vigemin_dict[key]["pos"]) / (seq_size-k+1)

print('Random minimizer density :', 2/(k-m+1))
print('Mean density for the keys :', np.mean([vigemin_dict[key]["density"] for key in keys]))

###########################################################################

t=time.time()
print('\nParsing the sequence to compute the heuristic...')

heuristic_dict = {"pos": set(), "minimizers": dict()}
heuristic_duplicated_dict = {"pos": set(), "minimizers": dict()}

for i in range(seq_size-k+1):
    sys.stdout.write("\rProgression: %i/%i " % (i + 1, seq_size-k+1))

    kmer = S[i:i+k]

    ### COMPUTING THE HEURISTIC

    tuple_list = []

    # le bloc qui suit est parallélisable
    for key in keys:
        idx_vigemin = find_vigemin(kmer,key)
        vigemin = kmer[idx_vigemin:idx_vigemin+m]

        oracle_val = VigeminCountingFunction(vigemin,key).kmer(k)

        tuple_list.append((oracle_val,idx_vigemin,vigemin,key))

    _, min_idx, min_vigemin, min_key = sorted(tuple_list,key= lambda tup: tup[0])[0]

    heuristic_dict["pos"].add(min_idx+i)
    heuristic_duplicated_dict["pos"].add(min_idx+i)

    if min_vigemin not in heuristic_dict["minimizers"].keys():
        heuristic_dict["minimizers"][min_vigemin] = 0

    heuristic_dict["minimizers"][min_vigemin]+=1

    if (min_key,min_vigemin) not in heuristic_duplicated_dict["minimizers"].keys():
        heuristic_duplicated_dict["minimizers"][(min_key,min_vigemin)] = 0

    heuristic_duplicated_dict["minimizers"][(min_key,min_vigemin)]+=1

t = time.time()-t
print('\nDone in %f s., hence %f s/kmer.' % (t,t/(seq_size-k+1)))

print('Heuristic number of buckets:',len(heuristic_dict["minimizers"].keys()))
print('Heuristic number of buckets (with duplication):',len(heuristic_duplicated_dict["minimizers"].keys()))

heuristic_dict["buckets"] = sorted(heuristic_dict["minimizers"].values(),reverse=True)
heuristic_duplicated_dict["buckets"] = sorted(heuristic_duplicated_dict["minimizers"].values(),reverse=True)

heuristic_dict['density'] = len(heuristic_dict["pos"]) / (seq_size-k+1)
heuristic_duplicated_dict['density'] = len(heuristic_duplicated_dict["pos"]) / (seq_size-k+1)

print('Random minimizer density :', 2/(k-m+1))
print('Density of the heuristic (duplicated or not):',heuristic_dict["density"])

###########################################################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

###########################################################################

fig, ax = plt.subplots(figsize=(12, 8))

for key in keys:
    ax.plot(vigemin_dict[key]["buckets"],color='C0',alpha=0.3)

ax.plot([],[],color='C0',label=str(N_keys)+' keys')

ax.plot(heuristic_dict["buckets"],color='C1',label='Heuristic')
ax.plot(heuristic_duplicated_dict["buckets"],color='C2',label='Heuristic (duplicated)')

ax.set_yscale("log",base=10)
ax.set_xlabel("Bucket rank",fontsize=fontsize)
ax.set_ylabel("Bucket size",fontsize=fontsize)

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.yaxis.get_offset_text().set_fontsize(fontsize)
ax.legend(fontsize=fontsize)

# plt.show()
fig.savefig('pikmin_benchmark_k=' + str(k) + '_m=' + str(m) +'_N_keys=' + str(N_keys) + 'seq_size=' + str(seq_size) + '.pdf',bbox_inches="tight")