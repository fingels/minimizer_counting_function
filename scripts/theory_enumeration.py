from src.lib import *
import time
import sys
import pickle

#####################################

k = 31
m = 10
compute_bounds = False

#####################################

greater_letters_dic = number_of_greater_letters()

minimizer_dic = {}

bound_up = {}
bound_low = {}

equal_bounds =0
tight_bounds = 0

somme = 0

t = time.time()
for i in range(4**m):
    sys.stdout.write("\rProgression: %i/%i" % (i+1,4**m))
    minimizer = int_to_kmer(i,m)

    obj = MinimizerCountingFunction(minimizer, number_of_greater_letters_dic=greater_letters_dic)

    if compute_bounds:
        bound_up[minimizer] = obj.kmer_upper_bound(k)
        bound_low[minimizer] = obj.kmer_lower_bound(k)

        if bound_up[minimizer]==bound_low[minimizer]:
            N= bound_up[minimizer]
            equal_bounds+=1
        else:
            N = obj.kmer(k)

            if N==bound_up[minimizer] or N==bound_low[minimizer]:
                tight_bounds+=1
    else:
        N = obj.kmer(k)

    somme+=N

    minimizer_dic[minimizer]= N

    if compute_bounds:
        assert bound_up[minimizer] >= N >= bound_low[minimizer], minimizer

assert somme == 4**k

t = time.time()-t

print('\nDone in %s seconds' % t)

if compute_bounds:
    print('Encountered %i tight bounds, representing %.0f%% of all minimizers; among which they were %i equal bounds, for %.0f%% of all minimizers.' % (
        equal_bounds+tight_bounds, 100 * (equal_bounds+tight_bounds) / 4 ** m,equal_bounds, 100*equal_bounds/4**m))
#
# with open('../Data/enumeration_k='+str(k)+'_m='+str(m)+'.p', 'wb') as f:
#     pickle.dump(minimizer_dic, f)
#
# if compute_bounds:
#     bounds = {'up':bound_up, 'low':bound_low}
#     with open('../Data/bounds_k=' + str(k) + '_m=' + str(m) + '.p', 'wb') as f:
#         pickle.dump(bounds, f)

