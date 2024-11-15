import sys
import pickle
import time
from lib import *
# import os
# print('Number of CPUs : {}'.format(os.cpu_count()))

#
# minimizer = 'ACA'
# m=len(minimizer)
# k=10
#
# greater_letters_dic = number_of_greater_letters()
#
# obj = PreProcess(minimizer,number_of_greater_letters_dic=greater_letters_dic)
#
# relmat = obj.autocorrelation_matrix
# minj = obj.prefix_letters_vectors
# a_max = obj.a_max
# prefix_max_size = obj.antemer_max_prefix_size
# postmer_max = obj.postmer_max_size
#
# print(postmer_max)

# print(obj.antemer_upper_bound(k))
# print(obj.antemer(k))
# print(obj.antemer_lower_bound(k))
# print('\n')
# print(obj.postmer_upper_bound(k+m)[m:])
# print(obj.postmer(k+m)[m:])
# print(obj.postmer_lower_bound(k+m)[m:])
# print('\n')


#
# for k in range(11):
#     print(k+m,obj.kmer_lower_bound(k+m),obj.kmer(k+m),obj.kmer_upper_bound(k+m))

# for k in range(11):
#     print(k+m, upper_bound(k + m, minimizer, postmer_max, prefix_max_size, greater_letters_dic, minj, a_max), number_of_kmers(k + m, minimizer, postmer_max, prefix_max_size, greater_letters_dic, minj, a_max), lower_bound(k + m, minimizer, postmer_max, prefix_max_size, greater_letters_dic, minj, a_max))

# print('\n')
# print('ANTEMERS')
#
# for alpha in range(6):
#
#     candidates = []
#     for i in range(4**alpha):
#         kmer = int_to_kmer(i,alpha)
#         candidates.append(kmer+minimizer)
#
#     # print(len(candidates))
#
#     acceptables = []
#
#     for kmer in candidates:
#         min_index = find_minimizer(kmer,m)
#         found_minimizer = kmer[min_index:min_index+m]
#         if found_minimizer==minimizer and min_index==alpha:
#             acceptables.append(kmer)
#
#     print(alpha, antemer(alpha, minimizer, prefix_max_size),len(acceptables))
#
# print('\n')
# print('POSTMER')
#
# for beta in range(m,9):
#
#     candidates = []
#     for i in range(4**(beta-m)):
#         kmer = int_to_kmer(i,beta-m)
#         candidates.append(minimizer+kmer)
#
#     # print(len(candidates))
#
#     acceptables = []
#
#     for kmer in candidates:
#         min_index = find_minimizer(kmer,m)
#         found_minimizer = kmer[min_index:min_index+m]
#         if found_minimizer>=minimizer:
#             acceptables.append(kmer)
#
#     print(beta, postmer(beta, minimizer,m),len(acceptables))



k = 31
m = 8

greater_letters_dic = number_of_greater_letters()

minimizer_dic = {}

bound_up = {}
bound_low = {}

tight_bound =0

somme = 0

t = time.time()
for i in range(4**m):
    sys.stdout.write("\rProgression: %i/%i" % (i+1,4**m))
    minimizer = int_to_kmer(i,m)

    obj = PreProcess(minimizer,number_of_greater_letters_dic=greater_letters_dic)

    bound_up[minimizer] = obj.kmer_upper_bound(k)
    bound_low[minimizer] = obj.kmer_lower_bound(k)

    if bound_up[minimizer]==bound_low[minimizer]:
        N= bound_up[minimizer]
        tight_bound+=1
    else:
        N = obj.kmer(k)

    somme+=N

    minimizer_dic[minimizer]= N

    # assert bound_up[minimizer] >= N and N >= bound_low[minimizer]

assert somme == 4**k

t = time.time()-t

print('\nDone in %s seconds' % t)

with open('../Data/data_k='+str(k)+'_m='+str(m)+'.p', 'wb') as f:
    pickle.dump(minimizer_dic, f)

complexity = 4**m * ((2*k-m)*m+m*m)

print('Estimated constant : %s seconds' % (t/complexity))

print('Encountered %i tight bounds, representing %.0f%% of all minimizers' % (tight_bound, 100*tight_bound/4**m))

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.exporters
import math

names = sorted(minimizer_dic.keys())
# names = sorted(bound_up.keys())
# names = sorted(bound_low.keys())

frequencies = []
up_bound = []
low_bound = []

normalization_factor = (k-m) + math.log(k-m+1,4)
# normalization_factor = 1

for name in names:
    frequencies.append(math.log(minimizer_dic[name],4)/normalization_factor)
    up_bound.append(math.log(bound_up[name],4)/normalization_factor)
    low_bound.append(math.log(bound_low[name], 4)/normalization_factor)

x_coordinates = [i/4**m for i in range(len(names))]

win = pg.plot()
# win.setWindowTitle()
win.setBackground('w')

pen = pg.mkPen(color=(0, 0, 0))
pen_red = pg.mkPen(color=(255, 0, 0))
pen_blue = pg.mkPen(color=(0, 0, 255))

data = pg.PlotDataItem(x_coordinates, frequencies,pen=pen)
up_bound = pg.PlotDataItem(x_coordinates,up_bound,pen=pen_blue)
low_bound = pg.PlotDataItem(x_coordinates,low_bound,pen=pen_red)
hline = pg.InfiniteLine(pos=(k - m)/normalization_factor, angle=0, pen=pen_red)

win.addItem(data)
win.addItem(up_bound)
win.addItem(low_bound)
win.addItem(hline)

exporter = pg.exporters.ImageExporter(win.plotItem)
exporter.export('../Images/data_k='+str(k)+'_m='+str(m)+'.png')

QtGui.QApplication.instance().exec_()




