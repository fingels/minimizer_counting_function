from random import randint
from functools import lru_cache
import sys
import pickle
import time
import numpy as np

def int_to_bin(num):
    binary = []
    while num != 0:
        bit = num % 2
        binary.append(bit)
        num = num // 2
    binary.reverse()
    return binary

def kmer_to_bin(kmer):
    dic = {'A':[0,0], "C":[0,1], "G":[1,0], "T":[1,1]}

    bin = []
    for char in kmer:
        bin+=dic[char]

    return bin

def int_to_kmer(num,k):
    binary = int_to_bin(num)
    while len(binary)<2*k:
        binary.insert(0,0)

    dic = {(0, 0): "A", (0, 1): "C", (1, 0): "G", (1, 1): "T"}

    s=''
    for i in range(k):
        bits = tuple(binary[2*i:2*i+2])
        s+=dic[bits]

    assert len(s)==k
    return s

def generate_k_mers(k,n=1):
    nmax = 2**(2*k)-1
    if n==1:
        return int_to_kmer(randint(0,nmax),k)
    else:
        list =  []
        for i in range(n):
            list.append(int_to_kmer(randint(0,nmax),k))

def find_minimizer(string: str, m: int):
    '''
    Returns the index where the smallest (w.r.t. the lexicographical order) substring of size `m` of `string` begins.
    '''
    assert len(string) >= m, 'String too short'

    index_min = 0
    for i in range(1, len(string) - m + 1):
        for j in range(m):  # on parcoure string[i:i+m] tant que les lettres coîncident avec string[index_min:index_min+m]
            if string[i + j] < string[index_min + j]:
                index_min = i
                break
            elif string[index_min + j] < string[i + j]:
                break

    return index_min

def number_of_greater_letters(alphabet = { 'A', 'T', 'C', 'G'}):
    dic = {}
    for a in alphabet:
        dic[a] = len([s for s in alphabet if s>a])
    return dic

def number_of_greater_words(string:str,alphabet = { 'A', 'T', 'C', 'G'}):
    # greater or equal
    greater_letters_dic = number_of_greater_letters(alphabet)

    m  = len(string)
    sum = 1 # pour le préfixe=m
    for i in range(0,m):
        ai = string[i]
        prod = greater_letters_dic[ai]*len(alphabet)**(m-i-1)
        sum+= prod
    return sum

def relation_matrix(string:str,alphabet = { 'A', 'T', 'C', 'G'}):
    m = len(string)

    antemer_max_prefix_size = m
    postmer_max_size = float('Inf')
    minj = {}

    mat = []
    for i in range(m):
        mat.append(["="]+[None]*i)

        minj[i+1] = {a:float('Inf') for a in alphabet}

    minj[1] = {a:0 for a in alphabet}
    minj[1][string[0]]=2

    for j in range(1,m):
        flag = False
        symb = None

        for i in range(j,m):
            if not flag:
                cand = string[j:i + 1]
                ref = string[:i - j + 1]
                if cand == ref:
                    symb = "="
                    minj[i + 1][string[i - j + 1]] = min(minj[i + 1][string[i - j + 1]], j + 1)
                elif cand < ref:
                    symb = "<"
                    flag = True
                    antemer_max_prefix_size = min(antemer_max_prefix_size,i+1)
                    postmer_max_size = min(postmer_max_size, j + 1)
                else:
                    symb = '>'
                    flag = True

            mat[i][j]=symb

        for a in alphabet:
            if minj[j+1][a]==float('Inf'):
                minj[j+1][a] = (j+2) * (a==string[0])

    # print('\n')
    # for j in range(m):
    #     print(j+1,mat[j])
    #

    return mat, minj, antemer_max_prefix_size,postmer_max_size

# def minj_dic(relmat,string:str,alphabet = { 'A', 'T', 'C', 'G'}):
#
#     m=len(string)
#
#     dic = {}
#
#     for i in range(m):
#
#         dic[i+1] = { a:float('Inf') for a in alphabet}
#
#         for j in range(1,i+1):
#             if relmat[i][j]=="=":
#                 dic[i+1][string[i-j+1]]=min(dic[i+1][string[i-j+1]],j+1)
#                 break
#
#         for a in alphabet:
#             if dic[i+1][a]==float('Inf'):
#                 if a==string[0]:
#                     dic[i+1][a]=i+2
#                 else:
#                     dic[i+1][a]=0
#
#     return dic

def antemer_lower_bound(alpha, minimizer, prefix_max_size, greater_letters_dic, minj, alphabet = {'A', 'T', 'C', 'G'}):
    array = [0]*(alpha+1)

    a_max = {}
    for i in range(1,prefix_max_size):
        a_max[i]=''
        for s in alphabet:
            if minj[i][s] != 0:
                a_max[i] = max(s, a_max[i])

    array[0] = 1

    for j in range(1,alpha+1):

        array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

        for i in range(1,prefix_max_size):

            alphabet_partition_A = []
            alphabet_partition_B = []

            for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
                if minj[i][s] == 0:
                    alphabet_partition_A.append(s)
                else:
                    alphabet_partition_B.append(s)

            if j-(i+1)>=0:
                array[j] += len(alphabet_partition_A) * array[j-(i+1)]

    # print('B',array)

    return array

def antemer_upper_bound(alpha, minimizer, prefix_max_size, greater_letters_dic, minj, alphabet = {'A', 'T', 'C', 'G'}):
    array = [0]*(alpha+1)

    a_max = {}
    for i in range(1,prefix_max_size):
        a_max[i]=''
        for s in alphabet:
            if minj[i][s] != 0:
                a_max[i] = max(s, a_max[i])

    array[0] = 1

    for j in range(1,alpha+1):

        array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

        for i in range(1,prefix_max_size):

            alphabet_partition_A = []
            alphabet_partition_B = []

            for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
                if minj[i][s] == 0:
                    alphabet_partition_A.append(s)
                else:
                    alphabet_partition_B.append(s)

            if j-(i+1)>=0:
                array[j] += len(alphabet_partition_A) * array[j-(i+1)]

            for a in alphabet_partition_B:
                if j-minj[i][a]>=0:
                    array[j] += array[j - minj[i][a] + 1] - greater_letters_dic[minimizer[0]] * array[j - minj[i][a]]
                elif j-minj[i][a]+1>=0:
                    array[j] += array[j - minj[i][a] + 1]

    # print('B',array)

    return array

def antemer(alpha,minimizer,prefix_max_size,greater_letters_dic,minj,alphabet = { 'A', 'T', 'C', 'G'}):
    m = len(minimizer)

    array_prefix = []
    array = [0]*(alpha+1)

    for i in range(prefix_max_size):
        array_prefix.append([0]*(alpha+1))

    a_max = {}
    for i in range(1,prefix_max_size):
        a_max[i]=''
        for s in alphabet:
            if minj[i][s] != 0:
                a_max[i] = max(s, a_max[i])

    for j in range(alpha+1):
        for i in range(prefix_max_size):

            if i>j:
                array_prefix[i][j]=0
            elif j==0:
                array_prefix[i][j]=1
            elif i==0:
                array_prefix[i][j]=greater_letters_dic[minimizer[0]] * array[j-1]
            elif i==j:
                prod = 1
                for l in range(i):
                    # print(prefix,j+1,relmat[prefix-1][j],'|',m,prefix-j+1,relmat[m-1][prefix-j])
                    prod *= (relmat[i - 1][l] == ">") or (
                                relmat[i - 1][l] == '=' and relmat[m - 1][i -l] == '<')
                array_prefix[i][j]=prod
            else:

                alphabet_partition_A = []
                alphabet_partition_B = []
                for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
                    if minj[i][s] == 0:
                        alphabet_partition_A.append(s)
                    else:
                        alphabet_partition_B.append(s)


                assert len(alphabet_partition_B)<=1

                A = len(alphabet_partition_A) * array[j-(i+1)]
                B = 0
                for a in alphabet_partition_B:
                    for new_prefix in range(i - minj[i][a] + 2, prefix_max_size):
                        B += array_prefix[new_prefix][j - minj[i][a] + 1]

                array_prefix[i][j]= A + B

        array[j]= sum([array_prefix[i][j] for i in range(prefix_max_size)])

    # ##############################################
    # #
    # for i in range(prefix_max_size):
    #     print(i,array_prefix[i])
    #
    # print('S',array)

    return array

def postmer_lower_bound(beta, minimizer, greater_letters_dic, minj, alphabet={'A', 'T', 'C', 'G'}):
    m = len(minimizer)

    array = [0]*(beta+1)
    array_prefix = [0]*(beta+1)

    for j in range(m):
        array[j]= len(alphabet)**j

    array[m]=number_of_greater_words(minimizer)
    array_prefix[m] =1

    for j in range(m+1,beta+1):

        array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

        for i in range(1,m+1):
            tilde_minj = {}

            a_max = ''
            for a in alphabet:
                tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
                if tilde_minj[a] != 0:
                    a_max= max(a_max,a)

            if i == m:
                aip =''
            else:
                aip = minimizer[i]

            alphabet_partition_A = []
            alphabet_partition_B = []
            for s in [a for a in alphabet if a > aip and a >= a_max]:
                if tilde_minj[s] == 0:
                    alphabet_partition_A.append(s)
                else:
                    alphabet_partition_B.append(s)

            array[j] += len(alphabet_partition_A) * array[j-(i+1)]

            if i==m:
                array_prefix[j] += len(alphabet_partition_A) * array[j-(m+1)]

    # print('B',array)

    # return array
    return array_prefix

def postmer_upper_bound(beta, minimizer, greater_letters_dic, minj, alphabet={'A', 'T', 'C', 'G'}):
    m = len(minimizer)

    array = [0]*(beta+1)
    array_prefix = [0]*(beta+1)

    for j in range(m):
        array[j]= len(alphabet)**j

    array[m]=number_of_greater_words(minimizer)
    array_prefix[m] = 1

    for j in range(m+1,beta+1):

        array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

        for i in range(1,m+1):
            tilde_minj = {}

            a_max = ''
            for a in alphabet:
                tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
                if tilde_minj[a] != 0:
                    a_max= max(a_max,a)

            if i == m:
                aip =''
            else:
                aip = minimizer[i]

            alphabet_partition_A = []
            alphabet_partition_B = []
            for s in [a for a in alphabet if a > aip and a >= a_max]:
                if tilde_minj[s] == 0:
                    alphabet_partition_A.append(s)
                else:
                    alphabet_partition_B.append(s)

            array[j] += len(alphabet_partition_A) * array[j-(i+1)]



            for a in alphabet_partition_B:
                array[j] += array[j-tilde_minj[a]+1] - greater_letters_dic[minimizer[0]] * array[j - tilde_minj[a]]

            if i==m:
                array_prefix[j] += len(alphabet_partition_A) * array[j-(m+1)]
                for a in alphabet_partition_B:
                    array_prefix[j] += array[j-tilde_minj[a]+1] - greater_letters_dic[minimizer[0]] * array[j - tilde_minj[a]]

    # print('B',array)

    # return array
    return array_prefix

def postmer(beta,minimizer,greater_letters_dic,minj,alphabet={ 'A', 'T', 'C', 'G'}):
    m = len(minimizer)
    array_prefix = []
    array = [0]*(beta+1)
    for i in range(m+1):
        array_prefix.append([0]*(beta+1))

    for j in range(beta+1):
        for i in range(m+1):
            if i>j:
                array_prefix[i][j]=0
            elif j==0:
                array_prefix[i][j]=1
            elif j<m:
                if i==0:
                    array_prefix[i][j]= (len(alphabet)-1) * array[j-1]
                elif i==j:
                    array_prefix[i][j]=1
                else:

                    alphabet_partition_A = []
                    alphabet_partition_B = []
                    for s in [a for a in alphabet if a != minimizer[i]]:
                        if minj[i][s] == 0:
                            alphabet_partition_A.append(s)
                        else:
                            alphabet_partition_B.append(s)

                    A = len(alphabet_partition_A) * array[j - (i+1)]
                    B = 0
                    for a in alphabet_partition_B:
                        for new_prefix in range(i - minj[i][a] + 2, j - minj[i][a] + 2):
                            B += array_prefix[new_prefix][j-minj[i][a]+1]

                    array_prefix[i][j]= A + B

            elif j==m:
                if i==m:
                    array_prefix[i][j]=1
                else:
                    array_prefix[i][j]= greater_letters_dic[minimizer[i]] * len(alphabet) ** (j-(i+1))

            else:
                if i==0:
                    array_prefix[i][j] = greater_letters_dic[minimizer[0]] * array[j-1]
                else:

                    tilde_minj = {}

                    a_max = ''
                    for a in alphabet:
                        tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
                        if tilde_minj[a] != 0:
                            a_max= max(a_max,a)

                    if i == m:
                        aip =''
                    else:
                        aip = minimizer[i]

                    alphabet_partition_A = []
                    alphabet_partition_B = []
                    for s in [a for a in alphabet if a > aip and a >= a_max]:
                        if tilde_minj[s] == 0:
                            alphabet_partition_A.append(s)
                        else:
                            alphabet_partition_B.append(s)

                    assert len(alphabet_partition_B) <= 1

                    A = len(alphabet_partition_A) * array[j - (i+1)]
                    B = 0
                    for a in alphabet_partition_B:
                        for new_prefix in range(i - tilde_minj[a] + 2, m + 1):
                            B += array_prefix[new_prefix][j-tilde_minj[a]+1]

                    array_prefix[i][j]= A + B

        array[j] = sum([array_prefix[i][j] for i in range(m+1)])

    ##############################################
    #
    # for i in range(m+1):
    #     print(i,array_prefix[i])
    #
    # print('S',array)

    return array_prefix[m]
    # return array

def number_of_kmers(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
    m = len(minimizer)

    beta_max = min(postmer_max-2,k-m)

    antemer_array = antemer(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)

    postmer_array = postmer(beta_max+m,minimizer,greater_letters_dic,minj)

    somme = 0

    for beta in range(0,beta_max+1):
        # somme += antemer(k-m-beta,minimizer,prefix_max_size) * postmer(beta+m,minimizer,m)

        somme += antemer_array[k - m - beta] * postmer_array[beta+m]

    return somme

def upper_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
    m = len(minimizer)

    beta_max = min(postmer_max-2,k-m)

    antemer_array = antemer_upper_bound(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)
    postmer_array = postmer_upper_bound(beta_max+m,minimizer,greater_letters_dic,minj)

    somme = 0

    for beta in range(0,beta_max+1):
        somme += antemer_array[k-m-beta] * postmer_array[beta+m]

    return min((beta_max+1)*4**(k-m),somme)

def lower_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
    m = len(minimizer)

    beta_max = min(postmer_max-2,k-m)

    antemer_array = antemer_lower_bound(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)
    postmer_array = postmer_lower_bound(beta_max+m,minimizer,greater_letters_dic,minj)

    somme = 0

    for beta in range(0,beta_max+1):
        somme += antemer_array[k-m-beta] * postmer_array[beta+m]

    return max(1,somme)


minimizer = 'ACACAC'
m=len(minimizer)
k=15

greater_letters_dic = number_of_greater_letters()

relmat, minj,prefix_max_size,postmer_max = relation_matrix(minimizer)

print(relmat)
print(minj)
print(prefix_max_size,postmer_max)
print('\n')
print(antemer_upper_bound(k,minimizer,prefix_max_size,greater_letters_dic,minj))
print(antemer(k,minimizer,prefix_max_size,greater_letters_dic,minj))
print(antemer_lower_bound(k,minimizer,prefix_max_size,greater_letters_dic,minj))
print('\n')
print(postmer_upper_bound(k+m,minimizer,greater_letters_dic,minj)[m:])
print(postmer(k+m,minimizer,greater_letters_dic,minj)[m:])
print(postmer_lower_bound(k+m,minimizer,greater_letters_dic,minj)[m:])
print('\n')
print(upper_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj))
print(number_of_kmers(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj))
print(lower_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj))

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

#
# k = 31
# m = 6
#
# greater_letters_dic = number_of_greater_letters()
#
# minimizer_dic = {}
#
# bound_up = {}
# bound_low = {}
#
# tight_bound =0
#
# somme = 0
#
# t = time.time()
# for i in range(4**m):
#     sys.stdout.write("\rProgression: %i/%i" % (i+1,4**m))
#     minimizer = int_to_kmer(i,m)
#
#     relmat, minj,prefix_max_size,postmer_max = relation_matrix(minimizer)
#
#     bound_up[minimizer] = upper_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj)
#     bound_low[minimizer] = lower_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj)
#
#     if bound_up[minimizer]==bound_low[minimizer]:
#         N= bound_up[minimizer]
#         tight_bound+=1
#     else:
#         N = number_of_kmers(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj)
#
#     somme+=N
#
#     minimizer_dic[minimizer]= N
#
#     # assert bound_up[minimizer] >= N and N >= bound_low[minimizer]
#
# assert somme == 4**k
#
# with open('/home/fingels/Documents/dev/kmer_count_vectors/kmer_factorisation/data_k='+str(k)+'_m='+str(m)+'.p', 'wb') as f:
#     pickle.dump(minimizer_dic, f)
#
# t = time.time()-t
#
# print('\nDone in %s seconds' % t)
#
# complexity = 4**m * ((2*k-m)*m+m*m)
#
# print('Estimated constant : %s seconds' % (t/complexity))
#
# print('Encountered %i tight bounds, representing %.0f%% of all minimizers' % (tight_bound, 100*tight_bound/4**m))
#
# import pyqtgraph as pg
# from pyqtgraph.Qt import QtCore, QtGui
# import pyqtgraph.exporters
# import math
#
# names = sorted(minimizer_dic.keys())
# # names = sorted(bound_up.keys())
# # names = sorted(bound_low.keys())
#
# frequencies = []
# up_bound = []
# low_bound = []
#
# for name in names:
#     frequencies.append(math.log(minimizer_dic[name],4))
#     up_bound.append(math.log(bound_up[name],4))
#     low_bound.append(math.log(bound_low[name], 4))
#
# x_coordinates = list(range(len(names)))
#
# win = pg.plot()
# # win.setWindowTitle()
# win.setBackground('w')
#
# pen = pg.mkPen(color=(0, 0, 0))
# pen_red = pg.mkPen(color=(255, 0, 0))
# pen_blue = pg.mkPen(color=(0, 0, 255))
#
# data = pg.PlotDataItem(x_coordinates, frequencies,pen=pen)
# up_bound = pg.PlotDataItem(x_coordinates,up_bound,pen=pen_blue)
# low_bound = pg.PlotDataItem(x_coordinates,low_bound,pen=pen_red)
# hline = pg.InfiniteLine(pos=k - m, angle=0, pen=pen_red)
#
# win.addItem(data)
# win.addItem(up_bound)
# win.addItem(low_bound)
# win.addItem(hline)
#
# exporter = pg.exporters.ImageExporter(win.plotItem)
# exporter.export('/home/fingels/Documents/dev/kmer_count_vectors/kmer_factorisation/data_k='+str(k)+'_m='+str(m)+'.png')
#
# QtGui.QApplication.instance().exec_()




