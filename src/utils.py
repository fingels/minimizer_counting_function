from random import randint

def xor_char(a,b):
    '''Implement XOR between two letters from the alphabet A,T,C,G'''
    x = min(a,b)
    y = max(a,b)

    if x==y:
        return 'A'
    elif x=='A':
        return y
    elif (x,y)==('C','G'):
        return 'T'
    elif (x,y)==('C','T'):
        return 'G'
    elif (x,y)==('G','T'):
        return 'C'
    else:
        return ' '

def xor_word(x,y):
    '''Compute the word  = (x XOR y)'''

    assert len(x)==len(y), str(x)+'|'+str(y)
    z = ''

    for i in range(len(x)):
        z+= xor_char(x[i],y[i])

    return z

def int_to_bin(num):
    '''Convert an int to a binary array'''
    binary = []
    while num != 0:
        bit = num % 2
        binary.append(bit)
        num = num // 2
    binary.reverse()
    return binary

def int_to_kmer(num,k):
    '''Convert an int to the corresponding k-mer'''
    binary = int_to_bin(num)
    while len(binary)<2*k:
        binary.insert(0,0)

    dic = {(0, 0): "A", (0, 1): "C", (1, 0): "G", (1, 1): "T"}

    s=''
    for i in range(k):
        bits = (binary[2*i],binary[2*i+1])
        s+=dic[bits]

    assert len(s)==k
    return s


def generate_random_k_mers(k,n=1):
    '''
    Generates n random k-mers
    '''
    nmax = 2**(2*k)-1
    if n==1:
        return int_to_kmer(randint(0,nmax),k)
    else:
        list =  []
        for i in range(n):
            list.append(int_to_kmer(randint(0,nmax),k))
        return list

def find_vigemin(string:str, key:str,m:int):
    '''
    Returns the index where the smallest substring of size `m`, xored against `key`, of `string` begins.
    '''
    assert len(string) >= m, 'String too short'
    assert len(key)==m, 'Key should be the same length as m'

    min_vigemer = 'T' * m
    index_min = float('Inf')

    for i in range(len(string)-m+1):
        if xor_word(string[i:i+m],key) < min_vigemer:
            min_vigemer = xor_word(string[i:i+m],key)
            index_min = i

    return index_min

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

def dna_reverse(letter):
    dic = {'A':'T', 'C':'G','G':'C', 'T':'A'}
    return dic[letter]

def number_of_greater_letters(alphabet = {'A', 'T', 'C', 'G'}):
    dic = {}
    for a in alphabet:
        dic[a] = len([s for s in alphabet if s>a])
    dic[' ']=len(alphabet)
    return dic


def number_of_greater_words(string:str,alphabet = {'A', 'T', 'C', 'G'}):
    # greater or equal
    greater_letters_dic = number_of_greater_letters(alphabet)

    m  = len(string)
    sum = 1 # pour le préfixe=m
    for i in range(0,m):
        ai = string[i]
        prod = greater_letters_dic[ai]*len(alphabet)**(m-i-1)
        sum+= prod
    return sum