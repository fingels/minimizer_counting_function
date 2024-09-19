use std::{collections::HashMap,usize};
use substring::Substring;
use std::cmp;

pub fn autocorrelation_matrix(s:&String,alphabet:&Vec<char>) {

    let m = s.len();

    let mut antemer_max_prefix_size = m;
    let mut postmer_max_size = usize::MAX;
    let mut minj = HashMap::new();
    let mut mat = vec![];

    for i in 0..m {
        let mut row = vec![' '; i+1];
        row[0] = '=';
        mat.push(row);
        
        let mut minj_val = HashMap::new();
        for a in alphabet.iter() {
            if i==0 {
                if a==&s.chars().nth(0).unwrap() {
                    minj_val.insert(a,2);
                }
                else {
                    minj_val.insert(a,0);
                }
            }
            else {
                minj_val.insert(a,usize::MAX);
            }
        }
        minj.insert(i+1, minj_val);
    }

    for j in 1..m {
        let mut flag = false;
        let mut symb = ' ';

        for i in j..m {
            if !flag {

                let cand = s.substring(j,i+1).to_string();
                let sref = s.substring(0, i-j+1).to_string();
                if cand==sref {
                    symb = '=';
                    let a = &s.chars().nth(i-j+1).unwrap();
                    minj.get_mut(&(i+1)).unwrap().entry(&a).and_modify(|elt| *elt = cmp::min(*elt,j+1));

                    // minj[i + 1][string[i - j + 1]] =  min(minj[i + 1][string[i - j + 1]], j + 1)
                }
                else if cand < sref {
                    symb = '<';
                    flag = true;
                    antemer_max_prefix_size = cmp::min(antemer_max_prefix_size,i+1);
                    postmer_max_size = cmp::min(postmer_max_size,j+1);
                }
                else {
                    symb = '>';
                    flag = true;
                }
            }

            mat[i][j] = symb;

            for a in alphabet.iter(){
//             if minj[j+1][a]==float('Inf'):
//                 minj[j+1][a] = (j+2) * (a==string[0])
            }

        }


        //     # print('\n')
//     # for j in range(m):
//     #     print(j+1,mat[j])
//     #

//     return mat, minj, antemer_max_prefix_size,postmer_max_size
    }

    // printing functions

    println!("{postmer_max_size}, {antemer_max_prefix_size}");

    for row in mat.iter(){
        println!("{:?}", row);
    }

    for (i,value) in minj.into_iter() {
        println!("{i} : {:?}",value);
    }


}




    
// def antemer_lower_bound(alpha, minimizer, prefix_max_size, greater_letters_dic, minj, alphabet = {'A', 'T', 'C', 'G'}):
//     array = [0]*(alpha+1)

//     a_max = {}
//     for i in range(1,prefix_max_size):
//         a_max[i]=''
//         for s in alphabet:
//             if minj[i][s] != 0:
//                 a_max[i] = max(s, a_max[i])

//     array[0] = 1

//     for j in range(1,alpha+1):

//         array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

//         for i in range(1,prefix_max_size):

//             alphabet_partition_A = []
//             alphabet_partition_B = []

//             for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
//                 if minj[i][s] == 0:
//                     alphabet_partition_A.append(s)
//                 else:
//                     alphabet_partition_B.append(s)

//             if j-(i+1)>=0:
//                 array[j] += len(alphabet_partition_A) * array[j-(i+1)]

//     # print('B',array)

//     return array

// def antemer_upper_bound(alpha, minimizer, prefix_max_size, greater_letters_dic, minj, alphabet = {'A', 'T', 'C', 'G'}):
//     array = [0]*(alpha+1)

//     a_max = {}
//     for i in range(1,prefix_max_size):
//         a_max[i]=''
//         for s in alphabet:
//             if minj[i][s] != 0:
//                 a_max[i] = max(s, a_max[i])

//     array[0] = 1

//     for j in range(1,alpha+1):

//         array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

//         for i in range(1,prefix_max_size):

//             alphabet_partition_A = []
//             alphabet_partition_B = []

//             for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
//                 if minj[i][s] == 0:
//                     alphabet_partition_A.append(s)
//                 else:
//                     alphabet_partition_B.append(s)

//             if j-(i+1)>=0:
//                 array[j] += len(alphabet_partition_A) * array[j-(i+1)]

//             for a in alphabet_partition_B:
//                 if j-minj[i][a]>=0:
//                     array[j] += array[j - minj[i][a] + 1] - greater_letters_dic[minimizer[0]] * array[j - minj[i][a]]
//                 elif j-minj[i][a]+1>=0:
//                     array[j] += array[j - minj[i][a] + 1]

//     # print('B',array)

//     return array

// def antemer(alpha,minimizer,prefix_max_size,greater_letters_dic,minj,alphabet = { 'A', 'T', 'C', 'G'}):
//     m = len(minimizer)

//     array_prefix = []
//     array = [0]*(alpha+1)

//     for i in range(prefix_max_size):
//         array_prefix.append([0]*(alpha+1))

//     a_max = {}
//     for i in range(1,prefix_max_size):
//         a_max[i]=''
//         for s in alphabet:
//             if minj[i][s] != 0:
//                 a_max[i] = max(s, a_max[i])

//     for j in range(alpha+1):
//         for i in range(prefix_max_size):

//             if i>j:
//                 array_prefix[i][j]=0
//             elif j==0:
//                 array_prefix[i][j]=1
//             elif i==0:
//                 array_prefix[i][j]=greater_letters_dic[minimizer[0]] * array[j-1]
//             elif i==j:
//                 prod = 1
//                 for l in range(i):
//                     # print(prefix,j+1,relmat[prefix-1][j],'|',m,prefix-j+1,relmat[m-1][prefix-j])
//                     prod *= (relmat[i - 1][l] == ">") or (
//                                 relmat[i - 1][l] == '=' and relmat[m - 1][i -l] == '<')
//                 array_prefix[i][j]=prod
//             else:

//                 alphabet_partition_A = []
//                 alphabet_partition_B = []
//                 for s in [a for a  in alphabet if a > minimizer[i] and a >= a_max[i]]:
//                     if minj[i][s] == 0:
//                         alphabet_partition_A.append(s)
//                     else:
//                         alphabet_partition_B.append(s)


//                 assert len(alphabet_partition_B)<=1

//                 A = len(alphabet_partition_A) * array[j-(i+1)]
//                 B = 0
//                 for a in alphabet_partition_B:
//                     for new_prefix in range(i - minj[i][a] + 2, prefix_max_size):
//                         B += array_prefix[new_prefix][j - minj[i][a] + 1]

//                 array_prefix[i][j]= A + B

//         array[j]= sum([array_prefix[i][j] for i in range(prefix_max_size)])

//     # ##############################################
//     # #
//     # for i in range(prefix_max_size):
//     #     print(i,array_prefix[i])
//     #
//     # print('S',array)

//     return array

// def postmer_lower_bound(beta, minimizer, greater_letters_dic, minj, alphabet={'A', 'T', 'C', 'G'}):
//     m = len(minimizer)

//     array = [0]*(beta+1)
//     array_prefix = [0]*(beta+1)

//     for j in range(m):
//         array[j]= len(alphabet)**j

//     array[m]=number_of_greater_words(minimizer)
//     array_prefix[m] =1

//     for j in range(m+1,beta+1):

//         array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

//         for i in range(1,m+1):
//             tilde_minj = {}

//             a_max = ''
//             for a in alphabet:
//                 tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
//                 if tilde_minj[a] != 0:
//                     a_max= max(a_max,a)

//             if i == m:
//                 aip =''
//             else:
//                 aip = minimizer[i]

//             alphabet_partition_A = []
//             alphabet_partition_B = []
//             for s in [a for a in alphabet if a > aip and a >= a_max]:
//                 if tilde_minj[s] == 0:
//                     alphabet_partition_A.append(s)
//                 else:
//                     alphabet_partition_B.append(s)

//             array[j] += len(alphabet_partition_A) * array[j-(i+1)]

//             if i==m:
//                 array_prefix[j] += len(alphabet_partition_A) * array[j-(m+1)]

//     # print('B',array)

//     # return array
//     return array_prefix

// def postmer_upper_bound(beta, minimizer, greater_letters_dic, minj, alphabet={'A', 'T', 'C', 'G'}):
//     m = len(minimizer)

//     array = [0]*(beta+1)
//     array_prefix = [0]*(beta+1)

//     for j in range(m):
//         array[j]= len(alphabet)**j

//     array[m]=number_of_greater_words(minimizer)
//     array_prefix[m] = 1

//     for j in range(m+1,beta+1):

//         array[j] = greater_letters_dic[minimizer[0]] * array[j-1]

//         for i in range(1,m+1):
//             tilde_minj = {}

//             a_max = ''
//             for a in alphabet:
//                 tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
//                 if tilde_minj[a] != 0:
//                     a_max= max(a_max,a)

//             if i == m:
//                 aip =''
//             else:
//                 aip = minimizer[i]

//             alphabet_partition_A = []
//             alphabet_partition_B = []
//             for s in [a for a in alphabet if a > aip and a >= a_max]:
//                 if tilde_minj[s] == 0:
//                     alphabet_partition_A.append(s)
//                 else:
//                     alphabet_partition_B.append(s)

//             array[j] += len(alphabet_partition_A) * array[j-(i+1)]



//             for a in alphabet_partition_B:
//                 array[j] += array[j-tilde_minj[a]+1] - greater_letters_dic[minimizer[0]] * array[j - tilde_minj[a]]

//             if i==m:
//                 array_prefix[j] += len(alphabet_partition_A) * array[j-(m+1)]
//                 for a in alphabet_partition_B:
//                     array_prefix[j] += array[j-tilde_minj[a]+1] - greater_letters_dic[minimizer[0]] * array[j - tilde_minj[a]]

//     # print('B',array)

//     # return array
//     return array_prefix

// def postmer(beta,minimizer,greater_letters_dic,minj,alphabet={ 'A', 'T', 'C', 'G'}):
//     m = len(minimizer)
//     array_prefix = []
//     array = [0]*(beta+1)
//     for i in range(m+1):
//         array_prefix.append([0]*(beta+1))

//     for j in range(beta+1):
//         for i in range(m+1):
//             if i>j:
//                 array_prefix[i][j]=0
//             elif j==0:
//                 array_prefix[i][j]=1
//             elif j<m:
//                 if i==0:
//                     array_prefix[i][j]= (len(alphabet)-1) * array[j-1]
//                 elif i==j:
//                     array_prefix[i][j]=1
//                 else:

//                     alphabet_partition_A = []
//                     alphabet_partition_B = []
//                     for s in [a for a in alphabet if a != minimizer[i]]:
//                         if minj[i][s] == 0:
//                             alphabet_partition_A.append(s)
//                         else:
//                             alphabet_partition_B.append(s)

//                     A = len(alphabet_partition_A) * array[j - (i+1)]
//                     B = 0
//                     for a in alphabet_partition_B:
//                         for new_prefix in range(i - minj[i][a] + 2, j - minj[i][a] + 2):
//                             B += array_prefix[new_prefix][j-minj[i][a]+1]

//                     array_prefix[i][j]= A + B

//             elif j==m:
//                 if i==m:
//                     array_prefix[i][j]=1
//                 else:
//                     array_prefix[i][j]= greater_letters_dic[minimizer[i]] * len(alphabet) ** (j-(i+1))

//             else:
//                 if i==0:
//                     array_prefix[i][j] = greater_letters_dic[minimizer[0]] * array[j-1]
//                 else:

//                     tilde_minj = {}

//                     a_max = ''
//                     for a in alphabet:
//                         tilde_minj[a] = minj[i][a] * (j >= m + minj[i][a] - 1)
//                         if tilde_minj[a] != 0:
//                             a_max= max(a_max,a)

//                     if i == m:
//                         aip =''
//                     else:
//                         aip = minimizer[i]

//                     alphabet_partition_A = []
//                     alphabet_partition_B = []
//                     for s in [a for a in alphabet if a > aip and a >= a_max]:
//                         if tilde_minj[s] == 0:
//                             alphabet_partition_A.append(s)
//                         else:
//                             alphabet_partition_B.append(s)

//                     assert len(alphabet_partition_B) <= 1

//                     A = len(alphabet_partition_A) * array[j - (i+1)]
//                     B = 0
//                     for a in alphabet_partition_B:
//                         for new_prefix in range(i - tilde_minj[a] + 2, m + 1):
//                             B += array_prefix[new_prefix][j-tilde_minj[a]+1]

//                     array_prefix[i][j]= A + B

//         array[j] = sum([array_prefix[i][j] for i in range(m+1)])

//     ##############################################
//     #
//     # for i in range(m+1):
//     #     print(i,array_prefix[i])
//     #
//     # print('S',array)

//     return array_prefix[m]
//     # return array

// def number_of_kmers(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
//     m = len(minimizer)

//     beta_max = min(postmer_max-2,k-m)

//     antemer_array = antemer(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)

//     postmer_array = postmer(beta_max+m,minimizer,greater_letters_dic,minj)

//     somme = 0

//     for beta in range(0,beta_max+1):
//         # somme += antemer(k-m-beta,minimizer,prefix_max_size) * postmer(beta+m,minimizer,m)

//         somme += antemer_array[k - m - beta] * postmer_array[beta+m]

//     return somme

// def upper_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
//     m = len(minimizer)

//     beta_max = min(postmer_max-2,k-m)

//     antemer_array = antemer_upper_bound(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)
//     postmer_array = postmer_upper_bound(beta_max+m,minimizer,greater_letters_dic,minj)

//     somme = 0

//     for beta in range(0,beta_max+1):
//         somme += antemer_array[k-m-beta] * postmer_array[beta+m]

//     return min((beta_max+1)*4**(k-m),somme)

// def lower_bound(k,minimizer,postmer_max,prefix_max_size,greater_letters_dic,minj):
//     m = len(minimizer)

//     beta_max = min(postmer_max-2,k-m)

//     antemer_array = antemer_lower_bound(k-m,minimizer,prefix_max_size,greater_letters_dic,minj)
//     postmer_array = postmer_lower_bound(beta_max+m,minimizer,greater_letters_dic,minj)

//     somme = 0

//     for beta in range(0,beta_max+1):
//         somme += antemer_array[k-m-beta] * postmer_array[beta+m]

//     return max(1,somme)
