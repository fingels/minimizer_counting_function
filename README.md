# MCF - Minimizer Counting Function

Provided an alphabet $\Sigma$, two integers $m\leq k$, and two words $w,\gamma\in\Sigma^m$, we define the quantity

$$\pi_k^\gamma(w) = |\lbrace x\in\Sigma^k : \min_\gamma(x) = w\rbrace|$$

where $\min_\gamma(x) = \text{argmin}_{(w'\in\Sigma^m) \wedge (w'\subseteq x)} \lbrace w'\oplus\gamma \rbrace$ (where the min is determined by the lexicographical order).

In other words, $\pi_k^\gamma(w)$ counts the number of $k$-mers (among all $|\Sigma|^k$ possible ones) that admits $w$ as their minimizer according to the order inducing by XORing $m$-mers against $\gamma$.

The purpose of this module is to compute $\pi_k^\gamma(w)$.

## Lexicographical minimizer counting function

The special case of $\gamma = \texttt{A}\cdots\texttt{A}$ - i.e. where $w'\oplus \gamma = w'$ relates to computing the function

$$\pi_k^\gamma(w) = |\lbrace x\in\Sigma^k : \min(x) = w\rbrace|$$

where $\min(x)$ designates the minimizer of $x$ for the lexicographical order, i.e. the smallest substring of size $m$ found in $x$.

The `LexMinimizerCountingFunction` class implements the equations and principles defined in the following paper to compute $\pi_k(w)$:

> _On the number of $k$-mers admitting a given lexicographical minimizer_
> F. Ingels, C. Marchet, M. Salson.
> (https://arxiv.org/abs/2412.17492)

It can be used as follows.

```
from src.lib import LexMinimizerCountingFunction

minimizer = 'ACACAA'
k = 10

obj = LexMinimizerCountingFunction(minimizer)

print(obj.kmer_lower_bound(k))
print(obj.kmer(k))
print(obj.kmer_upper_bound(k))

>>> 327
>>> 351
>>> 353
```

By default, the alphabet is `{'A','C','T','G'}`.  If you want to use another alphabet, use 
```
obj = LexMinimizerCountingFunction(minimizer,alphabet=my_alphabet)
```
Make sure that the minimizer is written in said alphabet. `my_custom_alphabet` should be a set of characters.

If you want to use a custom ordering of the letters instead of the classical one `A < B < ...`, define a dictionary where each key is a letter of the alphabet, and the value is the number of letters strictly greater that the key in the alphabet. For practical reasons, the dictionary must also have `my_dict[' ']=len(alphabet)`. 

As an example, the dictionary associated to the default alphabet is `{'A':3,'C':2,'G':1,'T':0,' ':4}`.

Then use 
```
obj = LexMinimizerCountingFunction(minimizer,number_of_greater_letters_dic=my_dict)
```

You can simultaneously use a custom alphabet and a custom order of the letters of that alphabet by providing both arguments to `MinimizerCountingFunction`. In the case of a custom alphabet only, the dictionary will be automatically computed using standard lexicographical order.

## Vigemin counting function

The `VigeminCountingFunction` class implements the solution to the general problem of computing $\pi_k^\gamma(w)$. It implements the equations and methods of the following paper, that are extensions of the lexicographical one:

> _Vigemers: on the number of $k$-mers sharing the same XOR-based minimizer_
> F. Ingels, A. Limasset, C. Marchet, M. Salson.
> (https://arxiv.org/abs/2602.03337)

It can be used as follows.

```
from src.lib import VigeminCountingFunction

minimizer = 'ACACAA'
key = 'CTGGGT'
k = 10

obj = VigeminCountingFunction(minimizer,key)

print(obj.kmer(k))

>>> 31

```

Unlike the previous class, no lower and upper bounds have been implemented (yet). While `VigeminCountingFunction` can be used with $\gamma = \texttt{A}\cdots\texttt{A}$, we suggest that you instead use `LexMinimizerCountingFunction` instead, as the code is more optimized to deal with the specificity of this case.

Finally, the only supported alphabet (yet) is the DNA alphabet `{'A','C','T','G'}`.


## Dependencies

The library itself has no dependencies except on `cmath` and `typing` which are natively included in Python.

To run tests, use `nose`.

To reproduce the figures of the papers and the associated analysis, make sure to have the following Python modules installed: `matplotlib`, `scipy`, `numpy`.

