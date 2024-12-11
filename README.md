# Minimizer Counting Function

Provided an alphabet $\Sigma$, two integers $m\leq k$ and a word $w\in\Sigma^m$, we define the quantity

$$\pi_k(w) = |\lbrace x\in\Sigma^k : \min(x) = w\rbrace|$$

where $\min(x)$ designates the minimizer of $x$ for the lexicographical order, i.e. the smallest substring of size $m$ found in $x$.

The purpose of this module is to compute $\pi_k(w)$, following the principles defined in the following paper:

> _On the number of $k$-mers admitting a given lexicographical minimizer_, F. Ingels, C. Marchet, M. Salson, Univ. Lille.
> [ref to arxiv]

### Dependencies

The library itself has no dependencies except on `cmath` and `typing` which are natively included in Python.

To run tests, use `nose`.

To reproduce the figures of the paper and the associated analysis, make sure to have the following Python modules installed: `matplotlib`, `scipy`, `numpy`.

### Usage of the module

```
from src.lib import MinimizerCountingFunction

minimizer = 'ACACAA'
k = 10

obj = MinimizerCountingFunction(minimizer)

print(obj.kmer_lower_bound(k))
print(obj.kmer(k))
print(obj.kmer_upper_bound(k))

>>> 327
>>> 351
>>> 353
```

By default, the alphabet is `{ 'A', 'C', 'T','G'}`.  If you want to use another alphabet, use 
```
obj = MinimizerCountingFunction(minimizer,alphabet=my_custom_alphabet)
```
Make sure that the minimizer is written in said alphabet. `my_custom_alphabet` should be a set of characters.

If you want to use a custom ordering of the letters instead of the classical one `A < B < ...`, define a dictionary where each key is a letter of the alphabet, and the value is the number of letters strictly greater that the key in the alphabet. For practical reasons, the dictionary must also have `my_dict[' ']=len(alphabet)`. 
Then use 
```
obj = MinimizerCountingFunction(minimizer,number_of_greater_letters_dic=my_dict)
```
