# Minimizer Counting Function

Provided an alphabet $\Sigma$, two integers $m\leq k$ and a word $w\in\Sigma^m$, we define the quantity

$$\pi_k(w) = |\lbrace x\in\Sigma^k : \min(x) = w\rbrace|$$

where $\min(x)$ designates the minimizer of $x$ for the lexicographical order, i.e. the smallest substring of size $m$ found in $x$.

The purpose of this module is to compute $\pi_k(w)$, following the principles defined in the following paper:

[REF TO ADD]

### Dependencies

Make sure to have the following python module installed:

`nose`, `matplotlib`

### Usage of the module

TODO