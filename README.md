# Minimizer Counting Function

Provided an alphabet $\Sigma$, two integers $m\leq k$ and a word $w\in\Sigma^m$, we define the quantity

$$\pi_k(w) = |\lbrace x\in\Sigma^k : \min(x) = w\rbrace|$$

where $\min(x)$ designates the minimizer of $x$ for the lexicographical order, i.e. the smallest substring of size $m$ found in $x$.

The purpose of this module is to compute $\pi_k(w)$, following the principles defined in the following paper:

> _On the number of $k$-mers admitting a given lexicographical minimizer_, F. Ingels, C. Marchet, M. Salson, Univ. Lille.
> [ref to arxiv]

### Dependencies

The library itself has no dependencies except on `cmath` which is natively included in Python.

To run tests, install `nose`.

To reproduce the figures of the paper and the associated analysis, make sure to have the following Python modules installed: `matplotlib`, `scipy`, `numpy`.

### Usage of the module

TODO