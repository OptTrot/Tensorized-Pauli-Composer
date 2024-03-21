## XZ code

$$(n_x, n_z)$$

### Pauli Algebra

$$P_i \cdot P_j = P_k$$

* $P_i = (n_{x}^i, n_{z}^i)$
* $P_j = (n_{x}^j, n_{z}^j)$
* $P_k = (n_{x}^k, n_{z}^k)$

$$n_{x}^k = n_{x}^i {}^\wedge n_{x}^j\\ n_{z}^k = n_{z}^i {}^\wedge n_{z}^j$$

where, ${}^\wedge$ is a XOR binary operator.

### Commutator

See, Reggio et al.

**Thorem**
For two given Pauli string, $P_i$ and $P_j$ and their xz-code representation. 

$$Com(P_i, P_j) = (Com(z_i, x_j) = Com(z_j, x_i))$$

where, $Com(P_i, P_j) = True$ if $[P_i, P_j] =0$ else $False$.

With binary representation of integer, $Com$ function could be implemented as $\&$ operator and "1" value count.

$$Count(n_x \& n_z, 1)$$

### Tensorized Decomposition

See Hantzko et al, Tensorized Pauli decomposition algorithm.

The algorithm decompose the given Hamiltonian only use iterative submatrix addtion and subtraction. The final result is a coefficient matrix whose location indicates corresponding Pauli string.

With xz-code, `nx, nz`, the index of the corresponding Pauli term is determined with single bit opeartor, where, `nr=nz` and `nc=nx^nz`. 
Reverse transformation is defined as `nz=nr` and `nx = nr^nc`.

Note that, the property of XOR operator is

```
b =  a^(a^b)
a =  b^(a^b)
```


## C-Python interfaces

* ctypes
* Cython
* cffi

There are some utils in numpy library.

https://numpy.org/doc/stable/reference/routines.ctypeslib.html#module-numpy.ctypeslib