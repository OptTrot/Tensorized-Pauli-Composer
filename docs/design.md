## XZ code

$$(n_x, n_z)$$

### Pauli Algebra

$$P_i \cdot P_j = P_k$$

* $P_i = (n_{x}^i, n_{z}^i)$
* $P_j = (n_{x}^j, n_{z}^j)$
* $P_k = (n_{x}^k, n_{z}^k)$

$$n_{x}^k = n_{x}^i {}^\wedge n_{x}^j\\ n_{z}^k = n_{z}^i {}^\wedge n_{z}^j$$

where, ${}^\wedge$ is a XOR binary operator.

### Tensor Product

$$P_i \otimes P_j = P_k$$

* $P_i = (n_{x}^i, n_{z}^i)$
* $P_j = (n_{x}^j, n_{z}^j)$
* $P_k = (n_{x}^k, n_{z}^k)$

$$n_{x}^k = (2^j \cdot n_{x}^i)\vee n_{x}^j$$
$$n_{z}^k = (2^j \cdot n_{z}^i)\vee n_{z}^j$$

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
In the original paper, the chasing the corresponding Pauli-term was achieved by adding characters to string variable.
However, XZ code, `nx, nz`, uniquely determines a corresponding ij index on coefficient matrix by single bit opeartor, where, `nr=nz` and `nc=nx^nz`. 
Reverse transformation is also defined as `nz=nr` and `nx = nr^nc`.

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