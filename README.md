# Tensorized Pauli Composer

This repository is an implementation of Pauli composition routine, Tensorized Pauli Composer,
to construct a Hamiltonian matrix defined by weighted Pauli polynomial.
The method is a combination of canonical matrix construction method based on symplectic representation of Pauli element,
and inversed Tensorized Pauli Decomposition([TPD](https://iopscience.iop.org/article/10.1088/1402-4896/ad6499)) method.

- `tp`: Compact implementation and utils.
- `data`: Benchmark data and source codes, initial implementation of OptTrot library.
- `preprint`: Tech document and theory.

To cite the software, please refer [Zenodo](https://zenodo.org/records/14245728) 
preprint.

```
Kim, H. (2024). TENSORIZED PAULI COMPOSER. Quantum Thread 2024, Seoul. Zenodo. https://doi.org/10.5281/zenodo.14245728
```

## TODO

- Diagonal decompositio, composition routines.
