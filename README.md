# Tensorized Pauli Composer

This repository is an implementation of Pauli composition routine, Tensorized Pauli Composer,
to construct a Hamiltonian matrix defined by weighted Pauli polynomial.
The method is a combination of canonical matrix construction method based on symplectic representation of Pauli element,
and inversed Tensorized Pauli Decomposition([TPD](https://iopscience.iop.org/article/10.1088/1402-4896/ad6499)) method.

- `tp`: Compact implementation and utils.
- `data`: Benchmark data and source codes, initial implementation of OptTrot library.
- `preprint`: Tech document and theory.


### Speed comparsion

A(xN.N): B took N.N times more slower than A.
In qubits < 12

|Opeartion|Faster|
|:-------:|:----:|
|o->mat   | PauliPoly(x2-1.2)| 
|mat->o   | Qiskit(x1.2)| 
|Addition | PauliPoly(x2) |
|Mat Mul  | PauliPoly(x Exponential) |


Single Pauli term

|Qiskit| TPC | Numpy Corresponding|
|:------:|:-----:|:-------------------:|
|27μs+-161ns| 1.31μs+-3.89ns | 11.7μs +- 102ns|


Pauli poly ( Matrix construction 12 qubits)
|Qiskit| TPC |
|:------:|:-----:|
|33.2 s ± 1.98 s|3.48 s ± 12.3 ms|



To cite the software, please refer [Zenodo](https://zenodo.org/records/14245728) 
preprint.

```
Kim, H. (2024). TENSORIZED PAULI COMPOSER. Quantum Thread 2024, Seoul. Zenodo. https://doi.org/10.5281/zenodo.14245728
```

## TODO

- Diagonal decomposition, composition routines.
