# Opttrot

Implmented routines about optimizing trotterization circuit of the given Hamiltonian, $H$.

* Basic Pauli-group algebra based on bit-string written in C.
* Fast Pauli-polynomial generation of the given Hamiltonian.
* Accelerated Hamiltonian decomposition.
* Pauli code and matrix bijective transformation.

## Z2xZ2 representation of Pauli-group.

$n$-length Pauli-string, for example, `IXXZY` is a Pauli-string of length 5,
is represented with $(n_x, n_z) \in \mathbb{Z}_2 \times \mathbb{Z}_2$.
This representation is based on the fast commuting determination algorithm suggested by Reggio et al.
They suggested the algorithm based on the 2 integer tuple, but such representation also has a great potential to implement fast algebric operations about Pauli-groups.

There is a limitation of the such integer tuple notation. In IEEE764 bit representation, 64bit integer consist of 63 bits to manipulated by integer operation, the first bit is used for integer sign.
Therefore, for $2^n$ bit system, we have capcacity for $2^{n -1}$. In 2024, common systems are 64bit and 32bit, however, $2^{31}$

* Integer tuple to Pauli-string
* Pauli-string to tuple
* Pauli-term synthesis
* General Hamiltonian to coefficient, xz-code and reverse.

## Interface for general optimizer 

Thereare many optimizer and solving algorithms for classic max-clique and 
quantum optimizer such as D-Wave and QuERA. 
*Opttrot* provides general interface to those optimizers, 
Currently, Networks, and D-Wave optimizers could be used in Opttrot.
