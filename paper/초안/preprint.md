

# Fast matrix reconstruction of general combination of Clifford algebra elements


Cliffod algebra now becomes dominant approach in various fields
to analyze complicated systems such as computer graphics, and robotics.
In physics, it has been a fundamental tool reprsent the algebra of 
quantum quantities and has a tremedous applications in mathematical physics.
<!-- Quantum field operator, computer graphics, robotices example references -->
With Clifford/Geometric algebra approach, 
the equations contain more information about geometry and the system
that did not appeared in convenience linear algebra system.
Even computationally hard problems also formulated by Cliiford equations.
<!-- Clifford algebra 로 NP hard 문제를 만든 논문이 있다. -->

However, in high dimension, calculation of Clifford algebra system 
requires significantly large resources in both time and space.
It is caused by the its higher dimension inherited in the algebra.

Even with the highly optimized framework,
matrix based method could show better time cost than the geometric based algorithm.
For example, in the case that rapidly applying same 
transformation for tremendous objects, the matrix based method 
has a benenfit. 
However, the initializing of the matrices was a bottelneck 
of the algorithm.

Therefore, practical calculation of the clifford algebra requires 
switching the representation by the situation to acclerate the overall
speed of the algorithm.
In the paper, we constructed an inverse algorithm to
generate a marix representation of general elements of Clifford algebra,
based on generalized Pauli element schema.
In general case, matrix reconstruction routines 
were *single element reconstruction* algorithm.
Meanwhile, in here, the method is for general elements.
For example, with quatornion system, $\mathbb{H} \approx Cl_3$.

$$M = a + bi + cj + dk \rightarrow  
\begin{bmatrix} 
a + b i & c+ di\\
-c + di & a- bi
\end{bmatrix}
$$

This is a method to construct a tensor product to the all elements 
in the given 

The paper consist of 3 parts. The first is reviewing a tensorized 
decomposition algorithm of hermit matrix into generalized Pauli matrices.
In second part, we will construct a coefficient matrix with only Pauli element index
in symplectic form. 
Lastly, we will show that the inverse algorithm could be used to
construct general Clifford algebra elements. 

## Matrix representation of Clifford algebra

## Hamiltonian decomposition problem

### Decomposition algorithm

### Symplectic representation of Pauli elements

### Inverse algorithm

The reason why we cannot use an inverse tensorized algorithm
is that we have no idea that 
where the given generalized Pauli element would be located 
in the coefficient matrix.
Since, inverse direction also has been ambiguous so that 
the original tensorized algorithm chase the each term in
the algorithm process, so that the required memory was higher 
than the written complexity.


Lemma: 


With the above theorem, we can directly extract the corresponding Pauli element
from the coefficient index and since the bitwise XOR operator is an inverse itself.
The coefficient index from Pauli element is also well defined.
With the conversion we can directly build a coefficient matrix 
from Pauli polynomial and reconstruct the overall matrix at once.

Naive

## Application to represent general Clifford algebra element

### CL(0, 2)

