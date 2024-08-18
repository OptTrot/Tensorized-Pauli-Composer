
from __future__ import annotations

from numbers import Number
from typing import *
from functools import reduce, product
#--------------------------------
import numpy as np


## Naive Pauli utils
class PauliUtils:
    I = np.eye(2)
    pauli_X = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_Y = complex(0, 1)*np.array([[0, -1], [1, 0]], dtype=complex)
    pauli_Z = np.array([[1, 0], [0, -1]], dtype=complex)
    p_basis = {"I":I, "X":pauli_X, "Y":pauli_Y, "Z":pauli_Z}
    p_map = {"I":(0,0), "X":(1, 0), "Y":(1,1), "Z":(0,1)}
    
    def __init__(self):
        pass
    # Basic utils
    def krons(*oper_list): # Operator Kronecker delta
        if len(oper_list) == 1:
            oper_list = oper_list[0]
        return reduce(np.kron, oper_list)
    def frobenius_inner(A, B): # Frobenius inner product.
        n, n2 = A.shape
        return np.trace((A.conj().T)@B)/(n)
    
    def pstr_to_matrix(pstr:str)->np.ndarray:
        """Convert Pauli string to corresponding matrix representation.

        Args:
            pstr (str): Pauli-string. For example, "IXYZ" is a Pauli-string of length 4.

        Returns:
            np.Ndarray: Corresponding matrix, in the example, I (x) X (x) Y (x) Z is returned. (x): <., .> is a kronecker delta.
        """
        result = []
        for p in pstr:
            result.append(PauliUtils.p_basis[p])
        return PauliUtils.krons(result)
    
    def pstr_to_xz_code(pstr:str)->Tuple[int, int]:
        """Convert Pauli string to xz family code.

        Args:
            pstr (str): Pauli string

        Returns:
            Tuple[int, int]: XZ family
        """
        num = 1
        x_num = 0 # Consider a bit represenation
        z_num = 0

        for p in reversed(pstr):
            nx, nz = PauliUtils.p_map[p]
            x_num += nx*num
            z_num += nz*num
            num += num
        return x_num, z_num
    def xz_code_to_pstr(ns:Tuple[int, int], l:int)->str:
        """Convert XZ family code to a corresponding Pauli string.

        Args:
            ns (Tuple[int, int]): XZ family of Pauli term. 
            l (int): Length of Pauli string. 

        Returns:
            str: Pauli string of length `l`.
        """
        assert l>0, "l must be positive integer and greater than 0."
        nx, nz = ns
        max_int_1 = 2**l
        assert (nx < max_int_1 and nz < max_int_1), "The given integers and the qubit dim are not matched."
        if nx==0:
            st = format(nz, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "Z")
            return st
        if nz==0:
            st = format(nx, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "X")
            return st

        st_x = format(nx, f"0{l}b")
        st_z = format(nz, f"0{l}b")
        result = []
        for x, z in zip(st_x, st_z):
            if x == z:
                if x =="1":
                    result.append("Y")
                else: 
                    result.append("I")
            elif x > z:
                result.append("X")
            else:
                result.append("Z")
        return "".join(result)

    def get_pauli_family_string(n, fam="Z"):
        return list(map(lambda x: "".join(x), product(f"I{fam}", repeat=int(n))))

    def get_pauli_family_matrix(n:int, fam="Z")->Iterable[np.matrix]:
        """Get pauli_family of `n` qubits of `fam` family. 

        Args:
            n (int): Number of qubits. The output matrices are :math:`2^n`.
            fam (str, optional): Type of Pauli-family of X, Y, or Z. Defaults to "Z".

        Returns:
            Iterable[np.matrix]: list of Pauli-matrices
        """
        G = PauliUtils.pauli_Z if fam=="Z" else (PauliUtils.pauli_X if fam=="X" else PauliUtils.pauli_Y)

        return list(map(PauliUtils.krons, product([PauliUtils.I, G], repeat=int(n))))

class Pauli:
    """For the scalable version, see Qiskit Pauli class.
    """
    def __init__(self, x:int, z:int, n:int, coef=1):
        assert isinstance(x, int), "x must be integer."
        assert isinstance(z, int), "z must be integer."
        assert isinstance(n, int), "z must be integer."
        assert x>0, "x must be positive integer."
        assert z>0, "z must be positive integer."
        assert n>0, "n must be positive integer."
        
        max_n = int(2**n)

        assert x<max_n, "x must be smaller than {}".format(max_n)
        assert z<max_n, "z must be smaller than {}".format(max_n)

        self.x = x
        self.z = z
        self.n = n
        self.string = PauliUtils.xz_code_to_pstr((self.x, self.z), self.n)
        self.coef = coef
    def __eq__(self, other:Pauli):
        if self.x == other.x and self.z == other.z:
            return True
        else:
            return False 
    def __str__(self):
        return f"{self.string}"
    def __repr__(self):
        return f"Pauli(xz=({self.x}, {self.z}), string={self.string}, n={self.n})"
    def __rmul__(self, other:Number):
        return PauliPoly([Pauli(self.x, self.z, self.n, self.coef * other)])
    def __add__(self, other:Pauli):
        if self == other:
            return Pauli(self.x, self.z, self.n, self.coef + other.coef)
        return PauliPoly([self, other])
    @staticmethod
    def ij_to_xz(i, j):
        return i^j, i 
    @property
    def coef_ij(self):
        return self.z, self.x^self.z

    def __matmul__(self, other:Pauli):
        n = self.n
        new_x = self.x^other.x
        new_z = self.z^other.z
        return Pauli(new_x, new_z, n)
    def otimes(self, other:Pauli):
        """Tensor product of P1 (x) P2

        Args:
            other (Pauli): Pauli class
        """
        n = self.n + other.n
        coef_m = 2**other.n
        new_x = (coef_m*self.x)|other.x
        new_z = (coef_m*self.z)|other.z
        return Pauli(new_x, new_z, n)
    def commute(self, other):
        n1 = bin(self.x&other.z).count("1")%2
        n2 = bin(self.z&other.x).count("1")%2
        return n1==n2
    def to_matrix(self) -> np.ndarray:
        return PauliUtils.pstr_to_matrix(self.string)

class PauliPoly:
    def __init__(self, p_list:Iterable[Pauli], coefs:Union[None, Iterable[Number]]=None):
        if isinstance(coefs, Iterable) and not isinstance(coefs, str):
            for i, coef in enumerate(coefs):
                p_list[i].coef = coef
        self._terms = p_list
    def __str__(self):
        return f"Pauli polynomial of {self._terms[0].n} qubit space."
    def __rmul__(self, other:Number):
        for i in range(len(self._terms)):
            self._terms[i].coef *= other
    def __add__(self, other:Union[Pauli, PauliPoly]):
        if isinstance(other, Pauli):
            p = PauliPoly([other])
        else:
            p = other
        mat = self.coef_matrix + p.coef_matrix
        return PauliPoly.from_coef_mat(mat)
    @classmethod
    def from_coef_mat(cls, mat, qubits, tol=1E-16):
        n, m = mat.shape
        p_list = []
        for i in range(n):
            for j in range(m):
                if np.fabs(mat[i, j]) <tol: 
                    continue
                x, z = Pauli.ij_to_xz(i, j)
                p_list.append(Pauli(x, z, qubits, mat[i, j]))
        return cls(p_list)
    @classmethod
    def from_hermit(cls, H):
        n_u2 = H.shape[0]
        n = int(np.log2(n_u2))
        coef_mat = PauliPoly.hermit_to_coef_mat(H)
        return cls(coef_mat, n)
    @property
    def terms(self):
        return [p.string for p in self._terms]
    @property
    def poly(self):
        return self._terms
    @property
    def coef_matrix(self):
        n = self._terms[0].n
        nn = 2**n
        mat = np.zeros((nn, nn), dtype=complex)
        for p in self._terms:
            mat[*p.coef_ij] = p.coef
        return mat
    def to_matrix(self)->np.matrix:
        """Coefficient matrix to Hermit matrix, inverse algorithm of Tensorized decomposition.

        Returns:
            np.matrix: Hermit matrix
        """
        mat = self.coef_matrix
        p_indexes = [(p.x, p.z) for p in self.poly]
        steps = int(np.log2(mat.shape[0]))
        unit_size= 1
        in_ps = [p for p in p_indexes]
        for step in range(steps):
            #print("========================================")
            psteps =[]
            dup = []
            #----
            for p in in_ps:
                i, j = p
                p_class = (i+j)%2 #0: IZ, 1:XY
                n, o = i%2, j%2 # (1), (2) determination
                l, m = (i+1-2*(n), j+ 1-(2*(o))) # get a corresponding location

                if (l,m) in dup: # Eliminate duplicated operation.
                    continue
                elif (i, j) in dup:
                    dup.append((l,m))
                    continue
                else:
                    dup.append((i,j))
                    dup.append((l,m))

                if n: # (2)
                    pair = ((l, m), (i, j))
                else: #(1)
                    pair = ((i, j), (l, m))

                #print((i,j), (l,m))
                #print(dup)
                #print("XY" if p_class else "IZ")
                #print("(2)" if n else "(1)")
                #print("size:", unit_size)
                #print(f"{unit_size*pair[0][0]}:{unit_size*pair[0][0]+unit_size},  {unit_size*pair[0][1]}:{unit_size*pair[0][1]+unit_size}")
                #print(f"{unit_size*pair[1][0]}: {unit_size*pair[1][0]+unit_size}, {unit_size*pair[1][1]}:{unit_size*pair[1][1]+unit_size}")

                coef = -1j if p_class else 1 # ture: XY, false: IZ

                r1i = unit_size*pair[0][0]
                r1f = r1i + unit_size
                c1i = unit_size*pair[0][1]
                c1f = c1i + +unit_size

                r2i = unit_size*pair[1][0]
                r2f = r2i + unit_size
                c2i = unit_size*pair[1][1]
                c2f = c2i + unit_size

                mat[r1i: r1f, c1i:c1f] += coef*mat[r2i: r2f, c2i:c2f]
                mat[r2i: r2f, c2i:c2f] = mat[r1i: r1f, c1i:c1f] -2*coef *mat[r2i: r2f, c2i:c2f]

                i >>=1
                j >>=1
                if (i, j) in psteps:
                    continue
                else:
                    psteps.append((i,j))
            #----
            in_ps = [p for p in psteps]
            unit_size *=2
        return mat
    @staticmethod
    def hermit_to_coef_mat(H):
        """Tensorized decomposition of hermit matrix into pauli terms.

        Args:
            H (np.matrix): Hermit matrix.

        Returns:
            np.matrix: Coefficient matrix of the given matrix.
        """
        n1, n2 = H.shape
        assert n1 == n2, "The given matrix must be a square matrix."
        n= int(np.log2(n1))
        l = n1
        for i in range(n):
            m = int(2**i) # Number of submatrix
            l = int(l/2) # Sub matrix size, square
            for j in range(m):
                for k in range(m):
                    num_i = j*(2*l) # Initial position of sub matrix row
                    num_j = k*(2*l) # Initial position of sub matrix column

                    # I-Z
                    H[num_i: num_i+l, num_j:num_j+l]        += H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] 
                    H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] = H[num_i: num_i+l, num_j:num_j+l] - 2*H[num_i+l: num_i+2*l, num_j+l:num_j+2*l]
                    # X-Y
                    H[num_i: num_i+l, num_j+l:num_j+2*l] += H[num_i+l: num_i+2*l, num_j:num_j+l] 
                    H[num_i+l: num_i+2*l, num_j:num_j+l] =  H[num_i: num_i+l, num_j+l:num_j+2*l] - 2*H[num_i+l: num_i+2*l, num_j:num_j+l]
                    H[num_i+l: num_i+2*l, num_j:num_j+l] *= -1j

        H *= (1/(2**n))
        return H