from typing import Tuple
import numpy as np
from scipy.sparse import coo_matrix
from .utils import pstr2ij_code, ij_code2_pstr

def ppoly2pauli_basis(ppoly:list[Tuple[str, float]], sparse=False)->np.matrix:
    pstr, _ = ppoly[0]
    N = int(2**len(pstr))
    mat = np.matrix(np.zeros((N, N)).astype(complex))
    for pstr, w in ppoly:
        i, j = pstr2ij_code(pstr)
        mat[i, j] = w
    return mat if not sparse else coo_matrix(mat)


def pauli_basis2ppoly(cmat:np.matrix):
    n = int(np.log2(cmat.shape[0]))
    s_mat = coo_matrix(cmat)
    ppoly = []
    for i, j, w in zip(s_mat.row, s_mat.col, s_mat.data):
        pstr = ij_code2_pstr((i,j), n)
        ppoly.append((pstr, w))
    return ppoly

#-------------------------------------------------------------------
def tpd(H:np.matrix)->np.matrix:
    """Tensorized matrix construction method from weighted Pauli sum.
    Args:
        H (np.matrix): Pauli_basis matrix of Pauli elements

    Returns:
        np.matrix: Restored Hamiltonian
    """
    #Tensrosized reconstruction method: O(n4^n)
    # Normal method: O(16^n)
    # mat = np.zeros(self.coef_matrix.shape) 
    # for p in self.poly:
    #   mat += p.coef*p.matrix
    #mat = self.coef_matrix
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
                
                # There is no problem of inplace operators.
                
                # I-Z
                H[num_i: num_i+l, num_j:num_j+l]        += H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] 
                H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] = H[num_i: num_i+l, num_j:num_j+l] - 2*H[num_i+l: num_i+2*l, num_j+l:num_j+2*l]
                # X-Y
                H[num_i: num_i+l, num_j+l:num_j+2*l] +=  H[num_i+l: num_i+2*l, num_j:num_j+l] 
                H[num_i+l: num_i+2*l, num_j:num_j+l]  =  H[num_i: num_i+l, num_j+l:num_j+2*l] - 2*H[num_i+l: num_i+2*l, num_j:num_j+l]
                H[num_i+l: num_i+2*l, num_j:num_j+l] *= 1j

    H *= (1/(2**n))
    return H

def tpd_xz(H:np.matrix)->np.matrix:
    """Modified Tensorized decomposition of hermit matrix into pauli terms.
    It generates a pauli_basis matrix whose row and column index is a symplectic tuple of Pauli terms.

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

                row1 = slice(num_i, num_i +l)
                row2 = slice(num_i+l, num_i +2*l)
                col1 = slice(num_j, num_j +l)
                col2 = slice(num_j+l, num_j +2*l)

                # Step 1: Swap the Z X
                H[row1, col2] = H[row1, col2] + H[row2, col2]
                H[row2, col2] = H[row1, col2] - H[row2, col2]
                H[row1, col2] = H[row1, col2] - H[row2, col2]

                # Step 2:
                H[row1, col1] = H[row1, col1] + H[row1, col2]
                H[row2, col1] = H[row2, col1] + H[row2, col2]

                # Step 3:
                H[row1, col2] = H[row1, col1] - 2* H[row1, col2]
                H[row2, col2] = H[row2, col1] - 2* H[row2, col2]

                # Step 4: Phase of Y
                H[row2, col2] = -1j*H[row2, col2]

    H *= (1/(2**n))
    return H

def itpd(ppoly:list[Tuple[str, float]])->np.matrix:
    return itpd_core(ppoly2pauli_basis(ppoly))

def itpd_core(mat :np.matrix)->np.matrix:
    """Tensorized matrix construction method from weighted Pauli sum.
    Args:
        mat (np.matrix): pauli_basis matrix of Pauli elements

    Returns:
        np.matrix: Restored Hamiltonian
    """
    _2n = mat.shape[0] # 2^n
    steps = int(np.log2(_2n))# n
    unit_size = 1
 
    for step in range(steps):
        step1 = step+1
        mat_size = int(2*(unit_size))
        indexes = np.arange(_2n/(2**step1)).astype(np.uint)
        indexes_ij = mat_size * indexes
        for i in indexes_ij:
            for j in indexes_ij:
                # (i, j)
                r1i     = i
                c1i     = j
                r1f2i   = int(r1i   + (unit_size))
                c1f2i   = int(c1i   + (unit_size))
                r2f     = int(r1f2i + (unit_size))
                c2f     = int(c1f2i + (unit_size))

                # Do not replace the below code to in-place operator += or *=.
                # Numba jit yieds different handling process in compile time. 
                # I - Z
                coef = 1
                mat[r1i: r1f2i, c1i:c1f2i] = mat[r1i: r1f2i, c1i:c1f2i] + coef*mat[r1f2i: r2f, c1f2i:c2f]
                mat[r1f2i: r2f, c1f2i:c2f] = mat[r1i: r1f2i, c1i:c1f2i] -2*coef *mat[r1f2i: r2f, c1f2i:c2f]
                # X -Y
                coef = -1j
                mat[r1i: r1f2i, c1f2i:c2f] = mat[r1i: r1f2i, c1f2i:c2f]  + coef*mat[r1f2i: r2f, c1i:c1f2i]
                mat[r1f2i: r2f, c1i:c1f2i] = mat[r1i: r1f2i, c1f2i:c2f] - 2*coef*mat[r1f2i: r2f, c1i:c1f2i]
                
        unit_size =int(2*unit_size)
    return mat

def itpd_eff(ppoly):
    pauli_basis = ppoly2pauli_basis(ppoly, sparse=True)
    return itpd_eff_core(pauli_basis.toarray(), np.stack([pauli_basis.row, pauli_basis.col]).T)


def itpd_eff_core(
    mat:np.matrix, 
    p_indexes:Tuple[Tuple[int, int], ...], 
    ):
    """Conversion routine from a coefficient matrix to original matrix, in computational basis 

    Args:
        mat (np.Matrix): A pauli_basis matrix of the given Pauli polynomial.
        p_indexes (Tuple[Tuple[int, int], ...]): Non zero term indexes.

    Returns:
        _type_: _description_
    """
    _2n = mat.shape[0] # 2^n
    steps = int(np.log2(_2n))# n
    unit_size= 1
    
    # Effective term chasing routine.
    in_ps = [p for p in p_indexes]

    for step in range(steps):

        psteps =[]
        dup = []
        #----
        for p in in_ps:
            i, j = p
            if (i, j) in dup: continue
            #p_class = (i+j)%2 #0: IZ, 1:XY, It was used in `coef` determination.
            n, o = i%2, j%2 # (1), (2) determination
            l, m = (i+1-2*(n), j+ 1-(2*(o))) # get a corresponding location
            
            dup.append((i,j))
            dup.append((l,m))
            
            #------------------
            # The blow commnets were replaced upper code "if (i, j) in dup: continue"
            #------------------
            #if (l,m) in dup: # Eliminate duplicated operation. 
            #    continue
            #elif (i, j) in dup: # At first, 
            #    dup.append((l,m))
            #    continue
            #else:
            #    dup.append((i,j))
            #    dup.append((l,m))
            if n: # (2)
                pair = ((l, m), (i, j))
            else: #(1)
                pair = ((i, j), (l, m))

            coef = -1j if (i+j)%2 else 1 # ture: XY, false: IZ
            
            r1i = unit_size*pair[0][0]
            r1f = r1i + unit_size
            c1i = unit_size*pair[0][1]
            c1f = c1i +unit_size
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