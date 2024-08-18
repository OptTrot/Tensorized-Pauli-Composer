import numpy as np
from typing import Tuple
from numba import jit


class TensorizedMethods:
    @staticmethod
    def mat2coefs_naive(H:np.matrix)-> np.matrix:
        """Tensorized decomposition of hermit matrix into pauli terms.
            See Hantzko et al, 2023.

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
                    H[num_i+l: num_i+2*l, num_j:num_j+l] *= 1j

        H *= (1/(2**n))
        return H
    @staticmethod
    def coefs2mat_naive(coefs_mat:np.matrix)->np.matrix:
        """Tensorized composition routine of coefficient matrix.

        Args:
            coefs_mat (np.matrix): Coefficient matrix whose element indicates each Pauli term weight.

        Returns:
            np.matrix: Matrix.

        """
        mat = coefs_mat
        _2n = mat.shape[0] # 2^n
        steps = int(np.log2(_2n))# n
        unit_size= 1
 
        for step in range(steps):
            step1 = step+1
            mat_size = int(2*(unit_size))
            indexes = np.arange(_2n/(2**step1)).astype(np.uint)
            indexes_ij = mat_size * indexes
            for i in indexes_ij:
                for j in indexes_ij:
                    # (i, j)
                    r1i     = i
                    r1f2i   = r1i + unit_size
                    c1i     = j
                    c1f2i   = c1i + +unit_size
                    r2f     = r1f2i + unit_size
                    c2f     = c1f2i + unit_size
                    # I - Z
                    coef = 1
                    mat[r1i: r1f2i, c1i:c1f2i] += coef*mat[r1f2i: r2f, c1f2i:c2f]
                    mat[r1f2i: r2f, c1f2i:c2f] = mat[r1i: r1f2i, c1i:c1f2i] -2*coef *mat[r1f2i: r2f, c1f2i:c2f]
                    # X -Y
                    coef = -1j
                    mat[r1f2i: r2f, c1i:c1f2i] += coef*mat[r1i: r1f2i, c1f2i:c2f]
                    mat[r1i: r1f2i, c1f2i:c2f] = mat[r1f2i: r2f, c1i:c1f2i] -2*coef *mat[r1i: r1f2i, c1f2i:c2f]
            
            unit_size *=2
        return mat
    @staticmethod
    def mat2coefs_np(H:np.matrix)-> np.matrix:
        """Tensorized decomposition of hermit matrix into pauli terms.
            See Hantzko et al, 2023.

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
                    H[num_i+l: num_i+2*l, num_j:num_j+l] *= 1j

        H *= (1/(2**n))
        return H