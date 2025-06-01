from typing import *
from itertools import product
from functools import reduce
import numpy as np
from numpy import kron

#-------------------
FLT_EPS = 1E-8
#--------------------

# Basic Pauli matrices
I = np.eye(2, dtype=complex)
X = np.matrix([[0, 1], [1, 0]], dtype=complex)
Y = np.matrix([[0, 1], [-1, 0]], dtype=complex) # omit -1j* phase
Z = np.matrix([[1, 0],[0, -1]], dtype = complex)
# Pauli matrix by name
PAULI_MATRICES = {"I": I,"X": X,"Y": Y,"Z": Z}
# SIMPLECTIC representation
PAULI_SIMPLECTIC = {"I":(0,0), "X":(1,0), "Y":(1,1), "Z":(0,1)}
SIMPLECTIC_PAULI={(0,0): "I",(1,0): "X",(1,1): "Y",(0,1): "Z"}

def krons(*oper_list): # Operator Kronecker delta
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)
def frobenius_inner(A, B): # Frobenius inner product.
    n, n2 = A.shape
    return np.trace((A.conj().T)@B)/(n)
def mat_diff(A, B):
     D = A-B
     return np.sqrt(frobenius_inner(D, D))
def pstr2matrix(pstr:str)->np.ndarray:
    """Convert Pauli string to corresponding matrix representation.
    Args:
        pstr (str): Pauli-string. For example, "IXYZ" is a Pauli-string of length 4.
    Returns:
        np.Ndarray: Corresponding matrix, in the example, I (x) X (x) Y (x) Z is returned. (x): <., .> is a kronecker delta.
    """
    result = []
    for p in pstr:
        result.append(PAULI_MATRICES[p])
    return krons(result)

def ppoly2mat(ppoly):
    pstr, _ = ppoly[0]
    N = int(2**len(pstr))
    mat = np.zeros((N, N)).astype(complex)
    for pstr, w in ppoly:
          mat += w*pstr2matrix(pstr)
    return mat
          

def pstr2xz_code(pstr:str)->Tuple[int, int]:
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
        nx, nz = PAULI_SIMPLECTIC[p]
        x_num += nx*num
        z_num += nz*num
        num += num
    return x_num, z_num
# Pauli string to matrix - Naive
def pstr2mat(pstr:str)->np.matrix:
        result = []
        for p in pstr:
            result.append(PAULI_MATRICES[p])
        phase = (-1j)**(pstr.count("Y")%4)
        return phase*krons(result)

def pstr2sym_code(pstr:str, sim_code:Union[dict, None]=None)->Tuple[int,int]:
        if sim_code is None:
            global PAULI_SIMPLECTIC
            pauli_sim_dict = PAULI_SIMPLECTIC
        else:
            pauli_sim_dict = sim_code
        num = 1

        x_num = 0 
        z_num = 0

        # x,z_num = 1*2^0 + 0*2^1 + 1*2^2 + ... 
        for p in reversed(pstr):
            nx, nz = pauli_sim_dict[p]
            x_num += nx*num
            z_num += nz*num
            num += num # 2*num
        return (x_num, z_num)

def pstr2ij_code(pstr:str):
     return sym_code2ij_code(*pstr2sym_code(pstr))
def sym_code2ij_code(x:int, z:int)->Tuple[int, int]:
    """(nx, nz) -> (i, j)

    Args:
        x (int): symplectic code of x
        z (int): symplectic code of z

    Returns:
        Tuple[int, int]: _description_
    """
    return  z, x^z
def ij_code2sym_code(i:int, j:int)->Tuple[int, int]:
    """(i, j) -> (x, z)

    Args:
        i (int): row index of canonical matrix
        j (int): column index of canonical matrix

    Returns:
        Tuple[int, int]: _description_
    """
    return i^j, i
def sym_code2pstr(ns:Tuple[int, int], l:int)->str:
        assert l>0, "l must be positive integer and greater than 0."
        nx, nz = ns
        max_int_1 = 2**l
        assert (nx < max_int_1 and nz < max_int_1), "The given integers and the qubit dim are not matched."
        if nx==0: # Z family
            st = format(nz, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "Z")
            return st
        if nz==0: # X family
            st = format(nx, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "X")
            return st
        # None of above
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
def ij_code2_pstr(ns:Tuple[int, int], l:int)->str:
     return sym_code2pstr(ij_code2sym_code(*ns), l)
# General Pauli terms
def get_pstrs(n:int):
     return list(map(lambda x: "".join(x), product(f"IXYZ", repeat=int(n))))
def pstrs2mats(pstrs:list[str]):
     return [pstr2mat(p) for p in pstrs]
def get_pauli_fam_terms(n, fam="Z"):
        return list(map(lambda x: "".join(x), product(f"I{fam}", repeat=int(n))))
def get_pauli_fam_mat(n, fam="Z"):
        return list(map(krons, product([I, PAULI_MATRICES[fam]], repeat=int(n))))
def pstr_commute(pa, pb):
    nx1, nz1 =  pstr2sym_code(pa)
    nx2, nz2 =  pstr2sym_code(pb)

    a = bin(nx1&nz2).count("1")%2
    b = bin(nx2&nz1).count("1")%2
    return a==b
