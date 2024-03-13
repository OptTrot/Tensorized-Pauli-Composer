from opttrot import *

I = np.eye(2)
pauli_X = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_Y = complex(0, 1)*np.array([[0, -1], [1, 0]], dtype=complex)
pauli_Z = np.array([[1, 0], [0, -1]], dtype=complex)
p_basis = {"I":I, "X":pauli_X, "Y":pauli_Y, "Z":pauli_Z}

# default bit functions
def krons(*oper_list): # Operator Kronecker delta
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)
def frobenius_inner(A, B): # Frobenius inner product.
    n, n2 = A.shape
    return np.trace((A.conj().T)@B)/(n)
#--------------------------------------------------------

def get_decomposition(pauli_basis:dict)->pd.DataFrame:
    """Convert Pauli term and coefficient dictionary to dataframe

    Args:
        pauli_basis (dct): Pauli polynomial of dictionary form. {"Pauli-term": coefficition}

    Returns:
        pandas.DataFrame: [Pstring", "type", "Z", "X", "Coef"]
    """
    p_dict = {}
    for p in pauli_basis.keys():
        nx, nz = pstr_to_xz_fam_code(p)
        num = 1 if nx>0 else 0
        num += num if nz>0 else 0
        p_dict[p] = (num, nz, nx, pauli_basis[p])
    df = pd.DataFrame.from_dict(
                                p_dict, 
                                orient="index",
                                columns = ["type", "Z", "X", "Coef"])
    df.reset_index(inplace=True, names="Pstring")
    return df

def pstr_to_matrix(pstr:str)->np.ndarray:
    """Convert Pauli string to corresponding matrix representation.

    Args:
        pstr (str): Pauli-string. For example, "IXYZ" is a Pauli-string of length 4.

    Returns:
        np.Ndarray: Corresponding matrix, in the example, I (x) X (x) Y (x) Z is returned. (x): <., .> is a kronecker delta.
    """
    result = []
    for p in pstr:
        result.append(p_basis[p])
    return krons(result)

#----------------------------------------------------------------------
# XZ Family code encoding of Pauli term.
# XZ fam code: Pauli term = (n_z:int, n_x:int) 
# This represents location of latin matrix of Pauli-term.
p_map = {"I":(0,0), "X":(1, 0), "Y":(1,1), "Z":(0,1)}

# This allow us to implement Z_3 group strucutre with 2 bit binary
# of each position of IEEE integer value.
# With this representation, we can implement Reggio et al method
# for fast determining of two Pauli term and Pauli term algebra.
# 
# ::For example, 
# "IXYYZ" = (IIZZZ)-(IXXXI)
# (0b00111) - (0b01110)
# 0b(IXYYZ) = (0b-0-0-1-1-1) + (0b0-1-1-1-0-)
# = 0b|00|10|11|11|01
# Addition of two Pauli term in XZ fam code
# p1 = [n1z, n1x]
# p2 = [n2z, n2x]
# p1 + p2 = [n1z^n2z, n1x^n2x], "^" XOR bit operation.


def pstr_to_xz_fam_code(pstr:str)->Tuple[int, int]:
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
        nx, nz = p_map[p]
        x_num += nx*num
        z_num += nz*num
        num += num
    return x_num, z_num

def xz_fam_code_to_pstr(ns:Tuple[int, int], l:int)->str:
    """Convert XZ family code to corresponding Pauli string.

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

def xz_fam_code_add(fam1, fam2):
    z_1, x_1 = fam1
    z_2, x_2 = fam2
    
    return [z_1^z_2, x_1^x_2]
# Bit operators for XZ codes=========================================
# Python standard routines are well optimized 
# It is no worth to use bit optmization
# Just remained as example of further implementation on compile language
# such as C. 
# ------------------------------------------------------
# pauli coef calculation
def pauli_xz_product_coef(x_int, z_int):
    return 1j**(bit_count(x_int&z_int))
def bit_count(n:int):
    #Brian Kernighan’s Algorithm.
    num = 0
    while n:
        n &= n-1
        num+=1
    return num

# Calculate string from integers
int_pchar = ["I", "Z", "X", "Y"]
def pstr_from_xz(x_int, z_int): # Same with `xz_fam_code_to_pstr` function
    z_modi = insert_zeros_in_gaps(z_int)
    x_modi = insert_zeros_in_gaps(x_int)
    x_modi <<= 1

    p_int = x_modi + z_modi 
    # In binary representation: (00)(10)(10)(11) form
    # 00:I, 10: X, 01: Z, 11: Y

    # Get length of str
    len_p = 0
    tem = p_int
    while tem:
        len_p +=1
        tem >>=1
    len_p += len_p&1
    pstr = len_p*['']
    i = 1
    while p_int >0:
        p = p_int & 3
        p_int >>= 2
        pstr[-i] = int_pchar[p]
        i+=1
    return "".join(pstr)

def insert_zeros_in_gaps(n):
    result = 0
    bit_position = 0

    while n > 0:
        # Isolate the rightmost bit
        rightmost_bit = n & 1
        # Shift the bit to its new position
        result |= rightmost_bit << (bit_position << 1)
        # Move to the next bit
        n >>= 1
        bit_position += 1

    return result
def commute_reggio(pa:Tuple[int, int], pb:Tuple[int, int])->bool:
    """Calculate commutation of two Pauli terms encoded in XZ family code.
    The result is a boolean value, `True` and `False`.
    The reference is "Reggio et al, Fast Partitioning of Pauli Strings into Commuting Families for Optimal Expectation Value Measurements of Dense Operators, arXiv, 2023-06-07".

    Args:
        pa (Tuple[int, int]): Pauli term 
        pb (Tuple[int, int]): Pauli term 

    Returns:
        bool: `True` is commute, and `False` is anti-commute.
    """
    nx_a, nz_a = pa
    nx_b, nz_b = pb
    
    a = bin(nx_a & nz_b).count("1")%2
    b = bin(nx_b & nz_a).count("1")%2
    return a==b
    
def commute_reggio_df(s):
    """Dataframe version of `commute_reggio` function.

    Args:
        s (_type_): _description_

    Returns:
        _type_: _description_
    """
    a = bin(s.iloc[0] & s.iloc[3]).count("1")%2
    b = bin(s.iloc[1] & s.iloc[2]).count("1")%2
    return a == b

def integer_order_map(int_list):
    sorted_unique = np.unique(np.array(int_list))
    return {num: idx for idx, num in enumerate(sorted_unique)}


def get_coef(x_str, z_str): 
    # i coefficient in construction of general pauli-element from XZ elements.
    # Use this function in python module
    # The below bitwise implementations are slower than the current function.
    # They are just for further C implementation.
    n = len(x_str)
    x_str = x_str.replace("X", "1")
    x_str = x_str.replace("I", "0")
    z_str = z_str.replace("Z", "1")
    z_str = z_str.replace("I", "0")
    
    x_int = int(x_str, 2)
    z_int = int(z_str, 2)
    return get_coef_bin(x_int, z_int, n)

def get_coef_bin(x_int:int, z_int:int, n:str):
    y_pos = x_int&z_int
    y_pos = format(x_int&z_int, f"0{n}b")
    z_pos = format((x_int|z_int) - x_int, f"0{n}b")
    x_pos = format((x_int|z_int) - z_int, f"0{n}b")

    g_str = []
    for x,y,z in zip(x_pos, y_pos, z_pos):
        if x==y and y==z:
            g_str.append("I")
        elif x== "1":
            g_str.append("X")
        elif y == "1":
            g_str.append("Y")
        else:
            g_str.append("Z")
    return 1j**y_pos.count("1"), "".join(g_str)

# Basis transformation weight by Kim
_k_weight = {
    "I" : {
        "I": 0,
        "Z": 0,
        "X": 1,
        "Y": 2
    },
    "Z" : {
        "I": 0,
        "Z": 0,
        "X": 1,
        "Y": 2
    },
    "X" : {
        "I": 1,
        "Z": 1,
        "X": 0,
        "Y": 3
    },
    "Y" : {
        "I": 2,
        "Z": 2,
        "X": 3,
        "Y": 0
    },
}

def get_basis_weight(s):
    s_str = s.iloc[0]
    t_str = s.iloc[1]
    w = 0
    for s_i, t_i in zip(s_str, t_str):
        w+= _k_weight[s_i][t_i]
    return w/len(s_str)


# Trotter circuit construction
import pennylane as qml
# Just basic Trotterization
# Just basic Trotterization
def evolve_circuit(pstr, on_wire:int, 
                   coeff:float, t:float, 
                   imaginary=False):
    """Return P evolution of exp(-i *t * coeff * P) or exp(- t * coeff*P)
    if `imaginary` is `True`.

    Args:
        pstr (_type_): Pauli string
        on_wire (int): Position of rotation gate
        coeff (float): Coefficient of Pauli term
        t (float): time
        imaginary (bool, optional): Evolution type REAL or IMAGINARY. Defaults to False.
    """
    act_wires=[]
    #basis_transform
    for i, s in enumerate(pstr):
        if s == "I":
            continue
        else: 
            act_wires.append(i)
            if s=="X":
                qml.Hadamard(wires=i)
            elif s=="Y":
                qml.adjoint(qml.S(wires=i))
                qml.Hadamard(wires=i)
    
    # CNOT
    if on_wire not in act_wires:
        on_wire = act_wires[0]
    for ai in act_wires:
        if on_wire == ai:
            continue
        qml.CNOT(wires=[ai, on_wire])
    if imaginary:
        dtau = t
        gamma = np.abs(coeff) # pure complex number 
        phi = 2*np.arccos(np.exp(-2*gamma*dtau))
        qml.ctrl(qml.RX, on_wire)(phi, wires=len(pstr))

    else:
        phi = coeff * t
        qml.RZ(phi, wires=on_wire, id=pstr)
    
    # CNOT
    for ai in reversed(act_wires):
        if on_wire == ai:
            continue
        qml.CNOT(wires=[ai, on_wire])
    # Reverse
    for i, s in enumerate(pstr):
        if s == "I" or s=="Z":
            continue
        elif s=="X":
            qml.Hadamard(wires=i)
        elif s=="Y":
            qml.Hadamard(wires=i)
            qml.S(wires=i)

            
    