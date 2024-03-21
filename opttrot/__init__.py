
from __future__ import annotations

VERSION = "0.0.1"

# Standard modules
from typing import *
from collections import OrderedDict
from itertools import combinations, combinations_with_replacement as re_combi, product
from functools import reduce
from pathlib import Path

from sys import float_info
FLOAT_EPS = 1E4 * float_info.min
float_tol = 1E-8


# Dependenct modules
import numpy as np
from scipy import linalg
import pandas as pd


from opttrot.utils import (
    commute_reggio_df, integer_order_map, 
    ternsorized_decomposition,
    frobenius_inner, krons, get_coef,
    get_decomposition, pstr_to_matrix,
    pstr_to_xz_fam_code,
    xz_fam_code_to_pstr,
    index_xzcode,
    pauli_X, pauli_Y, pauli_Z, I)


class Hamiltonian:
    def __init__(self, 
                 H:np.matrix, 
                 pauli_basis:Union[None, dict]=None,
                 tols=(1E4*float_tol , float_tol), 
                 commute_map = True):
        assert len(H.shape) ==2, f"H must be 2dim matrix. current: {H.shape}."
        n1, n2 = H.shape
        assert n1 == n2, f"Hamiltonian must be square matrix. Current:{(n1, n2)}."
        assert np.allclose(H, H.getH(), *tols), f"Hamiltonian must be a hermite matrix. Relative, absolute tolerance, {tols}."
        assert bin(n1)[2:].count("1") == 1, f"Dimension must be a 2^n. Current:{n1}."
        
        self.Hamiltonian = H
        
        if pauli_basis is None:
            pauli_basis = self.H_to_p_poly(H, tols[1])
        # None or Dataframe
        self.local_decomposition = get_decomposition(pauli_basis)  

        # Commute term
        self.commute_map = self.get_commuting_map() if commute_map else None
        self.commute_map_exist = True if commute_map else False
        self.qubit_num = len(bin(H.shape[0])[3:]) # Consider a 1 bit position of 2^n integer.
    #--------------------------------------------------------------
    def get_commuting_map(self):
        df = self.local_decomposition
        edge_df = pd.DataFrame(combinations(df["Pstring"].values ,2), columns=['source', 'target'])
        
        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="source", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zs", "X": "Xs"}, inplace=True)
        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="target", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zt", "X": "Xt"}, inplace=True)
        
        edge_df["commute"] = edge_df[["Zs", "Xs", "Zt", "Xt"]].apply(lambda x: int(commute_reggio_df(x)), axis=1)
        return edge_df
    def applying_weight_func(self, weight_func:Callable, columns:Iterable, name="Weight", inplace=False):
        """Calculate value based on the exist columns, `wieght_func` is a function to calculate the new value based on the exist values.
        `columns` is a column names or order on internal pandas dataframe. 
        The result would be saved in `name` column of `.commute_map` Pandas dataframe, if `inplace` is `True` else the result is returned by function. 
        If there is a column `name` then the column is replaced by the result, or new column is created. 
        Default value is "Weight".

        Example code:

        .. highlight:: python
        .. code-block:: python
            H_example = Hamiltonian(...)

            col_names = ["column1", "column2"]
            col_name_weight = "result"
            def weight_func(cols):
                c1 = cols.iloc[0]
                c2 = cols.iloc[1]
                ...
                return col_value
                
            H_example.applying_weight_func(weight_func, col_names, name=col_name_weight, inplace=False)

        Args:
            weight_func (Callable): _description_
            columns (_type_): _description_
            name (str, optional): _description_. Defaults to "Weight".
            inplace (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if not self.commute_map_exist:
            self.commute_map = self.get_commuting_map()
        if isinstance(columns[0], str): 
            name_series = self.commute_map.loc[:, columns].apply(weight_func, axis=1)
        elif isinstance(columns[0], int):
            name_series = self.commute_map.iloc[:, columns].apply(weight_func, axis=1)
        
        if inplace:
            self.commute_map[name] = name_series
        else:
            return name_series
    def save_as(self, filepath:Union[Path, str]): # In progress
        # Design a data model,
        # HDF? Pandas dataframe?
        raise NotImplementedError
        if isinstance(filepath, str):
            filepath = Path(filepath)
        pass
    
    #--------------------------------------------------------------
    @property
    def pauli_decomposition(self):
        return self.local_decomposition.loc[["Pstring", "Coef"]]
    @property
    def xz_family(self):
        return self.local_decomposition.loc[["Z", "X", "Coef"]]
    @property
    def latin_matrix(self):
        return self.local_decomposition.pivot(index="X", columns="Z", values="Coef") 
    @property
    def graph_edge(self):
        assert self.commute_map_exist, "No commute map exist. Execute <Hamiltonain>.get_commuting_map()."
        return self.commute_map[self.commute_map["source"] != self.qubit_num*("I")]
        
    #--------------------------------------------------------------
    @classmethod
    def from_latin_matrix(cls:Hamiltonian, 
                      l_matrix:np.matrix, 
                      xz_famileis:Tuple[Iterable, Iterable])->Hamiltonian: 
        # In progress
        pass
    @classmethod
    def from_pauli_polynomial(cls:Hamiltonian, 
                               p_poly:Union[dict, np.ndarray], 
                               p_coef:Union[None, np.ndarray]=None,
                               *args)-> Hamiltonian:
        if not isinstance(p_poly, dict):
            p_dict = {}
            for p, coef in zip(p_poly, p_coef):
                p_dict[p] = coef
        else: p_dict = p_poly

        keys = list(p_dict.keys())
        m_size = int(2**(p_dict[keys[0]]))
        H = np.zeros((m_size, m_size))
        
        for key in keys:
            H += p_dict[key] * pstr_to_matrix(key)
        return cls(np.asmatrix(H), p_poly, *args)
    @classmethod
    def from_data(cls:Hamiltonian, file_path)->Hamiltonian: 
        # In progress
        pass
    #------------------------------
    # Basic utils for hamiltonian analysis
    @staticmethod
    def p_poly_to_H(p_poly:dict):
        """Convert pauli-polynomial of dictionary form to total Hamiltonian matrix.
        The given polynomial must be a dictionary whose keys are pauli-terms and the values are coefficient.

        Args:
            pstrs (dict): _description_
        """
        n = len(list(p_poly.keys())[0])
        dim = int(2**n)
        shape = (dim, dim)
        result = np.asmatrix(np.zeros(shape, dtype=complex))
        for pstr in p_poly:
            coef = p_poly[pstr]
            result += coef*pstr_to_matrix(pstr)
        return result
    @staticmethod
    def H_to_p_poly(H, tol=float_tol, include_zeros=False)->dict:
        """Convert Hermit matrix to pauli-dict polynomial.

        Args:
            H (_type_): Hamiltonian Hermit matrix.
            tol (_type_, optional): Precision tolerance value. Defaults to float_tol.
            include_zeros (bool, optional): Including zero coefficient terms or not. Defaults to False.

        Returns:
            dict: Pauli term and coefficient dictionary.
        """
        n1, n2 = H.shape
        assert n1 == n2, "The given matrix must be a square matrix."
        poly = {}
        n = int(np.log2(n1))
        H_decom = ternsorized_decomposition(H)
        nr, nc = H_decom.shape

        for i in range(nr):
            for j in range(nc):
                coef = H_decom[i, j]
                coef = 0 if np.absolute(coef) < tol else coef
                if include_zeros or coef != 0:
                    poly[xz_fam_code_to_pstr(index_xzcode(i,j), n)] = coef
        #n = len(bin(H.shape[0])[3:])
        #p_mat, p_str = Hamiltonian.generate_pauli_terms(n)
        #for p_m, p_str in zip(p_mat, p_str):
        #    coef = frobenius_inner(p_m, H)
        #    coef = 0 if np.absolute(coef) < tol else coef
        #    if include_zeros:
        #        poly[p_str] = coef
        #    elif coef != 0:
        #        poly[p_str] = coef
        return poly
    @staticmethod
    def p_poly_to_latin(p_poly:dict, full=False)->Tuple[np.ndarray, list, list]:
        p_terms = list(p_poly.keys())
        x_fam = []
        z_fam = []
        for p in p_terms:
            nx, nz = pstr_to_xz_fam_code(p)
            x_fam.append(nx)
            z_fam.append(nz)
            
        x_fam_unique = np.unique(x_fam)
        z_fam_unique = np.unique(z_fam)
        x_l = x_fam_unique.size
        z_l = z_fam_unique.size
        
        x_l_map = integer_order_map(x_fam)
        z_l_map = integer_order_map(z_fam)
        
        latin_matrix = np.zeros(shape=(x_l, z_l), dtype=complex)
        for p, x_i, z_j in zip(p_poly.values(), x_fam, z_fam):
            xi, zi = x_l_map[x_i], z_l_map[z_j]

            latin_matrix[xi, zi] = p 
        return latin_matrix, x_fam_unique, z_fam_unique
    @staticmethod
    def generate_pauli_terms(
        qubit_num:int, 
        only:Literal["both", "string", "matrix"]="both")-> Union[Tuple[Iterable, Iterable], Iterable]:
        """Generate full set of pauli-terms in matrix and strings of `n` number of qubit system.

        Args:
            qubit_num (int): _description_
            only (Literal[&quot;both&quot;, &quot;string&quot;, &quot;matrix&quot;], optional): _description_. Defaults to "both".

        Returns:
            _type_: _description_
        """
        n = int(qubit_num)
        assert n >0, "The given argument must be a positive natural number."
        
        p_xs =  Hamiltonian.get_pauli_family_matrix(n, fam="X")
        p_zs =  Hamiltonian.get_pauli_family_matrix(n, fam="Z")
        p_xs_str = Hamiltonian.get_pauli_family_string(n, fam="X")
        p_zs_str = Hamiltonian.get_pauli_family_string(n, fam="Z")

        result = []
        if only=="both" or only=="matrix":
            p_g = []
            p_g_str =[]
            for x_i, x_str in zip(p_xs, p_xs_str):
                for z_j, z_str in zip(p_zs, p_zs_str):
                    g = x_i@z_j
                    g_coef, g_str = get_coef(x_str, z_str)
                    p_g.append(g_coef*g)
                    p_g_str.append(g_str)
            result.append(p_g) 
            if only =="both":
                result.append(p_g_str)
        elif only=="string":
            p_g_str = []
            for x_str in p_xs_str:
                for z_str in p_zs_str:
                    p_g_str.append(g_str)
            result.append(p_g_str)
        return result
    @staticmethod
    def get_pauli_family_string(n, fam="Z"):
        return list(map(lambda x: "".join(x), product(f"I{fam}", repeat=int(n))))
    @staticmethod
    def get_pauli_family_matrix(n:int, fam="Z")->Iterable[np.matrix]:
        """Get pauli_family of `n` qubits of `fam` family. 

        Args:
            n (int): Number of qubits. The output matrices are :math:`2^n`.
            fam (str, optional): Type of Pauli-family of X, Y, or Z. Defaults to "Z".

        Returns:
            Iterable[np.matrix]: list of Pauli-matrices
        """
        G = pauli_Z if fam=="Z" else (pauli_X if fam=="X" else pauli_Y)

        return list(map(krons, product([I, G], repeat=int(n))))
        