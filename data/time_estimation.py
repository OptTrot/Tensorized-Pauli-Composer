#Default routine
import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import Pauli, PauliList
import pennylane as qml
from pennylane.pauli import string_to_pauli_word, pauli_word_to_matrix
from paulicompsoser import pauli_composer as pc
import os
import sys
module_path = os.path.abspath(os.path.join('..\\..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from opttrot import pauli_mani as pm
import timeit


def naive(pauli_poly, ndim):
    mat= np.zeros((ndim, ndim), dtype=complex)
    # Naive tensor
    for pauli in pauli_poly.poly:
        mat += pauli.to_matrix()
    return mat
def pennylane_naive(p_dict, ndim):
    mat = np.zeros((ndim, ndim), dtype=complex)
    for p, coef in p_dict:
        mat += coef* qml.pauli.pauli_word_to_matrix(p)
    return mat
def qiskit_composer(coefs, qiskit_paulis, ndim):
    mat = np.zeros((ndim, ndim), dtype=complex)
    for p, coef in zip(qiskit_paulis, coefs):
        mat += coef* p.to_matrix()
    return mat
def qiskit_composer_list(qiskit_list, coefs, ndim):
    paulimats = qiskit_list.to_matrix(array=True)
    mat = np.zeros((ndim, ndim), dtype=complex)
    for coef, p in zip(coefs, paulimats):
        mat += coef*p 
    return mat
def pauli_composer(pauli_poly, ndim):
    p_composers = [pc.PauliComposer(p.string, p.coef) for p in pauli_poly.poly]
    mat = np.zeros((ndim, ndim), dtype=complex)
    for p in p_composers:
        mat += p.to_matrix()
    return mat