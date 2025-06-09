import numpy as np

def pdc(izstring):
    N = len(izstring)
    m = np.ones(2**N, dtype=np.int8)
    #m[0] = 1
    #_2l = 1
    for l, p in enumerate(reversed(izstring)):
        if p == "I":
            m[2**l: 2**(l+1)] = m[0: 2**l]
        else:
            m[2**l: 2**(l+1)] = -m[0: 2**l]
    return m
