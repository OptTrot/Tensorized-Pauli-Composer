## XZ code


## Tensorized Decomposition


With xz-code, `nx, nz`, the index of the corresponding Pauli term is determined with single bit opeartor, where, `nr=nz` and `nc=nx^nz`. 
Reverse transformation is defined as `nz=nr` and `nx = nr^nc`.

```
b = a^(a^b)
```


## C-Python interfaces

* ctypes
* Cython
* cffi

There are some utils in numpy library.

https://numpy.org/doc/stable/reference/routines.ctypeslib.html#module-numpy.ctypeslib