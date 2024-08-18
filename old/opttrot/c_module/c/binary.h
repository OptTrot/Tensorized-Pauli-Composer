// Macro
#define XZ_COEF_POW(x, z) bit_count(x&z)
#define ENCODE_XZCODE(x, z) (insert_zeros(x) <<1) + insert_zeros(z)

// Library

#include <stdio.h>

// typedef 
typedef unsigned long long ull;

size_t bit_count(ull);
ull  insert_zeros(ull);
ull  xz_coef_pow(unsigned int, unsigned int);
void decode_pcode(
    unsigned int, 
    size_t, 
    unsigned int *, 
    unsigned int * );

void print_binary(ull, size_t);
