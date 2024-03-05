/*----------------------*/
// Library
#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "binary.h"
/*----------------------*/
// Type definition
typedef struct {
    unsigned long long nx;
    unsigned long long nz;
    size_t length;
} PauliCode;

typedef struct {
    char * str;
    size_t length;
} PauliString;
/*----------------------*/
static const char PstringMap[4] = "IZXY";
/*----------------------*/

PauliString get_pstr(char [], size_t);
PauliCode pstr2code(PauliString);
PauliString code2pstr(PauliCode);
//void code2pstr_pointer(PauliCode, PauliString *);
void del_pstr(PauliString);

PauliCode add(PauliCode, PauliCode);
bool commute(PauliCode, PauliCode);

// Utils
PauliCode * PauliSet(size_t, size_t *);
PauliString * get_pstr_list(size_t, size_t);
void del_pstr_list(PauliString *, size_t);