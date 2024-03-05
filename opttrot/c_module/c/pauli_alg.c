#include "pauli_alg.h"

/*----------------------*/
PauliString get_pstr(char st[], size_t l){
    PauliString pst;
    pst.str =  (char*)malloc(l * sizeof(char));
    pst.length = l;

    for(size_t i =0 ; i <l; i++){
        switch(st[i]){
            case 'I':
            case 'Z':
            case 'X':
            case 'Y':
                pst.str[i] = st[i];
                continue;
            default:
                // Raise error
                break;
        }
    }
    return pst;
}

void del_pstr(PauliString pstr){free(pstr.str);}

PauliCode pstr2code(PauliString pst){
    // Convert the given pauli-string to xz-code
    PauliCode pcode;
    pcode.nx = 0;
    pcode.nz = 0;
    pcode.length = pst.length;

    unsigned int pint_x = 0;
    unsigned int pint_z = 0;
    for(size_t i=0; i<pst.length; i++){
        //I -> 00, 0
        //Z -> 01, 1
        //X -> 10, 2
        //Y -> 11, 3
        switch(pst.str[i]){
            case 'I':
                pint_x = 0;
                pint_z = 0;
                break;
            case 'Z':
                pint_x = 0;
                pint_z = 1;
                break;
            case 'X':
                pint_x = 1;
                pint_z = 0;
                break;
            case 'Y':
                pint_x = 1;
                pint_z = 1;
                break;
        }
        pcode.nx += pint_x; 
        pcode.nz += pint_z;
        pcode.nx <<= 1;
        pcode.nz <<= 1;
    }
    pcode.nx >>=1;
    pcode.nz >>=1;
    return pcode;
}

PauliString code2pstr(PauliCode pcode){
    // Convert the given xz-code to pauli string
    PauliString pstr;
    pstr.length = pcode.length;
    pstr.str = (char*)malloc(pcode.length * sizeof(char));
    ull a = insert_zeros(pcode.nx);
    a <<=1;
    ull b = insert_zeros(pcode.nz);
    ull c = a|b;
    ull end_int =0;

    for(size_t i = 0; i<pcode.length; i++){
        // Get end two bit
        end_int=c&3; 
        c >>= 2;
        pstr.str[pcode.length-1-i] = PstringMap[end_int];
    }
    return pstr;
}

/*
void code2pstr_pointer(PauliCode pcode, PauliString * pstr){
    (*pstr).length = pcode.length;
    (*pstr).str = (char*)malloc(pcode.length * sizeof(char));
    ull a = insert_zeros(pcode.nx);
    a <<=1;
    ull b = insert_zeros(pcode.nz);
    ull c = a|b;
    ull end_int =0;

    for(size_t i = 0; i<pcode.length; i++){
        // Get end two bit
        end_int=c&3; 
        c >>= 2;
        (*pstr).str[pcode.length-1-i] = PstringMap[end_int];
    }
}*/


PauliCode add(PauliCode p1, PauliCode p2){
    // raise error 
    if (p1.length != p2.length){
        // error implementation
        // no permitted to synethsis two differnet length string
    }
    PauliCode p3 = {
        p1.nx ^ p2.nx, 
        p1.nz ^ p2.nz, 
        p1.length};//xor bitwise
    return p3;
}
bool commute(PauliCode p1, PauliCode p2){
    // Reggio et al
    ull n1 = p1.nx & p2.nz;
    ull n2 = p2.nx & p1.nz;

    bool a = bit_count(n1)&1;// modular 2
    bool b = bit_count(n2)&1; 
    return a == b;
}

/*-----------------*/
// Utils

// Pauli set generation of n qubits

PauliCode * PauliSet(size_t n, size_t * p_dim){
    size_t max_n = (int)(pow(2, n)); // 0b111...1 of length n
    * p_dim = max_n * max_n;
    size_t index = 0;
    PauliCode * pcode_list = (PauliCode *)malloc(max_n*max_n*sizeof(PauliCode));
    
    
    printf("\n Max_n:%zu,\n", max_n);

    for(size_t i=0; i<max_n;i++){
        for(size_t j=0; j<max_n;j++){
            index = max_n*i +j;
            pcode_list[index].nx = i;
            pcode_list[index].nz = j;
            pcode_list[index].length = n;
        }
    }
    return pcode_list;
}

PauliString * get_pstrs_list(size_t n, size_t qubits){
    PauliString * pstrs = (PauliString *)malloc(n*sizeof(PauliString));
    char * p_strs = (char *)malloc(n* qubits* sizeof(char));

    for(size_t i=0; i < n; i++){
        pstrs[i].str = &p_strs[i*qubits];
        pstrs[i].length = qubits;
    }
    return pstrs;
}
void del_pstr_list(PauliString * pstrs, size_t n){
    for(size_t i=0; i < n; i++){
        free(pstrs[i].str);
    }
    free(pstrs);
}



