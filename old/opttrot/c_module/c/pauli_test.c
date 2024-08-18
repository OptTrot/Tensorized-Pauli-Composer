#include "pauli_alg.h"


int main(int args, char **argvs){
    //Test the routines
    char ch[3] = "ZXY";
    PauliString pstr = get_pstr(ch, 3);
    PauliCode pcode = pstr2code(pstr);
    PauliString pstr2 = code2pstr(pcode);

    printf("Pcode of %s, length:%zu\n", pstr.str, pstr.length);
    printf("%llu, %llu, length: %zu\n", pcode.nx, pcode.nz, pcode.length);
    printf("Pstr: %s, length:%zu\n", pstr2.str, pstr2.length);

    del_pstr(pstr);
    del_pstr(pstr2);

    size_t n = 4;
    printf("%zu-quit Pauli set generation.", n);
    size_t p_dim = 0;
    PauliCode * plist = PauliSet(n, &p_dim);
    PauliString * pslist = (PauliString *)malloc(p_dim * sizeof(PauliString));
    
    for(size_t i=0; i<(int)(pow(4,n));i++){
        PauliString pst = code2pstr(plist[i]);
        printf("%zu, %s, (%llu, %llu) \n", i, pst.str, plist[i].nx, plist[i].nz);
        del_pstr(pst);
    }
    free(plist);

    printf("Pauli string synthesis:\n");

    // Pauli term algebra.
    PauliCode pcode_list[2] = {
        {5, 3, 5}, 
        {2, 7, 5}
        };

    PauliCode pcode3 = add(pcode_list[0], pcode_list[1]);
    
    PauliString pstr001 = code2pstr(pcode_list[0]);
    PauliString pstr002 = code2pstr(pcode_list[1]);
    PauliString pstr003 = code2pstr(pcode3);

    printf("\n%s + %s = %s\n", 
    pstr001.str,  
    pstr002.str,  
    pstr003.str);

    printf("%llu, %llu\n", pcode_list[0].nx, pcode_list[0].nz);
    printf("%llu, %llu\n", pcode_list[1].nx, pcode_list[1].nz);
    printf("%llu, %llu\n", pcode3.nx, pcode3.nz);
    return 0;
}