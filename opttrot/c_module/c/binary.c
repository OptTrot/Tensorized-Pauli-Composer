// Macro
#define XZ_COEF_POW(x, z) bit_count(x&z)
#define ENCODE_XZCODE(x, z) (insert_zeros(x) <<1) + insert_zeros(z)


/*----------------------*/
// Library
#include "binary.h"
/*----------------------*/
size_t bit_count(ull n){
    //Brian Kernighan’s Algorithm
    int num = 0;
    while (n){n = n&(n-1);num++;}
    return num;
}

ull xz_coef_pow(unsigned int x, unsigned int z){
    return XZ_COEF_POW(x, z);
}

ull insert_zeros(ull n){ 
    // Insert zeros to nested gaps 
    // 110 -> 010100 = (01)(01)(00)
    int result = 0;
    int bit_position = 0;
    int r_most_bit = 0;
    while (n>0){
        r_most_bit = n&1; //  (...101)&(...001) = (...001), 1 or 0
        result = result | (r_most_bit << (bit_position <<1));
        n = n>>1;
        bit_position++;
    }
    return result;
}
unsigned int encode_xzcode(unsigned int x, unsigned int z){
    return ENCODE_XZCODE(x, z);
}

void decode_pcode(
    unsigned int p_int, 
    size_t p_len, 
    unsigned int * x, 
    unsigned int * z){
    *x = 0;
    *z = 0;

    unsigned int mask =1;
    int bit_position =0;
    while (p_int > 0) {
        if (bit_position % 2 == 0) {
            // Even bit position: Add to x
            *x |= (p_int & 1) << (bit_position >> 1);
        } else {
            // Odd bit position: Add to z
            *z |= (p_int & 1) << (bit_position >> 1);
        }
        p_int >>= 1; // Shift right to get the next bit
        bit_position++;
    }
}

void print_binary(ull n, size_t l){
    size_t i= 0;
    while (n || i<l) {
    if (n & 1) printf("1");
    else printf("0");
    n >>= 1;
    i++;
    }
}

/*
int main(int args, char **argvs){
    ull num[6] = {7, 12, 10, 42, 291, 1035};

    for(size_t i=0; i <6; i ++){
        printf("%llu\n", num[i]);
        print_binary(num[i], 10);
        printf("\n");
        printf("bit_count routine: %zu", bit_count(num[i]));
        printf("\n");
        printf("Insert zeros\n");
        print_binary(insert_zeros(num[i]), 20);
        printf("\n");
    }
    return 0;
}
*/