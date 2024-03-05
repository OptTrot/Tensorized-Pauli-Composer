#include "cmath.h"

complex get_complex(double x, double y){
    complex c = {.real=x, .img=y};
    return c;
}
complex c_conj(complex c){
    c.img = -c.img;
    return c;
}
double c_radi(complex c){
    return sqrt(pow(c.real,2) + pow(c.img,2));
}
double c_args(complex c){
    double r = c_radi(c);
    return asin(c.real/r);
}
complex c_inverse(complex c){
    double r2 = pow(c_radi(c),2);
    complex conj = c_conj(c);
    conj.real /=r2;
    conj.img /=r2;
    return conj;
}
double c_real(complex c){
    return c.real;
}
double c_img(complex c){
    return c.img;
}

complex c_add(const complex c1, const complex c2){
    complex c3 = {
        .real = c1.real + c2.real, 
        .img  = c1.img + c2.img
        };
    return c3;
}
complex c_sub(const complex c1, const complex c2){
    complex c3 = {
        .real = c1.real - c2.real, 
        .img  = c1.img - c2.img
        };
    return c3;
}
complex c_mul(const complex c1, const complex c2){
    complex c3 = {
        .real = c1.real*c2.real - c1.img*c2.img, 
        .img  = c1.real*c2.img + c1.img * c2.real
        };
    return c3;
}
complex c_div(complex c1, complex c2){
    double r2 = pow(c_radi(c2), 2);
    complex c2_conj = c_conj(c2);
    complex c3 = c_mul(c1, c2_conj);
    c3.real/=r2;
    c3.img/=r2;
    return c3;
}

complex c_add_real(complex c, double r){
    c.real += r;
    return c;
}
complex c_add_int(complex c, int i ){
    c.real += i;
    return c;
}

complex c_sub_real(complex c, double r){
    c.real -= r;
    return c;
}
complex c_sub_int(complex c, int i ){
    c.real -= i;
    return c;
}
