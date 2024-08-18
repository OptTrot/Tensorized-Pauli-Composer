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

//--------------------------------------

void matrix_free(matrix m){
    free(m.data)
}
matrix zeros(size_t nr, size_t nc){
    matrix m = {
        .nr = nr, 
        .nc = nc,
        .data = (double *)malloc(nr*nc * sizeof(double))
        };
    for(size_t i = 0; i <nr, i++){
        for(size_t j = 0; j <nc, j++){
            m.data[i][j] = 0;
        }
    }
    return m;
}
matrix eyes(size_t n){
    matrix m = {
        .nr = n, 
        .nc = n,
        .data = (double *)malloc(nr*nc * sizeof(double))
        };

    for(size_t i = 0; i <m.nr, i++){
        for(size_t j = 0; j <m.nc, j++){
            m.data[i][j] = (int)(i==j);
        }
    }
    return m;
}


void mat_scale(double s, matrix m1, matrix m2){
    for(size_t i = 0; i <m1.nr, i++){
        for(size_t j = 0; j <m1.nc, j++){
            m2.data[i][j] = s*m1.data[i][j];
        }
    }
}

void mat_add(matrix m1, matrix m2, matrix m3){
    //Error m1, m2, m3 dimension check;

    for(size_t i = 0; i <m.nr, i++){
        for(size_t j = 0; j <m.nc, j++){
            m3.data[i][j] = m1.data[i][j] + m2.data[i][j];
        }
    }
}
void mat_sub(matrix m1, matrix m2, matrix m3){
    //Error m1, m2, m3 dimension check;

    for(size_t i = 0; i <m.nr, i++){
        for(size_t j = 0; j <m.nc, j++){
            m3.data[i][j] = m1.data[i][j] - m2.data[i][j];
        }
    }
}

void mat_mul(matrix m1, matrix m2, matrix m3){
    //Error m1, m2, m3 dimension check;

    for(size_t i = 0; i <m.nr, i++){
        for(size_t j = 0; j <m.nc, j++){
            m3.data[i][j] = m1.data[i][j] * m2.data[i][j];
        }
    }
}
void mat_div(matrix m1, matrix m2, matrix m3){
    //Error m1, m2, m3 dimension check;

    for(size_t i = 0; i <m1.nr, i++){
        for(size_t j = 0; j <m1.nc, j++){
            m3.data[i][j] = m1.data[i][j] / m2.data[i][j];
        }
    }
}

matrix mat_prod(matrix m1, matrix m2, matrix m3){
    size_t i, j, k;

    for(i =0; i< m1.nr; i++){
        for(j =0; j< m2.nc; j++){
            double s = 0;
            for(k =0; k< m1.nc; k++){
                s+= (m1[i][k] * m2[k][j]);
            }
            m3[i][j] = s;
        }
    }
    return m3
}

double mat_trace(matrix m){
    double tr = 0.;
    for(size_t i = 0; i < m.nr; i++){
        tr += m.data[i][i];
    }
    return tr;
}
matrix mat_transpose(matrix m1, matrix m2){
    for(size_t i = 0; i <m1.nr, i++){
        for(size_t j = 0; j <m1.nc, j++){
            m2.data[j][i] = m1.data[i][j]
        }
    }
    return m2;
}
double mat_norm(matrix);
