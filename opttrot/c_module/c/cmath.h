#include <stdio.h>
#include <math.h>

typedef struct{
    double real;
    double img;
} complex;

typedef struct
{   size_t nr;
    size_t nc;
    double * r_data;
    double * i_data;
} complex_matrix;

typedef struct
{   size_t nr;
    size_t nc;
    double * data;
} matrix;

complex get_complex(double, double);
complex c_conj(complex);
double c_radi(complex);
double c_args(complex);
complex c_inverse(complex);
double c_real(complex);
double c_img(complex);

complex c_add(complex, complex);
complex c_sub(complex, complex);
complex c_mul(complex, complex);
complex c_div(complex, complex);

complex c_add_real(complex, double);
complex c_add_int(complex, int);

complex c_sub_real(complex, double);
complex c_sub_int(complex, int);

//--------------------------------------

void matrix_free(matrix);
matrix zeros(size_t nr, size_t nc);
matrix eyes(size_t n);

void mat_add(matrix, matrix, matrix);
void mat_sub(matrix, matrix, matrix);
void mat_scale(matrix, matrix, matrix);
void mat_mul(matrix, matrix, matrix);
void mat_div(matrix, matrix, matrix);

void mat_prod(matrix, matrix, matrix);

double mat_trace(matrix);
double mat_transpose(matrix, matrix);
double mat_norm(matrix);

void c_matrix_free(complex_matrix);
complex_matrix c_zeros(size_t nr, size_t nc);
complex_matrix c_eyes(size_t n);
complex_matrix c_from_matrix(matrix);
complex_matrix c_mat_complex_scale(matrix);

void c_mat_add(complex_matrix, complex_matrix, complex_matrix);
void c_mat_sub(complex_matrix, complex_matrix, complex_matrix);
void c_mat_scale(complex_matrix, complex_matrix, complex_matrix);
void c_mat_mul(complex_matrix, complex_matrix, complex_matrix);
void c_mat_div(complex_matrix, complex_matrix, complex_matrix);

void c_mat_prod(complex_matrix, complex_matrix, complex_matrix);

double c_mat_trace(complex_matrix);
double c_mat_transpose(complex_matrix, complex_matrix);
void c_mat_dagger(complex_matrix, complex_matrix);
