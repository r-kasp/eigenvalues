#ifndef MATRIX_H
#define MATRIX_H

int read_matrix(double *a, int n, char *name);
void init_matrix(double *a, int n, int k);
void init_vector(double *b, double *a, int n);
void print_matrix(double *a, int m, int n, int r);
void copy_matrix(double *x, double *a, int n);
void copy_vector(double *c, double *b, int n);
double f(int s, int n, int i, int j);
double norma_vector(double *a, int n);
double norma_matrix(double *a, int m, int n);

int to_partly_triangular_reflection(double *a, int n, double eps, double norm);
int qr_rotate(double *a, double *cos, double *sin, int n, int size, double eps, double norm);
int find_eigenvalues(double *egn, double *a, int n, double eps, double norm, double & t1, double & t2, int & it);


enum RETURN_CODES
{
    SUCCESS,
    ERROR_READ,
    ERROR_OPEN,
    SING_MATRIX,
    ERROR,
    DIVISION_ZERO,
    NEGATIVE_SQUARE
};

#endif
