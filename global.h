#ifndef GLOBAL_SOR
#define GLOBAL_SOR
#include <stdio.h>

typedef struct _matrix {
    double *vals;
    int *col_ind;
    int *row_start_ind;
    int nelements;
    int n;
} matrix;

typedef struct _values {
    double *v;
    int n;
} values;

void test();
double mul_matrix_row_ind(matrix *m, double *vector, int row_ind);
int double_equals(double a, double b);
void forward_subst(matrix *m, double *x, double *b);
matrix get_iterative_upper_matrix(matrix *a, double w);
matrix get_iterative_lower_matrix(matrix *a, double w);
double* get_iterative_vector(matrix *a, double w, double *b);
void mul_matrix_row(matrix *m, double *vector, double *out, int start_row, int end_row);
void add_vector(double *x, double *b, int n);
matrix create_matrix_from_file(FILE *in);
void solve(matrix *a, double *x, double *b);
values create_vector_from_file(FILE *in);
#endif
