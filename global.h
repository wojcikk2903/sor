#ifndef GLOBAL_SOR
#define GLOBAL_SOR

typedef struct _matrix {
    double *vals;
    int *col_ind;
    int *row_start_ind;
    int nelements;
    int n;
} matrix;

void test();
double mul_matrix_row_ind(matrix *m, double *vector, int row_ind);
int double_equals(double a, double b);
void forward_subst(matrix *m, double *x, double *b);
matrix get_iterative_upper_matrix(matrix *a, double w);
matrix get_iterative_lower_matrix(matrix *a, double w);
#endif
