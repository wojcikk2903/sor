#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"


const double DELTA = 1e-3;
const int TRUE = 1;
const int FALSE = 0;

typedef struct _matrix {
    double *vals;
    int *col_ind;
    int *row_start_ind;
    int nelements;
    int n;
} matrix;

double mul_matrix_row_ind(matrix *m, double *vector, int row_ind)
{
    int cur_row_start = m->row_start_ind[row_ind];
    int cur_row_end = row_ind < m->n-1 ? m->row_start_ind[row_ind+1] : m->nelements;
    double result = 0.0;

    for (int i = cur_row_start; i < cur_row_end; i++)
        result += m->vals[i]*vector[m->col_ind[i]];

    return result;
}

void forward_subst(matrix *m, double *x, double *b)
{
    int el_ind = 0;
    int col = 0;
    x[0] = b[0]/m->vals[el_ind];

    for (int i = 1; i < m->n; i++)
    {
        el_ind = m->row_start_ind[i];
        col = m->col_ind[el_ind];
        double accumulated = 0.0;
        while (col < i)
        {
            accumulated += m->vals[el_ind]*x[col];
            el_ind++;
            col = m->col_ind[el_ind];
        }
        x[i] = (b[i]-accumulated)/m->vals[el_ind];
    }
}

int find_diag_el_ind(matrix *a, int row)
{
    int el_ind = a->row_start_ind[row];
    while (a->col_ind[el_ind] < row)
        el_ind++;

    return el_ind;
}

matrix get_iterative_matrix1(matrix *a, double w)
{
    matrix it_matrix;
    matrix am = *a;
    it_matrix.vals = malloc(am.nelements*sizeof(double));
    it_matrix.col_ind = malloc(am.nelements*sizeof(int));
    it_matrix.row_start_ind = malloc(am.n*sizeof(int));

    int it_el_ind = 0;
    for (int i = 0; i < am.n; i++)
    {
        int diag_el_ind = find_diag_el_ind(a, i);
        double diag_el_value = am.vals[diag_el_ind];
        
        it_matrix.row_start_ind[i] = it_el_ind;
        for (int j = am.row_start_ind[i]; j < diag_el_ind; j++)
        {
            it_matrix.vals[it_el_ind] = -w*am.vals[j]/diag_el_value;
            it_matrix.col_ind[it_el_ind] = am.col_ind[j];
            it_el_ind++;
        }
        it_matrix.vals[it_el_ind++] = 1-2*w;
    }

    return it_matrix;
}

int double_equals (double a, double b)
{
    double mod_diff = a < b ? b-a : a-b;
    return mod_diff < DELTA;
}

void test()
{
    printf("\n\n-------------------------------\n\n\tMnożenie przez wektor\n\n-------------------------------\n\n");
    matrix m;
    double vals[] = {2.0, 4.0, 8.0, 5.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double vector[] = {1.0, 1.0};

    double res1 = mul_matrix_row_ind(&m, vector, 0);
    double res2 = mul_matrix_row_ind(&m, vector, 1);

    int all_passed = TRUE;
    if (!double_equals(6.0, res1))
        all_passed = FALSE;
    if (!double_equals(13.0, res2))
        all_passed = FALSE;
    printf("res1 = %.2f\n", res1);
    printf("res2 = %.2f\n", res2);

    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

    printf("\n------------------------------\n\tPodstawianie w przód\n---------------------------------\n");
    double vals2[] = { 2.0, 1.0, 3.0};
    int col_ind2[] = {0, 0, 1};
    int row_start_ind2[] = {0, 1};
    m.vals = vals2;
    m.n = 2;
    m.nelements = 3;
    m.col_ind = col_ind2;
    m.row_start_ind = row_start_ind2;

    double vector2[] = {1.0, 1.0};
    double b2[] = {2.0, 8.0};
    forward_subst(&m, vector2, b2);
    printf("x[0] = %.2f\n", vector2[0]);
    printf("x[1] = %.2f\n", vector2[1]);
    all_passed = TRUE;
    if (!double_equals(1.0, vector2[0]))
        all_passed = FALSE;
    if (!double_equals(2.33, vector2[1]))
        all_passed = FALSE;

    double vals3[] = { 3.0, 2.0, 1.0, 1.0, 5.0, 8.0};
    int col_ind3[] = {0, 0, 1, 0, 1, 2};
    int row_start_ind3[] = {0, 1, 3};
    m.vals = vals3;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind3;
    m.row_start_ind = row_start_ind3;
    double vector3[] = {1.0, 1.0, 1.0};
    double b3[] = {1.0, 2.0, 3.0};
    forward_subst(&m, vector3, b3);
    printf("x2[0] = %.2f\n", vector3[0]);
    printf("x2[1] = %.2f\n", vector3[1]);
    printf("x2[2] = %.2f\n", vector3[2]);
    all_passed = TRUE;
    if (!double_equals(0.333, vector3[0]))
        all_passed = FALSE;
    if (!double_equals(1.333, vector3[1]))
        all_passed = FALSE;
    if (!double_equals(-0.5, vector3[2]))
        all_passed = FALSE;
    
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

    printf("\n--------------------------------\n\tSzukanie pierwszej macierzy do iteracji\n-----------------------\n");
    double vals4[] = { 3.0, 1.0, 2.0, 1.0, 5.0, 8.0};
    int col_ind4[] = {0, 1, 2, 0, 1, 2};
    int row_start_ind4[] = {0, 1, 3};
    m.vals = vals4;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind4;
    m.row_start_ind = row_start_ind4;
    matrix it_matrix1 = get_iterative_matrix1(&m, 1.5);
    printf("m1.vals[0] = %.2f\n", it_matrix1.vals[0]);
    printf("m1.vals[1] = %.2f\n", it_matrix1.vals[1]);
    printf("m1.vals[2] = %.2f\n", it_matrix1.vals[2]);
    printf("m1.vals[3] = %.2f\n", it_matrix1.vals[3]);
    printf("m1.vals[4] = %.2f\n", it_matrix1.vals[4]);
    all_passed = TRUE;
    if (!double_equals(-2.0, it_matrix1.vals[0]))
        all_passed = FALSE;
    if (!double_equals(-2.0, it_matrix1.vals[1]))
        all_passed = FALSE;
    if (!double_equals(-0.1875, it_matrix1.vals[2]))
        all_passed = FALSE;
    if (!double_equals(-0.9375, it_matrix1.vals[3]))
        all_passed = FALSE;
    if (!double_equals(-2.0, it_matrix1.vals[4]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");
}

int main(int argc, char *argv[])
{
    int nproc;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    #ifdef DEBUG
    if (rank == 0) 
        test();
    #endif
    MPI_Finalize();
    return 0;
}
