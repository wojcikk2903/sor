#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "global.h"

double mul_matrix_row_ind(matrix *m, double *vector, int row_ind)
{
    int cur_row_start = m->row_start_ind[row_ind];
    int cur_row_end = row_ind < m->n-1 ? m->row_start_ind[row_ind+1] : m->nelements;
    double result = 0.0;

    for (int i = cur_row_start; i < cur_row_end; i++)
        result += m->vals[i]*vector[m->col_ind[i]];

    return result;
}

void mul_matrix_row(matrix *m, double *vector, double *out, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
        out[i] = mul_matrix_row_ind(m, vector, i);
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

matrix get_iterative_lower_matrix(matrix *a, double w)
{
    matrix it_matrix;
    matrix am = *a;
    it_matrix.vals = malloc(am.nelements*sizeof(double));
    it_matrix.col_ind = malloc(am.nelements*sizeof(int));
    it_matrix.row_start_ind = malloc(am.n*sizeof(int));
    it_matrix.n = am.n;

    int it_el_ind = 0;
    for (int i = 0; i < am.n; i++)
    {
        int diag_el_ind = find_diag_el_ind(a, i);
        double diag_el_value = am.vals[diag_el_ind];
        
        it_matrix.row_start_ind[i] = it_el_ind;
        for (int j = am.row_start_ind[i]; j < diag_el_ind; j++)
        {
            it_matrix.vals[it_el_ind] = w*am.vals[j]/diag_el_value;
            it_matrix.col_ind[it_el_ind] = am.col_ind[j];
            it_el_ind++;
        }
        // it_matrix.vals[it_el_ind++] = 1+w;
        it_matrix.vals[it_el_ind] = 1;
        it_matrix.col_ind[it_el_ind] = i;
        it_el_ind++;
    }
    it_matrix.nelements = it_el_ind;
    // it_matrix.row_start_ind[am.n-1] = it_el_ind-1;

    return it_matrix;
}

matrix get_iterative_upper_matrix(matrix *a, double w)
{
    matrix it_matrix;
    matrix am = *a;
    it_matrix.vals = malloc(am.nelements*sizeof(double));
    it_matrix.col_ind = malloc(am.nelements*sizeof(int));
    it_matrix.row_start_ind = malloc(am.n*sizeof(int));
    it_matrix.n = am.n;

    int it_el_ind = 0;
    for (int i = 0; i < am.n-1; i++)
    {
        int diag_el_ind = find_diag_el_ind(a, i);
        double diag_el_value = am.vals[diag_el_ind];
        
        it_matrix.row_start_ind[i] = it_el_ind;
        // it_matrix.vals[it_el_ind] = 1-2*w;
        it_matrix.vals[it_el_ind] = 1-w;
        it_matrix.col_ind[it_el_ind] = i;
        it_el_ind++;
        for (int j = diag_el_ind+1; j < am.row_start_ind[i+1]; j++)
        {
            it_matrix.vals[it_el_ind] = -w*am.vals[j]/diag_el_value;
            it_matrix.col_ind[it_el_ind] = am.col_ind[j];
            it_el_ind++;
        }
    }
    // it_matrix.vals[it_el_ind] = 1-2*w;
    it_matrix.vals[it_el_ind] = 1-w;
    it_matrix.col_ind[it_el_ind] = am.n-1;
    it_matrix.nelements = it_el_ind+1;
    it_matrix.row_start_ind[am.n-1] = it_el_ind;

    return it_matrix;
}

double* get_iterative_vector(matrix *a, double w, double *b)
{
    double *it_vec = malloc(a->n*sizeof(double));
    for (int i = 0; i < a->n; i++)
    {
        int diag_el_ind = find_diag_el_ind(a, i);
        double diag_el_val = a->vals[diag_el_ind];
        it_vec[i] = w*b[i]/diag_el_val;
    }
    return it_vec;
}

void add_vector(double *x, double *b, int n)
{
    for (int i = 0; i < n; i++)
        x[i] += b[i];
}

void iteration_step(matrix lower_matrix, matrix upper_matrix, double *x, int xn, double *b)
{

}

int main(int argc, char *argv[])
{
#ifndef DEBUG
    int nproc;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Finalize();
#endif

#ifdef DEBUG
    test();
#endif
    return 0;
}
