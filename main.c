#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "global.h"

int tag = 99;

typedef struct _range {
    int begin;
    int end;
} range;

void print_vector(double *x, int n)
{
    for (int i = 0; i < n; i++)
        printf("%f\n", x[i]);
}

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
        it_matrix.vals[it_el_ind] = 1;
        it_matrix.col_ind[it_el_ind] = i;
        it_el_ind++;
    }
    it_matrix.nelements = it_el_ind;

    return it_matrix;
}

matrix get_iterative_upper_matrix(matrix *a, double w)
{
    matrix *it_matrix_upper = malloc(sizeof(*it_matrix_upper));
    matrix it_matrix = *it_matrix_upper;
    matrix am = *a;
    it_matrix.col_ind = malloc(am.nelements*sizeof(int));
    it_matrix.row_start_ind = malloc(am.n*sizeof(int));
    it_matrix.vals = malloc(am.nelements*sizeof(double));
    it_matrix.n = am.n;

    int it_el_ind = 0;
    for (int i = 0; i < am.n-1; i++)
    {
        int diag_el_ind = find_diag_el_ind(a, i);
        double diag_el_value = am.vals[diag_el_ind];
        
        it_matrix.row_start_ind[i] = it_el_ind;
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

void solve(matrix *a, double *x, double *b)
{
    double w = 1.0;
    int steps = 1000;
    double *it_vec = get_iterative_vector(a, w, b);
    matrix upper = get_iterative_upper_matrix(a, w);
    matrix lower = get_iterative_lower_matrix(a, w);
    double *z = malloc(a->n*sizeof(*z));

    for (int i = 0; i < steps; i++)
    {
        mul_matrix_row(&upper, x, z, 0, a->n);
        add_vector(z, it_vec, a->n);
        forward_subst(&lower, x, z);
    }
}

range compute_range(int nrows, int nchild, int child_id)
{
    range r;
    int rows_per_child = nrows / nchild;
    int rest  = nrows % nchild;
    r.begin = child_id * rows_per_child + 
        (child_id < rest ? child_id : rest);
    r.end = r.begin + rows_per_child + (child_id < rest ? 1 : 0);
    return r;
}

void send_matrix(matrix send, int child_id)
{
    MPI_Send(&send.n, 1, MPI_INT, child_id, tag, MPI_COMM_WORLD);
    MPI_Send(&send.nelements, 1, MPI_INT, child_id, tag, MPI_COMM_WORLD);
    MPI_Send(send.vals, send.nelements, MPI_DOUBLE, child_id, tag, MPI_COMM_WORLD);
    MPI_Send(send.col_ind, send.nelements, MPI_INT, child_id, tag, MPI_COMM_WORLD);
    MPI_Send(send.row_start_ind, send.n, MPI_INT, child_id, tag, MPI_COMM_WORLD);
}

void send_values(double *v, int n, int child_id)
{
    MPI_Send(&n, 1, MPI_INT, child_id, tag, MPI_COMM_WORLD);
    MPI_Send(v, n, MPI_DOUBLE, child_id, tag, MPI_COMM_WORLD);
}

values receive_vector()
{
    MPI_Status status;
    values v;
    MPI_Recv(&v.n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
    v.v = malloc(v.n*sizeof(*v.v));
    MPI_Recv(v.v, v.n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status );
    return v;
}

void wait_for_results(double *z, int n, int nproc, int child_id)
{
    MPI_Status status;
    range r = compute_range(n, nproc, child_id);
    MPI_Recv(z+r.begin, r.end-r.begin, MPI_DOUBLE, child_id, tag, MPI_COMM_WORLD, &status);
}

void solve_parallel(matrix *a, double *x, double *b, int nproc)
{
    double w = 1.0;
    int steps = 2;
    double *it_vec = get_iterative_vector(a, w, b);
    matrix upper = get_iterative_upper_matrix(a, w);
    matrix lower = get_iterative_lower_matrix(a, w);
    double *z = malloc(a->n*sizeof(*z));


    for (int i = 1; i < nproc; i++) 
        send_matrix(upper, i);

    for (int s = 0; s < steps; s++)
    {
        for (int i = 1; i < nproc; i++) 
            send_values(x, a->n, i);
        range r = compute_range(a->n, nproc, 0);
        mul_matrix_row(&upper, x, z, r.begin, r.end);

        for (int i = 1; i < nproc; i++)
            wait_for_results(z, a->n, nproc, i);

        add_vector(z, it_vec, a->n);
        forward_subst(&lower, x, z);
    }
    print_vector(x, a->n);
}


matrix receive_matrix()
{
    matrix rec;
    MPI_Status status;
    MPI_Recv(&rec.n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
    MPI_Recv(&rec.nelements, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
    rec.vals = malloc(rec.nelements*sizeof(double));
    rec.col_ind = malloc(rec.nelements*sizeof(int));
    rec.row_start_ind = malloc(rec.n*sizeof(int));
    MPI_Recv(rec.vals, rec.nelements, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status );
    MPI_Recv(rec.col_ind, rec.nelements, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
    MPI_Recv(rec.row_start_ind, rec.n, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
    return rec;
}

double* initialize_x(int n, double initial_value)
{
    double *x = malloc(n*sizeof(double));
    for (int i = 0; i < n; i++)
        x[i] = initial_value;

    return x;
}



int main(int argc, char *argv[])
{
#ifndef DEBUG
    int nproc;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        FILE *fm = argc > 1 ? fopen(argv[1], "r") : fopen("matrixA.dat", "r"); 
        FILE *fv = argc > 2 ? fopen(argv[2], "r") : fopen("vectorB.dat", "r"); 
        matrix a = create_matrix_from_file(fm);
        values b = create_vector_from_file(fv);

        double *x = initialize_x(a.n, 0.0);

        solve_parallel(&a, x, b.v, nproc);
    }
    else
    {
        matrix m = receive_matrix();
        double *z = malloc(m.n*sizeof(*z));
        while (1){
            values v = receive_vector();
            range r = compute_range(v.n, nproc, rank);
            mul_matrix_row(&m, v.v, z, r.begin, r.end);
            MPI_Send(z+r.begin, r.end-r.begin, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
#endif

#ifdef DEBUG
    // test();
    FILE *fm = argc > 1 ? fopen(argv[1], "r") : fopen("matrixA.dat", "r"); 
    FILE *fv = argc > 2 ? fopen(argv[2], "r") : fopen("vectorB.dat", "r"); 
    matrix m = create_matrix_from_file(fm);
    values vector = create_vector_from_file(fv);

    double *x = initialize_x(m.n, 0.0);
    solve(&m, x, vector.v);
    print_vector(x, m.n);
#endif
    return 0;
}
