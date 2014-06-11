#include <stdio.h>
#include <stdlib.h>

#include "global.h"

#define EL_SIZE 100000

typedef struct _values {
    double *v;
    int n;
} values;

typedef struct _indices {
    int *ind;
    int n;
} indices;

typedef struct _element {
    double v;
    int col_ind;
    int row_ind;
} element;

values read_vals(FILE *in)
{
    double vals[EL_SIZE];
    int nel = 0;
    while (fgetc(in) != ']')
    {
        fscanf(in, "%lf", vals+nel);
        nel++;
    }
    nel--;

    values v;
    v.v = vals;
    v.n = nel;
    return v;
}

indices read_indices(FILE *in)
{
    int *v = malloc(EL_SIZE*sizeof(*v));
    int nel = 0;
    while (fgetc(in) != ']')
    {
        fscanf(in, "%i", v+nel);
        v[nel]--;
        nel++;
    }
    nel--;
    indices ind;
    ind.ind = v;
    ind.n = nel;
    return ind;
}

void omit_beginning(FILE *in)
{
    while (fgetc(in) != '[')
        ;
}

indices assign_col_ind_to_each(int nel, indices col_start)
{
    indices col_ind;
    col_ind.ind = malloc(nel*sizeof(int));
    col_ind.n = nel;
    int column_index = 0;
    for (int i = 0; i < nel; i++)
    {
        if (column_index < nel-1)
        {
            if (i == col_start.ind[column_index+1])
                column_index++;
        }

        col_ind.ind[i] = column_index;
    }
    return col_ind;
}

int compare_elements(const void *el1, const void *el2)
{
    element *e1 = (element*) el1;
    element *e2 = (element*) el2;
    if (e1->row_ind == e2->row_ind)
        return e1->col_ind > e2->col_ind;
    else
        return e1->row_ind > e2->row_ind;
}

indices compress_row_ind(int *row_ind, int nel, int nrow)
{
    int *row_start = malloc(nrow*sizeof(*row_start));
    int prev_row_ind = row_ind[0];
    int ri = 1;
    row_start[0] = 0;
    for (int i = 1; i < nel; i++)
    {
        if (row_ind[i] != prev_row_ind)
        {
            row_start[ri++] = i;
            prev_row_ind = row_ind[i];
        }
    }
    indices row_start_ind;
    row_start_ind.n = nrow;
    row_start_ind.ind = row_start;
    return row_start_ind;
}

matrix create_matrix_from_file(FILE *in)
{
    omit_beginning(in);
    values v = read_vals(in);
    omit_beginning(in);
    indices row_ind = read_indices(in);
    omit_beginning(in);
    indices col_start = read_indices(in);
    indices col_ind = assign_col_ind_to_each(v.n, col_start);

    element *ccs_matrix = malloc(v.n*sizeof(*ccs_matrix));
    for (int i = 0; i < v.n; i++)
    {
        ccs_matrix[i].v = v.v[i];
        ccs_matrix[i].col_ind = col_ind.ind[i];
        ccs_matrix[i].row_ind = row_ind.ind[i];
    }
    qsort(ccs_matrix, v.n, sizeof(*ccs_matrix), compare_elements);

    double *new_val = malloc(v.n*sizeof(*new_val));
    int *new_col_ind = malloc(row_ind.n*sizeof(*new_col_ind));
    int *new_row_ind = malloc(row_ind.n*sizeof(*new_row_ind));
    for (int i = 0; i < v.n; i++)
    {
        new_val[i] = ccs_matrix[i].v;
        new_col_ind[i] = ccs_matrix[i].col_ind;
        new_row_ind[i] = ccs_matrix[i].row_ind;
    }
    indices row_start_ind = compress_row_ind(new_row_ind, v.n, col_start.n);

    matrix m;
    m.vals = new_val;
    m.col_ind = new_col_ind;
    m.row_start_ind = row_start_ind.ind;
    m.nelements = v.n;
    m.n = col_start.n;
    return m;
}

matrix create_matrix_from_file_nieuzywane(FILE *in)
{
    omit_beginning(in);
    values v = read_vals(in);
    omit_beginning(in);
    indices row_ind = read_indices(in);
    omit_beginning(in);
    indices col_start = read_indices(in);

    double new_values[EL_SIZE];
    int new_row_ind[EL_SIZE];
    // int new_col_start[EL_SIZE];

    for (int i = 0; i < col_start.n-1; i++)
    {
        for (int j = col_start.ind[i]; j < col_start.ind[i+1]; j++)
        {
            int new_el_ind = col_start.ind[row_ind.ind[j]];
            new_values[new_el_ind] = v.v[j];
            new_row_ind[new_el_ind] = i;
        }
    }

    matrix m;
    m.vals = new_values;
    m.col_ind = new_row_ind;
    // m.row_start_ind = new_col_start
    m.nelements = v.n;
    m.n = col_start.n;
    return m;
}
