#include <stdio.h>
#include "global.h"

const double DELTA = 1e-3;
const int TRUE = 1;
const int FALSE = 0;

int double_equals (double a, double b)
{
    double mod_diff = a < b ? b-a : a-b;
    return mod_diff < DELTA;
}

void mul_matrix_row_ind_test()
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

}

void find_upper_iterative_matrix_test()
{
    printf("\n--------------------------------\n\tSzukanie górnej macierzy do iteracji\n-----------------------\n");
    matrix m;
    double vals4[] = { 3.0, 1.0, 2.0, 1.0, 5.0, 8.0};
    int col_ind4[] = {0, 1, 2, 0, 1, 2};
    int row_start_ind4[] = {0, 1, 3};
    m.vals = vals4;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind4;
    m.row_start_ind = row_start_ind4;
    matrix it_matrix1 = get_iterative_upper_matrix(&m, 1.5);
    printf("m1.vals[0] = %.2f\n", it_matrix1.vals[0]);
    printf("m1.vals[1] = %.2f\n", it_matrix1.vals[1]);
    printf("m1.vals[2] = %.2f\n", it_matrix1.vals[2]);
    printf("m1.vals[3] = %.2f\n", it_matrix1.vals[3]);

    int all_passed = TRUE;
    if (!double_equals(-2.0, it_matrix1.vals[0]))
        all_passed = FALSE;
    if (!double_equals(-2.0, it_matrix1.vals[1]))
        all_passed = FALSE;
    if (!double_equals(-3.0, it_matrix1.vals[2]))
        all_passed = FALSE;
    if (!double_equals(-2.0, it_matrix1.vals[3]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

}

void find_upper_iterative_matrix_simple_test()
{
    printf("\n--------------------------------\n\tSzukanie górnej macierzy do iteracji 2x2\n-----------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;
    matrix it_matrix1 = get_iterative_upper_matrix(&m, 1.0);
    printf("m1.vals[0] = %.2f\n", it_matrix1.vals[0]);
    printf("m1.vals[1] = %.2f\n", it_matrix1.vals[1]);
    printf("m1.vals[2] = %.2f\n", it_matrix1.vals[2]);
    printf("m1.n = %i\n", it_matrix1.n);
    printf("m1.nelements = %i\n", it_matrix1.nelements);
    printf("m1.col_ind[0] = %i\n", it_matrix1.col_ind[0]);
    printf("m1.col_ind[1] = %i\n", it_matrix1.col_ind[1]);
    printf("m1.col_ind[1] = %i\n", it_matrix1.col_ind[2]);

    int all_passed = TRUE;
    if (!double_equals(0.0, it_matrix1.vals[0]))
        all_passed = FALSE;
    if (!double_equals(-0.1875, it_matrix1.vals[1]))
        all_passed = FALSE;
    if (!double_equals(0.0, it_matrix1.vals[2]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

}


void find_lower_iterative_matrix_test()
{
    printf("\n--------------------------------\n\tSzukanie dolnej macierzy do iteracji\n-----------------------\n");
    matrix m;
    double vals4[] = { 3.0, 1.0, 2.0, 1.0, 5.0, 8.0};
    int col_ind4[] = {0, 1, 2, 0, 1, 2};
    int row_start_ind4[] = {0, 1, 3};
    m.vals = vals4;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind4;
    m.row_start_ind = row_start_ind4;
    matrix it_matrix1 = get_iterative_lower_matrix(&m, 1.5);
    printf("m1.vals[0] = %.2f\n", it_matrix1.vals[0]);
    printf("m1.vals[1] = %.2f\n", it_matrix1.vals[1]);
    printf("m1.vals[2] = %.2f\n", it_matrix1.vals[2]);
    printf("m1.vals[3] = %.2f\n", it_matrix1.vals[3]);
    printf("m1.vals[4] = %.2f\n", it_matrix1.vals[4]);

    int all_passed = TRUE;
    if (!double_equals(2.5, it_matrix1.vals[0]))
        all_passed = FALSE;
    if (!double_equals(2.5, it_matrix1.vals[1]))
        all_passed = FALSE;
    if (!double_equals(0.1875, it_matrix1.vals[2]))
        all_passed = FALSE;
    if (!double_equals(0.9375, it_matrix1.vals[3]))
        all_passed = FALSE;
    if (!double_equals(2.5, it_matrix1.vals[4]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

}

void find_lower_iterative_matrix_simple_test()
{
    printf("\n--------------------------------\n\tSzukanie dolnej macierzy do iteracji 2x2\n-----------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;
    matrix it_matrix1 = get_iterative_lower_matrix(&m, 1.0);
    printf("m1.vals[0] = %.2f\n", it_matrix1.vals[0]);
    printf("m1.vals[1] = %.2f\n", it_matrix1.vals[1]);
    printf("m1.vals[2] = %.2f\n", it_matrix1.vals[2]);
    printf("m1.n = %i\n", it_matrix1.n);
    printf("m1.nelements = %i\n", it_matrix1.nelements);
    printf("m1.col_ind[0] = %i\n", it_matrix1.col_ind[0]);
    printf("m1.col_ind[1] = %i\n", it_matrix1.col_ind[1]);
    printf("m1.col_ind[2] = %i\n", it_matrix1.col_ind[2]);

    int all_passed = TRUE;
    if (!double_equals(1.0, it_matrix1.vals[0]))
        all_passed = FALSE;
    if (!double_equals(-0.63636, it_matrix1.vals[1]))
        all_passed = FALSE;
    if (!double_equals(1.0, it_matrix1.vals[2]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

}


void find_iterative_vector_test()
{
    printf("\n-------------------------------\n\tSzukanie wektora do iteracji wD^-1*b\n------------------------\n");
    matrix m;
    double vals4[] = { 3.0, 1.0, 2.0, 1.0, 5.0, 8.0};
    int col_ind4[] = {0, 1, 2, 0, 1, 2};
    int row_start_ind4[] = {0, 1, 3};
    m.vals = vals4;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind4;
    m.row_start_ind = row_start_ind4;
    double b[] = {1.0, 2.0, 3.0};
    double *it_vec = get_iterative_vector(&m, 1.5, b);
    printf("it_vec[0] = %.2f\n", it_vec[0]);
    printf("it_vec[1] = %.2f\n", it_vec[1]);
    printf("it_vec[2] = %.2f\n", it_vec[2]);

    int all_passed = TRUE;
    if (!double_equals(0.5, it_vec[0]))
        all_passed = FALSE;
    if (!double_equals(3.0, it_vec[1]))
        all_passed = FALSE;
    if (!double_equals(0.5625, it_vec[2]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");
}

void find_iterative_vector_simple_test()
{
    printf("\n-------------------------------\n\tSzukanie wektora do iteracji wD^-1*b 2x2\n------------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {11.0, 13.0};
    double *it_vec = get_iterative_vector(&m, 1.0, b);
    printf("it_vec[0] = %.2f\n", it_vec[0]);
    printf("it_vec[1] = %.2f\n", it_vec[1]);

    int all_passed = TRUE;
    if (!double_equals(0.6875, it_vec[0]))
        all_passed = FALSE;
    if (!double_equals(-1.1818, it_vec[1]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");
}

void find_z_simple()
{
    printf("\n---------------------------------\n\tSzukanie z(1) 2x2\n----------------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {11.0, 13.0};
    double x[] = {1.0, 1.0};
    double z[2];
    double w = 1.0;

    matrix upper_iterative = get_iterative_upper_matrix(&m, w);
    double *iterative_vector = get_iterative_vector(&m, w, b);
    mul_matrix_row(&upper_iterative, x, z, 0, 2);
    add_vector(z, iterative_vector, m.n);

    printf("z[0] = %.2f\n", z[0]);
    printf("z[1] = %.2f\n", z[1]);

    int all_passed = TRUE;
    if (!double_equals(0.5, z[0]))
        all_passed = FALSE;
    if (!double_equals(-1.1818, z[1]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");
}

void find_next_x()
{
    printf("\n---------------------------------\n\tSzukanie x(1) 2x2\n----------------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {11.0, 13.0};
    double x[] = {1.0, 1.0};
    double z[2];
    double y[2];
    double w = 1.0;

    matrix upper_iterative = get_iterative_upper_matrix(&m, w);
    matrix lower_iterative = get_iterative_lower_matrix(&m, w);
    double *iterative_vector = get_iterative_vector(&m, w, b);
    mul_matrix_row(&upper_iterative, x, z, 0, 2);
    add_vector(z, iterative_vector, m.n);
    forward_subst(&lower_iterative, y, z);

    printf("y[0] = %.2f\n", y[0]);
    printf("y[1] = %.2f\n", y[1]);
    printf("z[0] = %.2f\n", z[0]);
    printf("z[1] = %.2f\n", z[1]);
    printf("lower_iterative[0] = %.2f\n", lower_iterative.vals[0]);
    printf("lower_iterative[1] = %.2f\n", lower_iterative.vals[1]);
    printf("lower_iterative[2] = %.2f\n", lower_iterative.vals[2]);
    printf("lower_iterative.n = %i\n", lower_iterative.n);
    printf("lower_iterative.nelements = %i\n", lower_iterative.nelements);
    printf("lower_iterative.col_ind[0] = %i\n", lower_iterative.col_ind[0]);
    printf("lower_iterative.col_ind[1] = %i\n", lower_iterative.col_ind[1]);
    printf("lower_iterative.col_ind[2] = %i\n", lower_iterative.col_ind[2]);

    int all_passed = TRUE;
    if (!double_equals(0.5, y[0]))
        all_passed = FALSE;
    if (!double_equals(-0.86364, y[1]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");
}

void solve_simple_system()
{
    printf("\n---------------------------------\n\tPrzykładowy układ 2x2\n----------------------------\n");
    matrix m;
    double vals[] = {16.0, 3.0, 7.0, -11.0};
    int col_ind[] = {0, 1, 0, 1};
    int row_start_ind[] = {0, 2};
    m.vals = vals;
    m.n = 2;
    m.nelements = 4;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {11.0, 13.0};
    double x[] = {1.0, 1.0};
    solve(&m, x, b);
    printf("x[0] = %.2f\n", x[0]);
    printf("x[1] = %.2f\n", x[1]);
}

void solve_system()
{
    printf("\n---------------------------------\n\tPrzykładowy układ\n----------------------------\n");
    matrix m;
    double vals[] = { 3.0, 1.0, 2.0, 1.0, 5.0, 8.0};
    int col_ind[] = {0, 1, 2, 0, 1, 2};
    int row_start_ind[] = {0, 1, 3};
    m.vals = vals;
    m.n = 3;
    m.nelements = 6;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {1.0, 2.0, 3.0};
    double x[] = {1.0, 1.0, 1.0};
    double z[3];
    double y[3];
    double w = 1.5;

    matrix lower_iterative = get_iterative_lower_matrix(&m, w);
    matrix upper_iterative = get_iterative_upper_matrix(&m, w);
    double *iterative_vector = get_iterative_vector(&m, w, b);
    mul_matrix_row(&upper_iterative, x, z, 0, 3);
    add_vector(z, iterative_vector, m.n);
    forward_subst(&lower_iterative, y, z);

    printf("y[0] = %.2f\n", y[0]);
    printf("y[1] = %.2f\n", y[1]);
    printf("y[2] = %.2f\n", y[2]);
}

void solve_system_hard()
{
    printf("\n---------------------------------\n\tPrzykładowy układ 4x4\n----------------------------\n");
    matrix m;
    double vals[] = {10.0, -1.0, 2.0, -1.0, 11.0, -1.0, 3.0, 2.0, -1.0, 10.0, -1.0, 3.0, -1.0, 8.0};
    int col_ind[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3};
    int row_start_ind[] = {0, 3, 7, 11};
    m.vals = vals;
    m.n = 4;
    m.nelements = 14;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {6.0, 25.0, -11.0, 15.0};
    double x[] = {0.0, 0.0, 0.0, 0.0};
    double z[4];
    double y[4];
    double w = 1.0;

    matrix lower_iterative = get_iterative_lower_matrix(&m, w);
    matrix upper_iterative = get_iterative_upper_matrix(&m, w);
    double *iterative_vector = get_iterative_vector(&m, w, b);
    mul_matrix_row(&upper_iterative, x, z, 0, 3);
    add_vector(z, iterative_vector, m.n);
    forward_subst(&lower_iterative, y, z);

    printf("y[0] = %.4f\n", y[0]);
    printf("y[1] = %.4f\n", y[1]);
    printf("y[2] = %.4f\n", y[2]);
    printf("y[3] = %.4f\n", y[3]);
}

void forward_subst_simple_test()
{
    printf("\n---------------------------------\n\tPodstawianie w przód 2x2\n----------------------------\n");
    matrix m;
    double vals[] = {1.0, -0.63636, 1.0};
    int col_ind[] = {0, 0, 1};
    int row_start_ind[] = {0, 1};
    m.vals = vals;
    m.n = 2;
    m.nelements = 3;
    m.col_ind = col_ind;
    m.row_start_ind = row_start_ind;

    double b[] = {0.5, -1.1818};
    double x[2];

    forward_subst(&m, x, b);

    printf("y[0] = %.2f\n", x[0]);
    printf("y[1] = %.2f\n", x[1]);

    int all_passed = TRUE;
    if (!double_equals(0.5, x[0]))
        all_passed = FALSE;
    if (!double_equals(-0.86364, x[1]))
        all_passed = FALSE;
    if (!all_passed)
        printf("Coś nie tak\n");
    else
        printf("\n----------------------------\n\tWszystko w porządku\n----------------------------\n");

}

void test()
{
    matrix m;
    int all_passed = TRUE;
    mul_matrix_row_ind_test();

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

    find_upper_iterative_matrix_test();
    find_lower_iterative_matrix_test();
    find_iterative_vector_test();
    solve_system();
    solve_simple_system();
    find_upper_iterative_matrix_simple_test();
    find_lower_iterative_matrix_simple_test();
    find_iterative_vector_simple_test();
    find_z_simple();
    find_next_x();
    forward_subst_simple_test();
    solve_system_hard();
}
