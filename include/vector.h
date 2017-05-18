#ifndef _VECTOR_H
#define _VECTOR_H

#include "mpi.h"

/* Sequential functions */

double v_dot(double* x, double* y, int n);
double v_norm(double *x, int n);

double * v_add(double *x, double *y, int n);
double * v_sub(double *x, double *y, int n);
double * v_mult(double *x, double s, int n);
double * v_div(double *x, double s, int n);
double * v_unorm(double *x, int n);

double * mv_mult(int **A,  double *x, int n);
double * v_create(double a, int n);

void v_print(double *x, int n);

// parallel functinos 
double par_v_dot(double *x, double *y, int n, MPI_Comm comm_world);
double par_v_norm(double *x, int n, MPI_Comm comm_world);
double * par_v_unorm(double *x, int n, MPI_Comm comm_world);
double * par_mv_mult(int **A, double *x, int my_n_rows, int n_rows, 
        MPI_Comm comm_world);

#endif
