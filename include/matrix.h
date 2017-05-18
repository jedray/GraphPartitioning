#ifndef _MATRIX_H
#define _MATRIX_H

#include "mpi.h" 

#define SIGN(x) (((x)%2) == 0 ? 1: -1)
#define ABS(x) ((x) > 0 ? (x) : (-(x))) 

// sequential functions
int norm1(int **A, int n_rows, int n_columns);
int ** m_neg_shift(int **A, int g, int n_rows);
void m_print(int **A, int n_rows, int n_columns);
double * power_iteration( int **A, double *X, int n);


// parallel functions
int par_norm1(int **A, int n_rows, int n_columns, MPI_Comm comm_world);
double * par_power_iteration(int **A, int my_n_rows, int n_columns, MPI_Comm comm_world);


#endif
