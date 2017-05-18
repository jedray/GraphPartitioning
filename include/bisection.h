#ifndef _BISECTION_H
#define _BISECTION_H


/* Merge sort is used to find quantiles */
void merge(double *u, double *t1, int n1, double *t2, int n2);
void merge_sort(double *t, int n);

/* Distributes a matrix by set of rows over the other processes */
void distribute_matrix(int **matrix, int n_rows, int ***my_matrix, 
        int *my_n_rows, MPI_Comm comm_world);

/* Bisect into two partitions */
int *bisect(double *my_fv, int my_n_rows, int n_rows, int k, int *k1, 
        int *k2, int *n1, int *n2, MPI_Comm comm_world);

/* Partition with respect to a given quantile and Feidler vector */
int *do_bisect(double * fv, double quantile, int quantile_count, int n_rows);

/* Spectral bisection */
void spectral_bisection(int **laplace, int *n_rows , int *n_columns, int k, 
        int *k1, int *k2, int *n1, int *n2, int **partition, 
        MPI_Comm comm_world);


/* Perform a k-way partitioning */
void recursive_spectral_bisection(int **laplace, int *n_rows, int *n_columns, 
        int k, int **partition, MPI_Comm comm_world);

/* Given a partition, split the laplacian into two small laplacians and build 
   the corresponding mappings */
void split_laplacian(int **laplace, int *partition, int n_rows, int n1, int n2, 
        int ***laplace1, int***laplace2, int **map1, int **map2);

/* Communicate a matrix between process 0 and process 1 */
void communicate_laplacian(int ***laplace, int n_rows, MPI_Comm comm_world);



#endif
