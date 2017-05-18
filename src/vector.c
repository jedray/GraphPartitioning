#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "vector.h"

/*****************************************************************************/
/*                        Sequential functions                               */
/*****************************************************************************/


// the scalar product of two vectors x and y of legth n 
double v_dot(double* x, double* y, int n){
	double sp = 0;
	int i;
	for(i=0; i<n; i++)
		sp += x[i]*y[i];
	return sp;
}

// the norm of a vector of length n
double v_norm(double *x, int n){	 
	 return sqrt(v_dot(x,x,n)) ; 
}

// Addition of vector x and y of length n
double * v_add(double *x, double *y, int n) {
	double *z = (double *) malloc(n*sizeof(double));
	int i;
	for(i=0; i<n; i++)
		z[i] = y[i] + x[i];
	return z;
}

// Substraction of vector y from vector x of lenght n
double * v_sub(double *x, double *y, int n){
	double *z = (double *) malloc(n*sizeof(double));
	int i;
	for(i=0; i<n; i++)
		z[i] = - y[i] + x[i];
	return z;
}

// Multiply a vector x of length n by a scalar s
double * v_mult(double *x, double s, int n){
	double *z = (double *) malloc(n*sizeof(double));
	int i;
	for(i=0; i<n; i++)
		z[i] = s * x[i];
	return z;
} 

// Divide a vector x of length n by a scalar s 
double * v_div(double *x, double s, int n){	
	return v_mult(x,1.0/s,n);
}


// The unit vector that corresponds to the vector x of legth n
double * v_unorm(double *x, int n){
	return v_div(x,v_norm(x,n),n);
}

// The result of the multiplication of A times x of length n 
double * mv_mult(int **A, double *x, int n){
	double *z = (double *) calloc(n,sizeof(double));
	int i,j;
	for(i = 0; i<n;i++){
		for(j = 0; j<A[i][0]; j++){
			z[i] += A[i][2*j+2] * x[A[i][2*j+1]];
		}
	}
	return z;
}


// Creates a vector of length n with values a
double * v_create(double a, int n){
    double *x = (double *) malloc(n*sizeof(double));
    int i;
    for(i=0;i<n;i++)
        x[i] = a;
    return x;
}

// Display a vector of length n
void v_print(double *x, int n){
	int i;
	for (i = 0; i < n; ++i){
		if(x[i] >= 0)
			printf("  %lf \n", x[i]);
		else 
			printf(" %lf \n",x[i]);
	}
}





/*****************************************************************************/
/*                              Parallel function                            */
/*****************************************************************************/

// The dot product of two vectors X = (x1,..,x,..,xp) and Y=(y1,..,y,..,yp) 
// where X and Y are distributed over p processes. All processes must call 
// the function and they all receive the dot product result of X and Y.
double par_v_dot(double* x, double* y, int n, MPI_Comm comm_world){
	double my_sp = 0, sp;
	int i;
	for(i=0; i<n; i++)
		my_sp += x[i]*y[i];
    MPI_Allreduce(&my_sp, &sp, 1, MPI_DOUBLE, MPI_SUM, comm_world);
	return sp;
}


// The norm of a vectors X = (x1,..,x,..,xp) where X is distributed over p 
// processes. All processes must call the function and they all receive the 
// norm of X.
double par_v_norm(double *x, int n, MPI_Comm comm_world){
    return sqrt(par_v_dot(x,x,n, comm_world));
}


 
// The calling process i will receive the i-th block of the unit vector 
// of the original vector X = (x1,..,xi,..xp). All processes must call the 
// function.  
double * par_v_unorm(double *x, int n, MPI_Comm comm_world){	
	return v_div(x,par_v_norm(x,n, comm_world),n);
}



// The calling process i will receive the i-th block of the result of  
// of the multiplication of A and X = (x1,..,xi,..,xp). All processes   
// must call the function.  
double * par_mv_mult(int **A, double *x, int my_n_rows, int n_rows, 
        MPI_Comm comm_world){

    int world_rank, world_size;
    int *rcv_counts, *displs ; 
    int p,q,r,displ = 0;
    double *X, *y;
    int i,j;
    
    MPI_Comm_rank(comm_world, &world_rank);
    MPI_Comm_size(comm_world, &world_size);
    
    rcv_counts = (int *) calloc(world_size,sizeof(int));
    displs= (int *) calloc(world_size,sizeof(int));    
    X = (double *) malloc(n_rows*sizeof(double));      // result of X = x1,..,xp
    y = (double *) calloc(my_n_rows,sizeof(double));   // result of A*X
    
    q = n_rows / world_size;
    r = n_rows % world_size;

    for(p=0;p<world_size;p++){
        displs[p] = displ;
        rcv_counts[p] = q + (world_size + r-1 - p)/world_size; 
        displ += rcv_counts[p];
    }

    // All processes gathers the vector X = (x1,..,xp)
    MPI_Barrier(comm_world);
    MPI_Allgatherv(x, my_n_rows, MPI_DOUBLE, X, rcv_counts, displs, MPI_DOUBLE, 
            comm_world);
    MPI_Barrier(comm_world);
    
	for(i = 0; i<my_n_rows;i++){
		for(j = 0; j<A[i][0]; j++){
			y[i] += A[i][2*j+2] * X[A[i][2*j+1]];
		}
	}

	return y;
}



