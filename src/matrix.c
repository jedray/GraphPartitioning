#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include "matrix.h"
#include "vector.h"
#include "mpi.h"


/*****************************************************************************/
/*                             Sequential Algorithms                         */
/*****************************************************************************/


// The 1-norm of the matrix A
int norm1(int **A, int n_rows, int n_colums){
    int i, j, max, sum;
    max = 0;
	// Process p finds the maximum sum each column wise sum
    for(i = 0; i<n_rows; i++){
		sum = 0;
		for( j=0; j<A[i][0]; j++)
			sum += ABS(A[i][2*j+2]);
		if(sum > max)
			max = sum;
	}
    return max;
}


// Creates the matrix  gI - A  where g is a scalar and I is the identity  
int ** m_neg_shift(int **A, int g, int n_rows){
	int i, j;
	int **B = (int **) malloc(n_rows*sizeof(int*));
	for(i=0;i<n_rows;i++)
		B[i] = (int *) calloc(A[i][0]*2 + 1,sizeof(int));
	for(i=0;i<n_rows;i++){
		B[i][0] = A[i][0];
		B[i][2] = g;
		for(j=0;j<A[i][0];j++){
			B[i][2*j+1] = A[i][2*j+1];
			B[i][2*j+2] += - A[i][2*j+2];
		}
	}
	return B;
}


// Display the matrix A
void m_print(int **A, int n_rows, int n_columns){
	int i,j;
	int **B = malloc(n_rows*sizeof(int *));
	for(i=0;i<n_rows;i++){
		B[i] = calloc(n_columns,sizeof(int));
		j = 0;
		while(j < A[i][0]){
			B[i][A[i][2*j+1]] = A[i][2*j+2];
			j++;
		}
	}
	for(i=0;i<n_rows;i++){
		for(j=0;j<n_columns;j++){
			if(B[i][j] >= 0)
                printf("  %d ", B[i][j]);
			else 
				printf(" %d ",B[i][j]);
		}
		printf("\n");
	}
}



// find the Fieldler vector that corresponds to the second eigenvalue 
// knowing the eigenvector that corresponds to the largest eigenvalue
double * power_iteration(int **A, double *X, int n_columns){
    
    int i; 
	double *y ; 
	double *projection;
	double *tmp;
	double *x;
    double c , cnew, error;
    double *u;
    
    // initial guess for the power method
    x = (double*) malloc(n_columns*sizeof(double));
    srand(time(NULL));
    for(i=0; i<n_columns; i++)
		x[i] = SIGN(rand());

    // project on the orthogonal space to X
    u = v_unorm(X,n_columns);      
	projection = v_mult(u,v_dot(x,u,n_columns), n_columns);
	tmp = v_sub(x,projection,n_columns);
	free(x);
	x = tmp;
	free(projection);
	
    // iterate until convergence condition
    i=0;
    c = 0;
    error = 1;
	while( error >= 0.000001 * c ){
        // the unit vector of x
		y = v_unorm(x,n_columns);                      
        free(x);
        // multiply A by x
		x = mv_mult(A,y,n_columns);                 
        // project x to the orthogonal space to X
        projection = v_mult(u,v_dot(x,u,n_columns),n_columns); 
        tmp = v_sub(x,projection,n_columns);           
        free(x);
        x = tmp;                         
        // scalar product of old x and new x, i.e. < Ax|x > 
        cnew = v_dot(x,y,n_columns);  
        // error convergence  
        error = ABS(cnew - c);        
        c = cnew;
        i++;
        free(projection);
        free(y);
    }
    free(u);

    //printf(" PM : number of steps %d \n", i);
    return x;
}


   




/*****************************************************************************/
/*                             Parallel functions                            */ 
/*****************************************************************************/


// All processes receives the 1-norm of the matrix A. All processes must call
// the function.
int par_norm1(int **A, int n_rows, int n_columns, MPI_Comm comm_world){
    int i, j, my_max, max, sum;
    my_max = 0;
	// Process p finds the max sum
    for(i = 0; i<n_rows; i++){
		sum = 0;
		for( j=0; j<A[i][0]; j++)
			sum += ABS(A[i][2*j+2]);
		if(sum > my_max)
			my_max = sum;
	}
    // Processes work together to find the 1-norm of A 
    MPI_Allreduce(&my_max, &max, 1, MPI_INT, MPI_MAX, comm_world);    
    return max;
}


// All processes work together to find the Fiedler vector. Process i will 
// receive the i-th block of the resulting Feidler vector. All processes    
// must call this function.
double * par_power_iteration(int **A, int my_n_rows, int n_columns, MPI_Comm comm_world){

    int i; 
	double *y ; 
	double *projection;
	double *tmp;
	double *x;
    double *X;
    double c , cnew, error;	
    double *u;
    
    // initial guess 
    x = (double*) malloc(my_n_rows*sizeof(double));
    srand(time(NULL));
    for(i=0; i<my_n_rows; i++)
		x[i] = SIGN(rand());

    // the i-th block the eigenvector corresponding to the largest eigenvalue
    X = v_create(1,my_n_rows);
	
    // project x on the orthogonal space to X
    u = par_v_unorm(X, my_n_rows, comm_world);      	
	projection = v_mult(u,par_v_dot(x,u,my_n_rows,comm_world), my_n_rows);
	tmp = v_sub(x,projection,my_n_rows);
	free(x);
	x = tmp;
	free(projection);

    // iterate the power method loop
    error = 1;
    c = 0;
	i=0;
	while( error >= 0.000001 * c ){
        // the unit vector of old x
		y = par_v_unorm(x,my_n_rows, comm_world);     
        free(x);
        // multply A by x 
		x = par_mv_mult(A,y,my_n_rows,n_columns,comm_world);  
        // project new x on the orthogonal space to X
        projection = v_mult(u,par_v_dot(x,u,my_n_rows,comm_world),my_n_rows); 
        tmp = v_sub(x,projection,my_n_rows);                 
        free(x);
        x = tmp;
        // the scalar product of old x and new x
        cnew = par_v_dot(x,y,my_n_rows, comm_world); 
        // the convergence error
        error = ABS(cnew - c);
        c = cnew;
        i++;
        free(projection);
        free(y);
    }
    free(u);

    //printf(" PM : number of steps %d \n", i);
    return x;
}

