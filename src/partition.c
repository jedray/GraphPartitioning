#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mpi.h"
#include "bisection.h"
#include "parser.h"
#include "vector.h"
#include "matrix.h"


int main(int argc, char *argv[]){

    int rank, size;
    MPI_Status status; 
    int n,m;
    int **adjacency = NULL;
    int **laplace = NULL;    
    int *partition;
    int i, k;
    double t1, t2;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0 && argc != 2){
		printf( "[USAGE]: ./do_par filename \n");
		exit(EXIT_FAILURE);
	}	

    MPI_Barrier(MPI_COMM_WORLD);
    t1= MPI_Wtime();
    
    if(rank == 0){
        parse(argv[1],&n,&m,&adjacency, &laplace);
    }
    
    k = 2;
    recursive_spectral_bisection(laplace,&n,&n,k, &partition, MPI_COMM_WORLD);
    
   /* 
    if(rank == 0){
        for(i=0;i<n;i++){
            printf("%d \n",partition[i]);
        }
        printf("\n");
    }
    */
    t2 = MPI_Wtime();    
    MPI_Barrier(MPI_COMM_WORLD);
    
    printf("process %d : time spent partitioning %f \n",rank, t2-t1);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(size-1 == rank)
        printf("\n ========================================= \n");
    MPI_Finalize();    
    return(EXIT_SUCCESS);
}

