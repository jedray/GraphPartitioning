#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mpi.h"
#include "matrix.h"
#include "vector.h"
#include "bisection.h"


// sorts a vector using merge sort algorihtm
void merge_sort(double *t , int n){
    if(n>1){
        /* divide and sort*/
        int n1,n2;
        double *t1, *t2;
        n1 = n/2;
        n2 = n - n1;
        t1 = t;
        t2 = t + n1;
        merge_sort(t1, n1);
        merge_sort(t2, n2);

        /* merge */
        double *u = (double *) malloc(n*sizeof(double));
        int k;
        
        merge(u,t1,n1,t2,n2);
        for(k=0; k<n; k++){
            t[k] = u[k];
        }
        free(u);
    }
}

void merge(double *u, double *t1, int n1, double*t2, int n2){
    int n, k, i, j;
    n = n1 + n2;
    i=0;
    j=0;
    for(k=0; k<n; k++){
        if (i >= n1){
            u[k] = t2[j];
            j++;
        } else if (j>= n2){
            u[k] = t1[i];
            i++;
        } else if (t1[i] > t2[j]) {
            u[k] = t2[j];
            j++;
        } else {
            u[k] = t1[i];
            i++;
        }
    }
}


// Distributes matrix over all processes of comm_world
void distribute_matrix(int **matrix, int n_rows, int ***my_matrix, 
        int *my_n_rows, MPI_Comm comm_world){
   
    // only process 0 knows matrix 
    
    int i, j, q, r;
    int world_size, world_rank, tag;
    MPI_Status status;
    int global_index;
    int local_n_rows;
    int local_n_columns;
       
    MPI_Comm_size(comm_world, &world_size);
    MPI_Comm_rank(comm_world, &world_rank);

    q = n_rows / world_size;
    r = n_rows % world_size;
       
    tag = 100;

    if(world_rank == 0){

        // process 0 sends chunks of the matrix to other processes 
        global_index = 0;
        for(i=1; i<world_size; i++){
            // # rows of the chunk to send 
            local_n_rows = q + (world_size +r-1-i)/world_size;
            MPI_Send(&local_n_rows, 1, MPI_INT, i, tag, comm_world);
            // position of the chunk
            global_index += q + (world_size +r-i)/world_size;
            for(j=0;j<local_n_rows;j++){
                // # columns of each row of the chunk
                local_n_columns = 2*(matrix[global_index+j][0])+1;
                MPI_Send(&local_n_columns, 1, MPI_INT, i, tag+1, comm_world);
                MPI_Send(matrix[global_index+j], local_n_columns, MPI_INT, i, 
                        tag+2, comm_world);
            }
        }
        
        // process 0 initializes its own matrix        
        // # rows of the chunk to save  
        *my_n_rows = q + (world_size+r-1)/world_size;
        *my_matrix = (int **) malloc((*my_n_rows)*sizeof(int *));        
        for(i = 0; i<*my_n_rows; i++){
            // # columns of each row in the chunk
            local_n_columns = 2*(matrix[i][0])+1;
            (*my_matrix)[i] = (int *) malloc(local_n_columns*sizeof(int));
            for( j=0; j<local_n_columns; j++)
               (* my_matrix)[i][j] = matrix[i][j];
        }

    } else {

        // process world_rank receives # rows of its chunk
        MPI_Recv(my_n_rows, 1, MPI_INT, 0, tag, comm_world, &status);
        (*my_matrix) = (int **) malloc((*my_n_rows)*sizeof(int*));
        for(j=0;j<*my_n_rows;j++){
            // process world_rank receives # columns of each row in its chunk
            MPI_Recv(&local_n_columns, 1, MPI_INT,0, tag+1, comm_world, &status);
            (*my_matrix)[j] = (int *) calloc(local_n_columns,sizeof(int));
            // process world_rank receives each row of its chunk 
            MPI_Recv((*my_matrix)[j], local_n_columns, MPI_INT, 0, tag+2, 
                    comm_world, &status);
        }
    }
   
}


// Processes give each a block of the Fiedler vector to process 0. Process 0 
// uses the complete Fiedler vector to construct a partition 
int *bisect(double *my_fv, int my_n_rows, int n_rows, int k, int *k1, int *k2, 
        int *n1, int *n2, MPI_Comm comm_world){
    
    int world_rank, world_size;
    int r,q;
    double *fv, *fv_sorted;
    int * rcv_counts ; 
    int * displs ; 
    int p,displ;
    int quantile_count; 
    double quantile;
    int *bisection = NULL;
        
    MPI_Comm_rank(comm_world, &world_rank);
    MPI_Comm_size(comm_world, &world_size);    
        
    if (world_rank == 0) 
        fv= (double *) malloc(n_rows*sizeof(double));
    
    rcv_counts= malloc(world_size*sizeof(int));
    displs= malloc(world_size*sizeof(int));
    
    q = n_rows / world_size;
    r = n_rows % world_size;
    displ = 0;
    
    for(p=0; p<world_size; p++){
        displs[p] = displ;
        rcv_counts[p] = q + (world_size + r- 1 -p)/world_size; 
        displ += rcv_counts[p];
    }
    
    // Process 0 gathers the parts and forms the field vector
    MPI_Barrier(comm_world);    
    MPI_Gatherv(my_fv,my_n_rows,MPI_DOUBLE,fv,rcv_counts,displs,
            MPI_DOUBLE,0,comm_world);    
    MPI_Barrier(comm_world);
    
    // determine k1, k2
    *k1 = k / 2;
    *k2 = k - (*k1);
    *n1 = ((*k1)*n_rows) / k;
    *n2 = n_rows - (*n1);
    // Process 0 finds the bisection
    if (world_rank == 0){
        
        // Process 0 sorts the Fiedler vector
        fv_sorted= v_mult(fv,1,n_rows);
        merge_sort(fv_sorted,n_rows);
        
        // Process 0 finds the quantile with respect to k
        quantile_count = *n1;
        quantile = fv_sorted[quantile_count];
        
        // Process 0 makes the partition 
        bisection = do_bisect(fv,quantile, quantile_count, n_rows);
    }

    return bisection;
}


// Given the Feidler vector, a quantile the function partitions into two sets 
// the nodes 
int * do_bisect(double * fv, double quantile, int quantile_count, int n_rows){
    int * bisection; 
    int i, count = 0;
    bisection = (int *) calloc(n_rows, sizeof(int));
    for(i=0; i<n_rows; i++){
        if ( fv[i] <= quantile && count < quantile_count){
            bisection[i] = 0;
            count++;
        }else{
            bisection[i] = 1;
        }
    }
    assert(count == quantile_count);
    return bisection;
}

// This is the so called spectral bisection
void spectral_bisection(int **laplace, int *n_rows , int *n_columns, int k, 
        int *k1, int *k2, int *n1, int *n2, int **partition, 
        MPI_Comm comm_world){
  
    // Only process 0 knows laplace, n_rows, n_columns
    
    int **my_laplace;
    int my_n_rows;
    double *my_fv;
    int g;
    int **my_B;
    int i;
    
    //printf("check broadcasting \n");
    // Process 0 broadcasts the dimensions of the laplacian matrix 
    MPI_Bcast(n_columns, 1, MPI_INT, 0, comm_world);
    MPI_Bcast(n_rows, 1, MPI_INT, 0, comm_world);    

    //printf("check distributing \n");
    // Process 0 distributes the laplacian matrix over the processes
    distribute_matrix(laplace,*n_rows,&my_laplace,&my_n_rows,comm_world);
    //printf("my_n_rows %d n_rows %d \n", my_n_rows,*n_rows);
     	
    //printf("check par norm1 \n");
    // Processes work togather to find a bound on the max eigenvalue
    g = par_norm1(my_laplace, my_n_rows, *n_rows, comm_world);  
 
    //printf("check neg shift \n");
    // Process p apply a transformation on its part of the laplacian
    my_B = m_neg_shift(my_laplace,g,my_n_rows);
   
    //printf("check power iteration \n");
    // Processes work together to find the Fiedler vector 
    my_fv = par_power_iteration(my_B,my_n_rows,*n_rows, comm_world);
    
    //printf("check bisect \n");
    // Process 0 obtaines the bisection
    (*partition) = bisect(my_fv, my_n_rows, *n_rows, k, k1, k2, n1, n2, 
            comm_world);
    
    // free allocated memory
    free(my_fv);
    for(i=0;i<my_n_rows;i++){
        free(my_B[i]);
        free(my_laplace[i]);
    }
    free(my_B);
    free(my_laplace);
}



void split_laplacian(int **laplace, int *partition, int n_rows, int n1, int n2, 
        int ***laplace1, int***laplace2, int **map1, int **map2){

    (*laplace1) = (int**) malloc(n1*sizeof(int*));
    (*map1) = (int *) calloc(n1,sizeof(int));
    
    (*laplace2) = (int**) malloc(n2*sizeof(int*));
    (*map2) = (int *) calloc(n2,sizeof(int));

    int *map = (int *) calloc(n_rows,sizeof(int));
    int local_length,i, j, k, l, r, gain;
    l = 0;
    r = 0;
    for(i=0;i<n_rows;i++){
        if(partition[i] == 0){
            map[i] = l;
            l++;
        }else{
            map[i] = r;
            r++;
        }
    }
    l=0;
    r=0;
    for(i=0;i<n_rows;i++){
        local_length= 0;        
        if(partition[i] == 0){
            (*map1)[l] = i;
            gain = 0;
            for(j=0;j<laplace[i][0];j++){
                if(partition[laplace[i][2*j+1]] == 0){
                    local_length++;
                }else{
                    gain += laplace[i][2*j+2]; 
                }
            }
            (*laplace1)[l] = malloc((2*local_length+1)*sizeof(int));
            (*laplace1)[l][0] = local_length;
            k = 0;
            for(j=0;j<laplace[i][0];j++){
                if(partition[laplace[i][2*j+1]] == 0){
                    (*laplace1)[l][2*k+1] = map[laplace[i][2*j+1]];
                    (*laplace1)[l][2*k+2] = laplace[i][2*j+2];
                    if(l == map[laplace[i][2*j+1]]){
                        (*laplace1)[l][2*k+2] += gain;
                    }
                    k++;
                }
            }
           l++;
        } else {
            (*map2)[r] = i;
            gain = 0;
            for(j=0;j<laplace[i][0];j++){
                if(partition[laplace[i][2*j+1]] == 1){
                    local_length++;
                }else{
                    gain += laplace[i][2*j+2]; 
                }
            }

           // printf("local lenght 2  %d \n",local_length);
            (*laplace2)[r] = malloc((2*local_length+1)*sizeof(int));
            (*laplace2)[r][0] = local_length;
            k = 0;
            for(j=0;j<laplace[i][0];j++){
                if(partition[laplace[i][2*j+1]] == 1){
                    (*laplace2)[r][2*k+1] = map[laplace[i][2*j+1]];
                    (*laplace2)[r][2*k+2] = laplace[i][2*j+2];
                    if(r== map[laplace[i][2*j+1]]){
                        (*laplace2)[r][2*k+2] += gain;
                    }
                    k++;
                }
            }   
            r++;
        }
    }
    assert(l == n1);
    assert(r == n2);
}

// use a spectral bisection to split the problem into two subproblems 
void recursive_spectral_bisection(int **laplace, int *n_rows, int *n_columns, 
        int k, int **partition, MPI_Comm comm_world){
    
    // if k = 1 then no partition is needed
    if (k <= 1){
        (*partition) = (int *) calloc(*n_rows,sizeof(int));
    // if k > 1 then do partition
    } else {

        int k1, k2, n1, n2, i;
        int world_rank, world_size, tag, color;
        int **laplace1, **laplace2;
        int *partition1, *partition2;
        int *map1, *map2;
        MPI_Status status;
        MPI_Comm sub_comm_world;
        
        
        partition1 = NULL;
        partition2 = NULL;
        laplace1 = NULL;
        laplace2 = NULL;
        tag = 400;
        
        MPI_Comm_rank(comm_world, &world_rank);
        MPI_Comm_size(comm_world, &world_size);


        // find a partition to the current graph
        spectral_bisection(laplace, n_rows, n_columns, k, &k1, &k2, &n1, &n2,
                partition, comm_world);
        assert(k1 + k2 == k);
        assert(n1 + n2 == *n_rows);

        // split the graph problem into two graph subproblems given the partition        
        if(world_rank == 0){
            // Process 0, given a partition, splits the graph into two small graphs
            // thus, two laplacians and two mappings of nodes 
            split_laplacian(laplace,*partition,*n_rows,n1,n2,&laplace1,&laplace2,
                    &map1, &map2);
        }
        
        // solve the two graph subproblems    
        if(world_size <= 1){
            // Process 0 is the only available process
            // Process 0 solves both the two graph subproblems
            recursive_spectral_bisection(laplace1,&n1,&n1,k1,&partition1,
                    comm_world);
            recursive_spectral_bisection(laplace2,&n2,&n2,k2,&partition2,
                    comm_world);
        }else{
            // Many processes are available
            // Split into two group of processes
            // Old process 0 is the new process 0 of the first group of processes
            // Old process 1 is the new process 0 of the second group of processes
            color = world_rank % 2;
            MPI_Comm_split(comm_world,color,world_rank,&sub_comm_world);
            
            // Process 0 sends laplace 2 to process 1
            if(world_rank <= 1){
                communicate_laplacian(&laplace2,n2,comm_world);
            }
            
            if(color == 0){
                // First group of processes solves the first graph subproblem
                recursive_spectral_bisection(laplace1,&n1,&n1,k1,&partition1,
                            sub_comm_world);
            }
            if(color == 1){
                // Second group of processes solves the second graph subproblem
                recursive_spectral_bisection(laplace2,&n2,&n2,k2,&partition2,
                        sub_comm_world);
            }
                
            if(world_rank == 1)
                MPI_Send(partition2,n2,MPI_INT,0,tag,comm_world);
            
            if(world_rank == 0){
                // Process 0 gathers the solution of both groups
                partition2 = (int *) calloc(n2,sizeof(int));
                MPI_Recv(partition2,n2,MPI_INT,1,tag,comm_world,&status);
            }
            
            
            
        }
        
        if(world_rank == 0){
            for(i=0;i<n1;i++)
                (*partition)[map1[i]] = partition1[i];
            for(i=0;i<n2;i++)
                (*partition)[map2[i]] = partition2[i] + k1;
        }
    }
}

void communicate_laplacian(int ***laplace, int n_rows, MPI_Comm comm_world){
     
   int world_rank, world_size, tag, my_n_columns, i; 
   MPI_Status status;
   MPI_Comm_rank(comm_world, &world_rank);
   MPI_Comm_size(comm_world, &world_size);    
   tag = 200; 
   if(world_size > 1){
       // Process 0 sends the laplace matrix
       if(world_rank == 0){
           for(i=0;i<n_rows;i++){
               my_n_columns = (*laplace)[i][0]*2 + 1;
               MPI_Send(&my_n_columns,1,MPI_INT,1,tag+2*i,comm_world);
               MPI_Send((*laplace)[i],my_n_columns,MPI_INT,1,tag+2*i+1,comm_world);
           }
       }
       // Process 1 receives the laplace matrix 
       if(world_rank == 1){
           (*laplace) = (int **) malloc(n_rows*sizeof(int*));
           for(i=0;i<n_rows;i++){
               MPI_Recv(&my_n_columns,1,MPI_INT,0,tag+2*i,comm_world,&status);
               (*laplace)[i] = (int *) calloc(my_n_columns,sizeof(int));
               MPI_Recv((*laplace)[i],my_n_columns,MPI_INT,0,tag+2*i+1,
                       comm_world,&status);
           }
       }
   } 
}
