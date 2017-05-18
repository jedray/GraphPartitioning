
/*****************************************************************************/
/*                          Matrix Representation                            */
/*****************************************************************************/

/* 
 *  The adjacency and the laplacian matrices are assumed to be sparse.
 *  therefore they are not stored as n by n matrice but instead as 
 *  n arrays of different lengths.
 * 
 *  The i-th array correspoding to the node i has the form:
 *
 *  |p| |n1|w1| |n2|w2| - - - |np|wp| 
 *
 *  where:   p = # adjacent nodes
 *          nj = j-th adjacent node
 *          wj = weight of edge i--j (non zero) 
 *
 *  example: 
 *
 *  5 by 5 adjacency matrix                   our matrix representation   
 *
 *         0 0 0 1 0                              array 0 : |1| |3|1| 
 *         0 0 1 0 1           <--------->        array 1 : |2| |2|1| |4|1|
 *         0 1 0 0 0           <--------->        array 2 : |1| |1|1|
 *         1 0 0 0 1           <--------->        array 3 : |2| |0|1| |4|1|
 *         0 1 0 1 0                              array 4 : |2| |1|1| |3|1|
 *
**/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "parser.h"


int parse(char *filename, int *n, int *m, int ***adjacency, int ***laplace) {	
	
	int type, i, j, tmp;  
	
	FILE *fp;
	fp = fopen(filename, "r");
	//printf(" %s \n", filename);
	
	char *line = NULL;
	char *linecpy; 
	int length;
	size_t len;
	int num;
	char *token = NULL;
		
	// read number of vertices and edges  
	getline(&line, &len, fp);
	token = strtok(line, " ");
	(*n)  = atoi(token);
	token = strtok(NULL, " ");
	(*m)  = atoi(token);
	token = strtok(NULL, " ");
	if (token != NULL)
		type = atoi(token);
	else 
		type = 0;


	// allocate memory for adjacency matrix
	(*adjacency) = (int **) malloc((*n)*sizeof(int*));
	(*laplace) = (int **) malloc((*n)*sizeof(int*));


	i = 0;
	while (getline(&line,&len,fp) != -1 && i < (*n)){
		// tokenize line
		length = strlen(line) + 1;
		linecpy = (char *) malloc(length*sizeof(char));
		strcpy(linecpy,line);
		// determine number of adjacent nodes
		token = strtok(linecpy, " ");
		int nb = 0;
		while (token != NULL){
			if(atoi(token) > 0)
				nb++;
			token = strtok(NULL, " ");
		}

		// alocate memory for the list of adjacent nodes 
		if(type == 0){
			(*adjacency)[i] = (int *) calloc(1+nb*2,sizeof(int));
			(*laplace)[i] = (int *) calloc(3+nb*2,sizeof(int));
		} else {
			(*adjacency)[i] = (int *) calloc(1+nb,sizeof(int));
			(*laplace)[i] = (int *) calloc(3+nb,sizeof(int));
		}
		
        // # adjacent nodes 
		(*adjacency)[i][0] = (type == 0) ? nb : nb/2;     
        // # adjacent nodes + current node
        (*laplace)[i][0] = (type == 0) ? nb+1 : (nb/2)+1; 
		

		(*laplace)[i][1] = i;                     // identifier of node i for L(i,i)     
		(*laplace)[i][2] = (type == 0) ? nb : 0;  // weight for L(i,i) is deg(i) if unweighted   

		// fill the adjacency lists
		token = strtok(line, " ");
		j=0;
		while (token != NULL && j < nb){
			if(atoi(token) > 0 ){
				(*adjacency)[i][2*j+1] = atoi(token) - 1;
				(*laplace)[i][(2*j+1)+2] = atoi(token) - 1;	
					
				if(type == 0){
					(*adjacency)[i][(2*j+1)+1] = 1;
					(*laplace)[i][(2*j+1)+3] = -1;
				}else{
					token = strtok(NULL, " ");
					if(token == NULL)
						exit(EXIT_FAILURE);
					(*adjacency)[i][(2*j+1)+1] = atoi(token);
					(*laplace)[i][(2*j+1)+3] = -atoi(token);
					(*laplace)[i][2] += atoi(token);
				}				
				j++;
			}
			token = strtok(NULL, " ");
			
		}
		i++;
		free(linecpy);
	}

	fclose(fp); 
	free(line);
    return 1;
}


