/*
 * Prim's Algorithm
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ==================== UTILS ======================
void fscanfEdgeList(FILE* file, int **adMatrix, int *nodesNmb) {
  // Function executed on processId = 0

  int edgesNmb, node1, node2, weight, matrixSize, nodes;

  // read nodes number and edges count
  fscanf(file, "nodes-%d,edges-%d,weights", &nodes, &edgesNmb);
  *nodesNmb = nodes;

  // alloc adMatrix as 1D array
  matrixSize = nodes * nodes;
  *adMatrix = (int*) malloc(matrixSize * sizeof(int));
  
  // initialize adMatrix
  for (int i = 0; i < matrixSize; (*adMatrix)[i] = 0, i++);

  // iterate over every edge from csv file
  for (int i = 0; i < edgesNmb; i++) {
      fscanf(file, "%d,%d,%d", &node1, &node2, &weight);
      (*adMatrix)[node1 * nodes + node2] = weight;
      (*adMatrix)[node2 * nodes + node1] = weight;
  }
}

void fprintfAdjMatrix(FILE* file, int* adjMatrix, int nodesNmb) {
  // Function prints adj matrix to output stream

  int index = 0;
  for (int row = 0; row < nodesNmb; row++) {
    for (int column = 0; column < nodesNmb; column++,index++) {
      fprintf(file, "%d ", adjMatrix[index]);
    }
    fprintf(file, "\n");
  }
}

// ==================== ALGORITM ======================


// ====================== MAIN ========================
int main( int argc, char *argv[] )
{
	// New
  FILE  *file;
	char  *filename = argv[1];
  int   *adMatrix = 0;
  int   nodesNmb = 0;

	// Old
  int n, myid, numprocs, i;
  double PI25DT = 3.141592653589793238462643;
  double mysum, pi, sum, inv;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	if (myid == 0) {

    // Read edge list from file
    file = fopen(filename, "r");
    if (file == 0) {
      printf("Cannot found input file %s!\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
		fscanfEdgeList(file, &adMatrix, &nodesNmb);
    fclose(file);

    fprintfAdjMatrix(stdout, adMatrix, nodesNmb);    
    free(adMatrix);
	}

  while (1) {
    if (myid == 0) {

     
      printf("Enter the number of Leibniz series terms: (0 quits) ");
      scanf("%d",&n);
    }

/***
   * Fill in, please ...
   */

      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


      if (n == 0)
        break;

      mysum = sum = inv = pi = 0.0;

      for (i = myid; i < n; i += numprocs) {
        inv = 1.0 / (2.0 * (double) i  + 1.0);
        if (i % 2 == 0)
          mysum += inv;	// i is even
        else
          mysum -= inv;	// i is odd
      }

      MPI_Reduce(&mysum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if (myid == 0) {
        pi = 4 * sum;
        printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
      }
    }
    MPI_Finalize();
    return 0;
}

