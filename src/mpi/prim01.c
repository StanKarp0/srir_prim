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
  
  FILE  *file;
	char  *filename = argv[1];
  int   *adMatrix = 0;
  int   nodesNmb = 0;
  int   processId, processNmb;

  // MPI Initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processNmb);
  MPI_Comm_rank(MPI_COMM_WORLD, &processId);

	if (processId == 0) {

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

  MPI_Finalize();
  return 0;
}

