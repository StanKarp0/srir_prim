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

  // allocation of adMatrix as 1D array
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

void fprintfAdMatrix(FILE* file, int* adMatrix, int rowNmb, int colNmb) {
  // Function prints adj matrix to output stream

  int index = 0;
  for (int row = 0; row < rowNmb; row++) {
    for (int column = 0; column < colNmb; column++,index++) {
      fprintf(file, "%d ", adMatrix[index]);
    }
    fprintf(file, "\n");
  }
}

// ==================== ALGORITM ======================
void primPartitionMatrix(int *adMatrixFull, int **adMatrixPartial, int nodesNmb, int *nodesProcesNmb, int processId, int processNmb) {
  // The partitioning of the adjecency matrix among processNmb processes. Page 445.
  // Last partition will be different size then others. nodes % processes can be different from 0.

  int matrixSize = nodesNmb * nodesNmb;
  int lastId = processNmb - 1;

  // node nuber
  int middleNodes = (int) ceil((float) nodesNmb / (float) processNmb);
  int middleSize = nodesNmb * middleNodes;

  // size in 1d array indexes
  int lastNodes = nodesNmb - (middleNodes * lastId);
  int lastSize = matrixSize - lastId * middleSize;
  int lastCommProcessNmb, lastCommProcessId;

  // Separate communicator to last process - different matrix size
  // Even when lastSize == middleSize logic stay the same for simplicity
  MPI_Comm lastComm;
  
  if (processNmb > 1) {
    // more then one procces - partition

    // int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
    //   Creates new communicators based on colors and keys 
    // Input Parameters
    //   comm - communicator (handle) 
    //   color - control of subset assignment (nonnegative integer). 
    //     Processes with the same color are in the same new communicator 
    //   key - control of rank assignment (integer) 
    // Output Parameters
    //   newcomm - new communicator (handle) 
    MPI_Comm_split(MPI_COMM_WORLD, processId / lastId, processId, &lastComm);
    MPI_Comm_size(lastComm, &lastCommProcessNmb);
    MPI_Comm_rank(lastComm, &lastCommProcessId);

    // lastComm size is 1 when processId == lastId and processNmb - 1 otherwise

    // http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/
    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //                void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //                MPI_Comm comm)
    // root - rank of sending process (integer) 
    
    if (processId != lastId) {
      *adMatrixPartial = (int*) malloc(middleSize * sizeof(int));
      *nodesProcesNmb = middleNodes;
      MPI_Scatter(adMatrixFull, middleSize, MPI_INT, 
                 *adMatrixPartial, middleSize, MPI_INT, 
                 0, lastComm);
    } else {
      *adMatrixPartial = (int*) malloc(lastSize * sizeof(int));
      *nodesProcesNmb = lastNodes;
      MPI_Scatter(adMatrixFull + (lastId * middleSize), lastSize, MPI_INT, 
                 *adMatrixPartial, lastSize, MPI_INT, 
                 0, lastComm);
    }

  } else {
    // only one process - copy everything
    *adMatrixPartial = (int*) malloc(matrixSize * sizeof(int));
    memcpy(*adMatrixPartial, adMatrixFull, matrixSize * sizeof(int));
  }

}

// ====================== MAIN ========================
int main( int argc, char *argv[] )
{
  FILE  *file;
  char  *filename = argv[1];
  int   *adMatrixFull = 0;
  int   *adMatrixPartial = 0;
  int   nodesNmb = 0, nodesNmbProcess = 0;
  int   processId, processNmb;

  // MPI Initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processNmb);
  MPI_Comm_rank(MPI_COMM_WORLD, &processId);

  // Read edge list from file
	if (processId == 0) {
    file = fopen(filename, "r");
    if (file == 0) {
      printf("Cannot found input file %s!\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
		fscanfEdgeList(file, &adMatrixFull, &nodesNmb);
    fclose(file);
    fprintfAdMatrix(stdout, adMatrixFull, nodesNmb, nodesNmb);    
	}

  // Partitioning of adjacency matrix
  primPartitionMatrix(adMatrixFull, &adMatrixPartial, nodesNmb, &nodesNmbProcess, processId, processNmb);

  printf("%d, %d, %d\n", processId, nodesNmb, nodesNmbProcess);

  // Free allocated resources
  if(processId == 0) {
    free(adMatrixFull);
  }
  if (adMatrixPartial != 0) {
    free(adMatrixPartial);
  }


  MPI_Finalize();
  return 0;
}

