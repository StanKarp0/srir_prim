/*
 * Prim's Algorithm
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


// ================== STRUCTURES ===================
typedef struct DWeight {
  // Structure is a part of d table - table of minimal weights to node
  int weight;
  int node; // global value
} DWeight;

typedef struct DEdge {
  // 
  int fromNode;
  int toNode;
  int weight;
} DEdge;

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

void fprintfDTable(FILE* file, DWeight* dTable, int rowNmb) {
  // Function prints d weights table to output stream

  for (int row = 0; row < rowNmb; row++) {
    fprintf(file, "%d|%d,%d\n", row, dTable[row].node, dTable[row].weight);
  }
}

void fprintfDEdges(FILE* file, DEdge* edes, int rowNmb) {
  // Function prints edges to output stream

  for (int row = 0; row < rowNmb; row++) {
    fprintf(file, "%d: %d -> %d - %d\n", row, edes[row].fromNode, edes[row].toNode, edes[row].weight);
  }
}

// ==================== ALGORITM ===================
void primPartitionMatrix(int *adMatrixFull, int nodesNmb, int processId, int processNmb, int **adMatrixPartial, int *nodesProcesNmb, int *startNode) {
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
  
  // define start node of each partition
  *startNode = middleNodes * processId;

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
      // TODO check if root should be 0 - this is different comm
      MPI_Scatter(adMatrixFull, middleSize, MPI_INT, 
                 *adMatrixPartial, middleSize, MPI_INT, 
                 0, lastComm);
    } else {
      *adMatrixPartial = (int*) malloc(lastSize * sizeof(int));
      *nodesProcesNmb = lastNodes;
    }

    // handle lastProcess. root not in lastComm so scatter can not be used
    // simple send and receive
    if (processId == 0) {
      MPI_Send(adMatrixFull + (lastId * middleSize), lastSize, MPI_INT, lastId, 0, MPI_COMM_WORLD);
    }
    else if (processId == lastId){
      MPI_Recv(*adMatrixPartial, lastSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  } else {
    // only one process - copy everything
    *adMatrixPartial = (int*) malloc(matrixSize * sizeof(int));
    memcpy(*adMatrixPartial, adMatrixFull, matrixSize * sizeof(int));
  }
}

void primPartitionDArray(int nodesNmbProcess, int nodesNmb, int firstNode, int* adMatrixPartial, DWeight** dTable) {
  // The partitioning of distance array among processes  
  // Function creates DTable for each process 

  *dTable = (DWeight*) malloc(nodesNmbProcess * sizeof(DWeight));

  for (int row = 0; row < nodesNmbProcess; row++) {
    (*dTable)[row].weight = adMatrixPartial[row * nodesNmb + firstNode];
    (*dTable)[row].node = firstNode;
  }
}

void primFindMinimum(int startNode, int nodesNmbProcess, DWeight* dTable, int* isAdded, DWeight *local) {
  // Function finds local minimum for partitioned adMatrix
  // Function returns -1 when connection not found 

  int globalNode  = 0;
  int localWeight = 0;
  local->weight   = INT_MAX;
  local->node     = -1;

  for (int localNode = 0; localNode < nodesNmbProcess; localNode++) {
    globalNode = localNode + startNode;

    if (!isAdded[globalNode]) {
      localWeight = dTable[localNode].weight;
      if (localWeight != 0 && (local->node == -1 || local->weight > localWeight)) {
        local->weight = localWeight;
        local->node   = globalNode;
      }
    }
  }
}

void primBroadcastSolution(int startNode, int nodesNmbProcess, DWeight *dTable, DWeight globalMin, DEdge* edge) {
  // Function brodcasts solution to every process. 

  int row = 0;
  int fromNode = -1, fromNodeGlobal;
  MPI_Bcast(&globalMin, 1, MPI_2INTEGER, 0, MPI_COMM_WORLD);

  // finding responisble partition
  if (startNode <= globalMin.node && globalMin.node < startNode + nodesNmbProcess) {
    for (; row < nodesNmbProcess && dTable[row].weight != globalMin.weight; row++);
    fromNode = dTable[row].node;
  } 

  // sync edges to other processes 
  MPI_Reduce(&fromNode, &fromNodeGlobal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&fromNodeGlobal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  edge->weight = globalMin.weight;
  edge->toNode = globalMin.node;
  edge->fromNode = fromNodeGlobal;
}

void primAlgorithm(int *adMatrix, int nodesNmb, int processId, int processNmb, DEdge* edges) {
  // Prim's Algorithm 

  int     *adMatrixPartial = 0;     // chunk of adMatrix 
  int     *isNodeAdded = 0;         // stores 1 if node is already in tree
  int     nodesNmbProcess = 0;      // nodes per process
  int     firstNode = 0;            // algorithm start node
  int     edgesNmb = nodesNmb - 1;  // edge count
  int     startNode = 0;            // start node number for partition

  DWeight *dTable = 0;              // weights array
  DWeight localMin;
  DWeight globalMin;

  // Partitioning of adjacency matrix and distance array
  primPartitionMatrix(adMatrix, nodesNmb, processId, processNmb, &adMatrixPartial, &nodesNmbProcess, &startNode);
  primPartitionDArray(nodesNmbProcess, nodesNmb, firstNode, adMatrixPartial, &dTable);

  // Algorithm initialization
  isNodeAdded = (int*) malloc(nodesNmb * sizeof(int));

  // Main iteration
  for (int index = 0; index < 1/*edgesNmb*/; index++) {

    DEdge *edge = &edges[index];

    // 1. Each process P_i computes d_i = min{dTable}
    primFindMinimum(startNode, nodesNmbProcess, dTable, isNodeAdded, &localMin);
    
    // 2. Global minimum d is then obtained by using all-to-one reduction operation
    //    and its stored in P_0 process. P_0 stores new vortex u
    MPI_Reduce(&localMin, &globalMin, 1, MPI_2INTEGER, MPI_MINLOC, 0, MPI_COMM_WORLD);
    
    // 3. Process P_0 broadcasts u one-to-all. The process responsible for u 
    //    marks u as belonging to tree.
    primBroadcastSolution(startNode, nodesNmbProcess, dTable, globalMin, edge);
    
    // 4. Each process updates the values od d[v] for its local vertices 

  }
  fprintfDEdges(stdout, edges, 3);


  free(isNodeAdded);
  free(adMatrixPartial);
  free(dTable);
}

// ===================== MAIN ======================
int main( int argc, char *argv[] )
{
  FILE  *file;
  char  *filename = argv[1];
  int   processId, processNmb;
  int   nodesNmb = 0;
  int   *adMatrix = 0;
  DEdge *edges;

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
		fscanfEdgeList(file, &adMatrix, &nodesNmb);
    fclose(file);
    fprintfAdMatrix(stdout, adMatrix, nodesNmb, nodesNmb);    
	}

  MPI_Bcast(&nodesNmb, 1, MPI_INT, 0, MPI_COMM_WORLD);
  edges = (DEdge*) malloc((nodesNmb - 1) * sizeof(DEdge));
  primAlgorithm(adMatrix, nodesNmb, processId, processNmb, edges);

  // Free allocated resources
  if(processId == 0) {
    free(adMatrix);
  }
  free(edges);

  MPI_Finalize();
  return 0;
}

