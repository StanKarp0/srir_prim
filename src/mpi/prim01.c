

/* calculate pi by Leibniz series
 * pi = 4 * SUM (n) { (-1)^n / 2n +1 }
 * Distribute terms of Leibnitz series among processes
 * in round-robin fashion.
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>
int main( int argc, char *argv[] )
{
    int n, myid, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double mysum, pi, sum, inv;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);


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

