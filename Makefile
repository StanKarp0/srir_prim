exec = prim01
src = src/mpi/prim01.c

build: $(src)
	mpicc -o $(exec) $(src) -lm

run6:
	mpiexec -n 1 $(exec) graphs/nodes6edges12.csv

run10:
	mpiexec -n 2 $(exec) graphs/nodes10edges20.csv

run20:
	mpiexec -n 3 $(exec) graphs/nodes20edges50.csv

run50:
	mpiexec -n 4 $(exec) graphs/nodes50edges100.csv

clean:
	rm $(exec)