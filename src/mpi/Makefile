exec = prim01
src = prim01.c

build: $(src)
	mpicc -o $(exec) $(src) -lm

run:
	mpiexec -n 3 $(exec) graphs/nodes10egdes20wages10.csv

run2:
	mpiexec -n 4 $(exec) graphs/nodes20egdes50wages10.csv

run3:
	mpiexec -n 2 $(exec) graphs/nodes6edges12.csv

clean:
	rm $(exec)