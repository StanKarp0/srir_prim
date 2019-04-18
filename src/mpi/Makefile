exec = prim01
src = prim01.c

build: $(src)
	mpicc -o $(exec) $(src) -lm

run:
	mpiexec -n 3 $(exec) graphs/nodes10egdes20wages10.csv

clean:
	rm $(exec)