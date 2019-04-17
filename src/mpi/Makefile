exec = prim01
src = prim01.c

build: $(src)
	mpicc -o $(exec) $(src) -lm

run:
	mpiexec -n 3 $(exec)

clean:
	rm $(exec)