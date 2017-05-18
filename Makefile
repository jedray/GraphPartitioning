DATA=data
INPUT=karate.graph
SRC=src/partition.c src/bisection.c src/parser.c src/matrix.c src/vector.c 

all: compile run


compile:
	mpiicc -Wall  $(SRC) -I include -o ./bin/do_par

run:
	srun --nodes=1 -A edu17.SF2568 -t 1 mpirun -np 20  ./bin/do_par $(DATA)/$(INPUT) 
