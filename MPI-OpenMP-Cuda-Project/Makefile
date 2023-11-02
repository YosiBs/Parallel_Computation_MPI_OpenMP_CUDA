CC = mpicc -fopenmp
CFLAGS = -Wall -Werror $(DEBUG)

all: serial_seq_alignment build_parallel

serial_seq_alignment: serial_seq_alignment.c
	mpicc $(CFLAGS) -o s_sa serial_seq_alignment.c -lm

build_parallel:
	# note: cFunctions.c is a "C" file (not C++) but it can be compiled as
	# a C++ file with mpicxx (which uses a C++ compiler). 
	# gcc -c cFunctions.c -o cFunctions.o 
	nvcc -gencode arch=compute_61,code=sm_61 -c cuda_part.cu -o cuda_part.o
	mpicxx -fopenmp -c parallel_seq_alignment.c -o parallel_seq_alignment.o
	#linking:
	mpicxx  -fopenmp -o p_sa cuda_part.o parallel_seq_alignment.o  -L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcudart
	



# run the serial program:
run_s:
	mpiexec -n 1 ./s_sa inputs/score-table-2.txt < inputs/input-1.txt > outputs/answer-s.txt

# run the parallel program:
run_p:
	mpiexec -n 4 ./p_sa inputs/score-table-2.txt < inputs/input-1.txt > outputs/answer-p.txt

.PHONY: clean

clean:
	rm -f s_sa serial_seq_alignment.o p_sa parallel_seq_alignment.o cuda_part *.o *.btr
