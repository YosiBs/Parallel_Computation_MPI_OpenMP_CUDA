CC = mpicc -fopenmp
CFLAGS = -Wall -Werror $(DEBUG)

all: serial_seq_alignment parallel_seq_alignment

serial_seq_alignment: serial_seq_alignment.c
	mpicc $(CFLAGS) -o s_sa serial_seq_alignment.c -lm

parallel_seq_alignment: parallel_seq_alignment.c
	mpicc $(CFLAGS) -o p_sa parallel_seq_alignment.c -lm

# run the serial program:
run_s:
	mpiexec -n 1 ./s_sa score-table-2.txt < input-2.txt

run_p:
	mpiexec -n 4 ./p_sa score-table-2.txt < input-2.txt

.PHONY: clean

clean:
	rm -f s_sa serial_seq_alignment.o p_sa parallel_seq_alignment.o
