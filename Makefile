CC = mpicc -fopenmp
CFLAGS = -Wall -Werror $(DEBUG)

all: serial_sequence_alignment

serial_sequence_alignment: serial_sequence_alignment.c
	mpicc $(CFLAGS) -o s_sa serial_sequence_alignment.c -lm

# run the serial program:
run_s:
	mpiexec -n 1 ./s_sa st1.txt < input.txt

.PHONY: clean

clean:
	rm -f s_sa serial_sequence_alignment.o 
