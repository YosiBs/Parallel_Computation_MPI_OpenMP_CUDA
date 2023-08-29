CC = mpicc -fopenmp
CFLAGS = -Wall -Werror $(DEBUG)

all: serial_sequence_alignment

serial_sequence_alignment: serial_sequence_alignment.c
	mpicc $(CFLAGS) -o serial_sequence_alignment serial_sequence_alignment.c -lm

genGraph: genGraph.c
	gcc $(CFLAGS) -o genGraph genGraph.c -lm

run_serial:
	mpiexec -n 1 ./serial_sequence_alignment

.PHONY: clean

clean:
	rm -f serial_sequence_alignment serial_sequence_alignment.o 
