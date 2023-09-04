#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

//#define DEBUG

#define NUMBER_OF_LETTERS 26
//Assumption: seq(i) will always have smaller len than seq(1)
#define SEQ_1_MAX_LEN 3000
#define SEQ_I_MAX_LEN 2000

typedef struct { 
    // score
    // offset
    // k
    // seq
    // seq_len
} Seq_Info;


int main()
{

    return EXIT_SUCCESS;
}