#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <mpi.h>
#include <time.h>
#include <math.h>



#define NUMBER_OF_LETTERS 26


int *score_table;

char* MS(char* seq, int k);
void fill_score_mat();
void init(int argc, char **argv);
void read_score_table(int argc, char **argv);
void printGraph();


int main(int argc, char **argv)
{

    
    printf("HI\n");
    init(argc, argv);



    return EXIT_SUCCESS;
}

/*-----------------------------------------------------------------------
    The input template:
        1.  Seq1 (String): (The main series).
            its the series that all the other serieses will be compared to.
            Assumptions : Seq1 maximum lenth is 3000.
        2.  number_of_sequences (int) : this will tell how many Strings 
            there is to read (Seq(2) , Seq(3) ... Seq(n)).
        3.  Seq(1..n) : all the serieses line by line.
        (i.e) : 
                ABBDAB
                4
                ADC
                BDE
                AAB
                CBA
-----------------------------------------------------------------------*/
void init(int argc, char **argv)
{
    read_score_table(argc, argv);
    


}


void read_score_table(int argc, char **argv) {
    /*
    int c;
    int w;
    int count_w = 0; // number of values read in so far
    */

    score_table = (int*)calloc(NUMBER_OF_LETTERS * NUMBER_OF_LETTERS , sizeof(int));
    if(score_table == NULL) { perror("malloc"); exit(1); }

    if(argc == 2)
    {
        // TODO: fill the score_table with values from the file in argv[1]
    }else{
        // Default option (diagonal = 1, the rest = 0) :
        int *next_entry = score_table;
        for(int i = 0 ; i < NUMBER_OF_LETTERS ; i++)
            next_entry[NUMBER_OF_LETTERS*i+i] = 1;
    }
    printGraph();
    
}


// can be used for debugging
void printGraph() {

    printf("Score Table:\n");
    for (int i = 0; i < NUMBER_OF_LETTERS; i++) 
    {

        for (int j = 0; j < NUMBER_OF_LETTERS; j++)
        {
            if(i == 0)
            {
                printf("%c  ", 'A'+j);
            }else if(j != 0)
            {
                printf("%d  ", score_table[NUMBER_OF_LETTERS*i+j]);
            }else{
                printf("%c  ", 'A'+i);
            }
        }
        putchar('\n');
    }
}