#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <time.h>
#include <math.h>



#define NUMBER_OF_LETTERS 26
//Assumption: seq(i) will always be smaller len than seq(1)
#define SEQ_1_MAX_LEN 3000
#define SEQ_I_MAX_LEN 2000

char* seq1;
char **seq_arr;

int lineno = 1;
int num_of_seqs;
int *score_table;

//Functions:
char* MS(char* seq, int k);
void fill_score_mat();
void init(int argc, char **argv);
void read_score_table(int argc, char **argv);
void printGraph();
int load_score_table_from_text_file( const char* fileName);
void skip_white_space();
int read_input_seq();




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


int main(int argc, char **argv)
{

    init(argc, argv);

    free(score_table);
    return EXIT_SUCCESS;
}

void init(int argc, char **argv)
{
    read_score_table(argc, argv);
    read_input_seq();


}





int read_input_seq()
{

    char* temp = (char*)malloc(SEQ_1_MAX_LEN * sizeof(char));
    if (temp == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    // Read SEQ1
    if (fgets(temp, SEQ_1_MAX_LEN, stdin) == NULL) {
        printf("Failed to read input_string.\n");
        return 1;
    }
    int temp_len = strlen(temp);
    if (temp_len > 0 && temp[temp_len - 1] == '\n') 
    {
        temp[temp_len - 1] = '\0';
    }
    char* seq1 = strdup(temp);
    int seq1_len = strlen(seq1);
    printf("~~> Seq1 = %s , len = %d\n", seq1, seq1_len);


    
/*
    // Read the num_of_seqs
    if (fscanf(stdin, "%d", &num_of_seqs) != 1) {
        printf("Failed to read input_number.\n");
        return 1;
    }

    // Read the input_array
    for (int i = 0; i < num_of_seqs; i++) {
        if (fgets(seq_arr[i], sizeof(SEQ_I_MAX_LEN, stdin) == NULL) {
            printf("Failed to read input_array[%d].\n", i);
            return 1;
        }
    }

    // Access the variables as needed
    printf("Input String: %s", input_string);
    printf("Input Number: %d\n", input_number);
    printf("Input Array:\n");
    for (int i = 0; i < input_number; i++) {
        printf("%s", input_array[i]);
    }
*/

    return 0;
}









































void read_score_table(int argc, char **argv) 
{
    score_table = (int*)calloc(NUMBER_OF_LETTERS * NUMBER_OF_LETTERS , sizeof(int));
    if(score_table == NULL) { perror("malloc"); exit(1); }
    if(argc == 2)
    {   
        // Read Score table from a text file.
        char* fileName = argv[1];
        load_score_table_from_text_file( fileName);
        

    }else{
        // Default option (diagonal = 1, the rest = 0).
        for(int i = 0 ; i < NUMBER_OF_LETTERS ; i++)
            score_table[NUMBER_OF_LETTERS*i+i] = 1;
    }
    printGraph();
}


int load_score_table_from_text_file( const char* fileName)
{
    FILE* fp;
    int count_v = 0;
    int value;
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        perror("Error opening file");
        exit(1);
    }
    for(int i = 0 ; i < NUMBER_OF_LETTERS * NUMBER_OF_LETTERS ; i++)
    {
        if(fscanf(fp , "%d", &value) == 1)
        {
            count_v++;
            score_table[i] = value;
        }else{
            printf("no more values\n");
            break;
        }
    }
    if (count_v != NUMBER_OF_LETTERS*NUMBER_OF_LETTERS) 
    {
        fprintf(stderr, "%d values appear in the text file (expected %d values, (padding with \"0\"))\n", 
        count_v, NUMBER_OF_LETTERS*NUMBER_OF_LETTERS);
    }
    fclose(fp);
    return 0;
}


// can be used for debugging
void printGraph() 
{
    for (int i = 0; i < NUMBER_OF_LETTERS; i++) 
    {
        for (int j = 0; j < NUMBER_OF_LETTERS; j++)
        {
            printf("%4d", score_table[NUMBER_OF_LETTERS*i+j]);
        }
        putchar('\n');
    }
}

    
void skip_white_space() {
   int c;
   while(1) {
       if ((c = getchar()) == '\n')
           lineno++;
       else if (isspace(c))
           continue;
       else if (c == EOF)
           break;
       else {
         ungetc(c, stdin);
         break;
       }
   }
}