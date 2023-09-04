#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define DEBUG

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
void MS(char* seq, int k);
void fill_score_mat();
void init(int argc, char **argv);
void read_score_table(int argc, char **argv);
void printGraph();
int load_score_table_from_text_file( const char* fileName);
void skip_white_space();
int read_input_seq();
char* handle_string_from_input();
void toUpperCase(char *str);
void find_score_offset_MS(int *seq_score, char* seq_temp);
int check_score(char* seq, int offset);
void Work();

int main(int argc, char **argv)
{

    init(argc, argv);
    Work();






    free(score_table);
    free(seq1);
    for(int i = 0 ; i < num_of_seqs ; i++)
    {
        free(seq_arr[i]);
    }
    free(seq_arr);
    return EXIT_SUCCESS;
}

void init(int argc, char **argv)
{
    read_score_table(argc, argv);
    read_input_seq();


}

void Work()
{
    int seq_score[3] = {0};

    //char temp[SEQ_I_MAX_LEN];

    for(int i = 0 ; i < num_of_seqs ; i++)
    {
        //temp = 
        find_score_offset_MS(seq_score, seq_arr[i]);

        //printf("Seq = %s -> ", seq_arr[i]);
        printf("Highest alignment score = %d, Offset = %d, K = %d \n", seq_score[0], seq_score[1], seq_score[2]);

    }

}


void find_score_offset_MS(int *seq_score, char* seq_temp)
{
    int seq1_len = strlen(seq1);
    int seq_temp_len = strlen(seq_temp);
    int temp_score = 0;

    //MS
    for(int k = 0 ; k <= seq_temp_len ; k++ )
    {
        printf("@@ 1\n");
        
        MS(seq_temp, k);
        //offset
        for(int offset = 0 ; offset <= seq1_len - seq_temp_len ; offset++)
        {
            printf("@@ 2\n");

            temp_score = check_score(seq_temp, offset);
            if(temp_score > seq_score[0])
            {
                seq_score[0] = temp_score; // alignment score
                seq_score[1] = offset; // offset
                seq_score[2] = k; // K
            }
        }
    }

}


void MS(char* seq, int k)
{
    printf("@@ 3\n");

    if(k == strlen(seq)){
		return;
	}else{
        seq[k] = (seq[k] - 'A' + 1) % 26 + 'A';
        }
}

int check_score(char* seq, int offset)
{
    
    int value = 0;
    int seq_len = strlen(seq);
	for (int i = 0 ; i < seq_len ; i++){
        printf("@@ 1\n");
		value += score_table[(seq1[offset+i] -'A')  * (NUMBER_OF_LETTERS*NUMBER_OF_LETTERS) + (seq[i]-'A')];
	}
    return value;
}











/*-----------------------------------------------------------------------
@brief  Function to make a string letters all upper case.

@param param1 -- char* str ~> the string to be upper cased
@return: void
-----------------------------------------------------------------------*/
void toUpperCase(char *str) {
    if (str == NULL)
        return;

    for (int i = 0; str[i]; i++) {
        str[i] = toupper((unsigned char)str[i]);
    }
}




/*-----------------------------------------------------------------------
@brief  This function reads the input from "stdin" and initializing
        the variables needed.
            The input template:
                1.  Seq1 (String): (The main series).
                    its the series that all the other serieses will be compared to.
                    Assumptions : Seq1 maximum lenth is 3000.
                2.  number_of_sequences (int) : this will tell how many Strings 
                    there is to read (Seq(2) , Seq(3) ... Seq(n)).
                3.  Seq(1..n) : all the serieses line by line.
                (i.e) : 
                        ABBDAB
                        3
                        ADC
                        BDE
                        AAB

@return : int (Success\Fail)
-----------------------------------------------------------------------*/
int read_input_seq()
{
    char* seq1 = handle_string_from_input();

    // Read the num_of_seqs
    if (fscanf(stdin, "%d", &num_of_seqs) != 1) {
        printf("Failed to read input_number.\n");
        return 1;
    }
#ifdef DEBUG
    int seq1_len = strlen(seq1);
    printf("~~> Seq1 = %s , len = %d\n", seq1, seq1_len);
    printf("~~> num_of_seq = %d \n", num_of_seqs);
#endif
    fgetc(stdin); // (read the "\n" that come after the number)
    seq_arr = (char**)malloc(num_of_seqs * sizeof(char*));

    // Read the seq_array
    for (int i = 0; i < num_of_seqs; i++) 
    {
        seq_arr[i] = handle_string_from_input();
    }
#ifdef DEBUG
    for (int i = 0; i < num_of_seqs; i++) {
        printf("~~> %s \n", seq_arr[i]);
    }
#endif
    return 0;
}



/*-----------------------------------------------------------------------
@brief  Function to Handle the string input:
        this function allocates memory for a 3000 long string and
        after that cut's down the unused memory.

@return: char*
-----------------------------------------------------------------------*/
char* handle_string_from_input()
{
    char* temp = (char*)malloc(SEQ_1_MAX_LEN * sizeof(char));
    if (temp == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    // Read SEQ
    if (fgets(temp, SEQ_1_MAX_LEN, stdin) == NULL) {
        printf("Failed to read input_string.\n");
        exit(1);
    }
    int temp_len = strlen(temp);
    if (temp_len > 0 && temp[temp_len - 1] == '\n') 
    {
        temp[temp_len - 1] = '\0';
    }
    char* str = strdup(temp);
    toUpperCase(str);

    return str;
}


/*-----------------------------------------------------------------------
@brief  Function to initialize the score table.
        the function checks if there is an argument passed to the program
        and read the values from there to the score table.
        -if there is no argument passed to it, the function will initialize
        the score table with "0" exept the diagonal line (will be "1").
        -if there is not enogh values in the txt file then the function
        will padd the score table with "0".

@param param1 -- int argc ~> how many arguments passed at execution
@param param2 -- char** argv ~> an array of strings of the arguments
@ return: void
-----------------------------------------------------------------------*/
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


/*-----------------------------------------------------------------------
@brief  Function to Print the score table.

@ return: void
-----------------------------------------------------------------------*/
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