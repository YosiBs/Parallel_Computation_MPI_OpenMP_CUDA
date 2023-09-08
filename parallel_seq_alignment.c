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
#define TABLE_SIZE 26*26
//Assumption: seq(i) will always have smaller len than seq(1)
#define SEQ_1_MAX_LEN 3000
#define SEQ_2_MAX_LEN 2000
#define MASTER 0
#define INITSCORE -1000
enum tags { STOP, WORK, INIT_PACK} tag;

typedef struct { 
    int score;
    int offset;
    int k;
    int seq_len;
    char seq[SEQ_2_MAX_LEN];
} Seq_Info;

typedef struct {
    int num_of_seqs;
    int seq1_len;
    int score_table[NUMBER_OF_LETTERS][NUMBER_OF_LETTERS];
    char seq1[SEQ_1_MAX_LEN];
} init_variables;


char* seq1;

init_variables init_pack;
Seq_Info* seq_arr;
int lineno = 1;
int start_end[2];
int chunk_size;



// Functions:
void MS(char* seq, int k);
void fill_score_mat();
void read_score_table(int argc, char **argv);
void print_score_table();
int load_score_table_from_text_file( const char* fileName);
void skip_white_space();
int read_input_seq();
char* handle_string_from_input();
void toUpperCase(char *str);
void find_score_offset_MS(int *seq_score, char* seq_temp);
int check_score(char* seq, int offset);
//void Work();
void init(int argc, char **argv,int myRank, int nprocs);
//void exit_safely();
void masterProcess(int nprocs);
void workerProcess(int nprocs, int myRank, int reminder);
void get_MPI_init_Datatype(MPI_Datatype* custom_type);
void get_MPI_Seq_Info_Datatype(MPI_Datatype* custom_type);


int main(int argc, char **argv)
{

    int myRank;
	int nprocs; 	
    int reminder = 0;
    
    /* Start up MPI */
	MPI_Init(&argc, &argv);
 	/* Find number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	/* Find process myRank */
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    

    init(argc, argv, myRank, nprocs);


    if(myRank == MASTER)
        masterProcess(nprocs);
    else
        workerProcess(nprocs, myRank, reminder);




    //exit_safely();
    MPI_Finalize();
    return EXIT_SUCCESS;
}

void masterProcess(int nprocs)
{
        double start_time = MPI_Wtime();
        //Work(nprocs, MASTER);
    
        printf("Runtime: %lf\n", MPI_Wtime() - start_time);
}


void workerProcess(int nprocs, int myRank, int reminder)
{
    //Work();
}



void init(int argc, char **argv,int myRank, int nprocs)
{

    MPI_Datatype mpi_init_datatype;
    MPI_Datatype mpi_seq_info_datatype;
    get_MPI_init_Datatype(&mpi_init_datatype);
    get_MPI_Seq_Info_Datatype(&mpi_seq_info_datatype);

    if(myRank == MASTER)
    {
        read_score_table(argc, argv);
        read_input_seq();
        
    }
    MPI_Bcast(&init_pack, 1, mpi_init_datatype, MASTER, MPI_COMM_WORLD);

#ifdef DEBUG
    //printf("~~>[%d] num of seqs = %d, seq1 len= %d seq1 = %s\n", myRank, init_pack.num_of_seqs, init_pack.seq1_len, init_pack.seq1);
    //print_score_table();
#endif
    
    chunk_size = 0;
    int reminder = init_pack.num_of_seqs%(nprocs-1);
    if(myRank == MASTER)
    {
        chunk_size = init_pack.num_of_seqs;
        start_end[0] = 0;
        start_end[1] = init_pack.num_of_seqs - 1;

        int elements_per_proc[nprocs-1];
        int displacements[nprocs-1];

        for(int i = 0 ; i < nprocs - 1 ;i++ )
        {
            elements_per_proc[i] = (init_pack.num_of_seqs/(nprocs-1)) + ((myRank == nprocs - 1) ? reminder : 0);
            displacements[i] = (i * (init_pack.num_of_seqs/(nprocs-1)) ) * sizeof(Seq_Info);
#ifdef DEBUG
            printf("~~%d. num_elements=%d, displacement=%d\n", i+1, elements_per_proc[i],displacements[i]);
#endif
        }
        MPI_Scatterv(seq_arr, elements_per_proc, displacements, mpi_seq_info_datatype, MPI_IN_PLACE, 0, mpi_seq_info_datatype, MASTER, MPI_COMM_WORLD);
#ifdef DEBUG
        printf("~~%d. AFTER SCATTER\n", myRank);
#endif











    }else{
        chunk_size = (init_pack.num_of_seqs/(nprocs-1)) + ((myRank == nprocs - 1) ? reminder : 0);
        start_end[0] = (myRank - 1) * (init_pack.num_of_seqs / (nprocs - 1)) + 1;
        start_end[1] = start_end[0] + chunk_size - 1; 

        seq_arr = (Seq_Info*)malloc(chunk_size * sizeof(Seq_Info));
        if (seq_arr == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
        }
        printf("~~%d. BEFORE SCATTER\n", myRank);
        MPI_Scatterv(NULL, NULL, NULL, mpi_seq_info_datatype, seq_arr, chunk_size * sizeof(Seq_Info), mpi_seq_info_datatype, MASTER, MPI_COMM_WORLD);


        printf("~~%d. AFTER SCATTER\n", myRank);


#ifdef DEBUG
        for(int i = 0 ; i < chunk_size ; i++)
            printf("~~>[%d] seq[%d] = %s\n", myRank, start_end[i] , seq_arr[i].seq);
#endif
    }





}





























void get_MPI_init_Datatype(MPI_Datatype* custom_type)
{
    int block_lengths[4] = {1, 1, TABLE_SIZE, SEQ_1_MAX_LEN};
    MPI_Aint displacements[4];
    MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_CHAR};

    displacements[0] = 0;
    displacements[1] = sizeof(int); 
    displacements[2] = 2 * sizeof(int);
    displacements[3] = (2 + TABLE_SIZE) * sizeof(int);
    MPI_Type_create_struct(4, block_lengths, displacements, types, custom_type);
    MPI_Type_commit(custom_type);
}


void get_MPI_Seq_Info_Datatype(MPI_Datatype* custom_type)
{
    int blocklengths[5] = {1, 1, 1, 1, SEQ_2_MAX_LEN};
    MPI_Aint displacements[5];
    MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR};

    for(int i = 0 ; i < 5 ; i++)
        displacements[i] = i * sizeof(int);

    MPI_Type_create_struct(5, blocklengths, displacements, types, custom_type);
    MPI_Type_commit(custom_type);
}





























































/*-----------------------------------------------------------------------
@brief ~~> void Work()
    The Function will initialize an int array with size 3 to "0".
    and will run the function "find_score_offset_MS" for each seq.
    

@return: void
-----------------------------------------------------------------------*/
/*
void Work()
{
    int seq_score[3] = {0}; // Score , offset and K, will be written to these array
    char temp[SEQ_2_MAX_LEN];
    for(int i = 0 ; i < num_of_seqs ; i++)
    {
        strcpy(temp, seq_arr[i]);
        find_score_offset_MS(seq_score, seq_arr[i]);


        printf("%d. Seq = %s -> ", i+1, temp);
        printf("Highest alignment score = %d, Offset = %d, K = %d \n", seq_score[0], seq_score[1], seq_score[2]);
    }
}
*/


/*-----------------------------------------------------------------------
@brief ~~> void find_score_offset_MS(int *seq_score, char* seq_temp)
    The Function will find the best score, offset, k.
    the answer will be applied to seq_score array.
        seq_score[0] = Best Score value.
        seq_score[1] = Best Offset value.
        seq_score[2] = Best k value.

@param param1 -- int *seq_score: an array of 3 elements. 
@param param1 -- char* seq_temp : the current seq to work on
@return: void
-----------------------------------------------------------------------*/
void find_score_offset_MS(int *seq_score, char* seq_temp)
{
    int seq_temp_len = strlen(seq_temp);
    int temp_score = 0;

    for(int k = seq_temp_len ; k >= 0 ; k-- )
    {
        MS(seq_temp, k);
        
        for(int offset = 0 ; offset <= init_pack.seq1_len - seq_temp_len ; offset++)
        {
            temp_score = check_score(seq_temp, offset);
#ifdef DEBUG
            printf("~~> score=%d \n", temp_score);
#endif
            if(temp_score > seq_score[0])
            {
                seq_score[0] = temp_score; // alignment score
                seq_score[1] = offset; // offset
                seq_score[2] = k; // K
            }
        }
    }
}
/*-----------------------------------------------------------------------
@brief ~~> void MS(char* seq, int k)
    The Function will change the current seq with a given K.

@param param1 -- char* seq : the current seq to work on. 
@param param1 -- int k : the current k.
@return: void
-----------------------------------------------------------------------*/
void MS(char* seq, int k)
{
    if(k == strlen(seq)){
		return;
	}else{
        seq[k] = (seq[k] - 'A' + 1) % 26 + 'A';
        }
}

/*-----------------------------------------------------------------------
@brief ~~> int check_score(char* seq, int offset)
    The Function will calculate the score for a specific K and Offset.

@param param1 -- char* seq : the current seq to work on. 
@param param1 -- int offset : the current offset. 
@return: int value : the best score for the current seq.
-----------------------------------------------------------------------*/
int check_score(char* seq, int offset)
{
    int value = 0;
    int seq_len = strlen(seq);
	for (int i = 0 ; i < seq_len ; i++)
    {
		value += init_pack.score_table[(init_pack.seq1[offset+i] - 'A')][(seq[i] - 'A')];
#ifdef DEBUG
        //printf("~~> table slot: %d * %d + %d\n",(seq1[offset+i] - 'A'),NUMBER_OF_LETTERS, (seq[i] - 'A'));
        printf("~~> value so far = %d\n",value);
#endif
    }
    return value;
}

/*-----------------------------------------------------------------------
@brief  ~~> void exit_safely()
    Function to free all the memory allocated in the program.

@return: void
-----------------------------------------------------------------------*/
/*
void exit_safely()
{
    free(score_table);
    free(seq1);
    for(int i = 0 ; i < num_of_seqs ; i++)
        free(seq_arr[i]);
    free(seq_arr);
}
*/



/*-----------------------------------------------------------------------
@brief  ~~> void toUpperCase(char *str)
    Function to make a string letters all upper case.

@param param1 -- char* str ~> the string to be upper cased
@return: void
-----------------------------------------------------------------------*/
void toUpperCase(char *str) 
{
    if (str == NULL)
        return;

    for (int i = 0; str[i]; i++) {
        str[i] = toupper((unsigned char)str[i]);
    }
}




/*-----------------------------------------------------------------------
@brief  ~~> int read_input_seq()
    This function reads the input from "stdin" and initializing
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
    seq1 = handle_string_from_input();
    init_pack.seq1_len = strlen(seq1);
    strcpy(init_pack.seq1, seq1);
    init_pack.seq1[init_pack.seq1_len] = '\0';
    
    // Read the num_of_seqs value
    if (fscanf(stdin, "%d", &init_pack.num_of_seqs) != 1) {
        printf("Failed to read input_number.\n");
        return 1;
    }
#ifdef DEBUG
    printf("~~> Seq1 = %s , len = %d\n", init_pack.seq1, init_pack.seq1_len);
    printf("~~> num_of_seq = %d \n", init_pack.num_of_seqs);
#endif
    fgetc(stdin); // (read the "\n" that come after the number)
    seq_arr = (Seq_Info*)malloc(init_pack.num_of_seqs * sizeof(Seq_Info));
    if (seq_arr == NULL) {
    fprintf(stderr, "Memory allocation failed.\n");
    exit(1);
    }
    // Read the seq_array
    for (int i = 0; i < init_pack.num_of_seqs; i++) 
    {
        char* temp = handle_string_from_input();
        seq_arr[i].seq_len = strlen(seq1);
        strcpy(seq_arr[i].seq, temp);
        seq_arr[i].seq[init_pack.seq1_len] = '\0';
        seq_arr[i].seq_len = strlen(seq_arr[i].seq);
        seq_arr[i].offset = 0;
        seq_arr[i].k = 0;
        seq_arr[i].score = INITSCORE;
    }

#ifdef DEBUG
    //Print all the seqs
    for (int i = 0; i < init_pack.num_of_seqs; i++) {
        printf("~~> %s \n", seq_arr[i].seq);
    }
#endif
    return 0;
}


/*-----------------------------------------------------------------------
@brief  ~~> char* handle_string_from_input()
    Function to Handle the string input:
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
@brief  ~~>void read_score_table(int argc, char **argv) 
    Function to initialize the score table.
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
    
    if(argc == 2)
    {   
        // Read Score table from a text file.
        char* fileName = argv[1];
        load_score_table_from_text_file( fileName);
        

    }else{
        // Default option (diagonal = 1, the rest = 0).
        for(int i = 0 ; i < NUMBER_OF_LETTERS ; i++)
            init_pack.score_table[i][i] = 1;
    }
    print_score_table();
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
    for(int i = 0 ; i < NUMBER_OF_LETTERS ; i++)
    {
        for(int j = 0; j < NUMBER_OF_LETTERS ; j++ )
        {

            if(fscanf(fp , "%d", &value) == 1)
            {
                count_v++;
                init_pack.score_table[i][j] = value;
            }else{
                printf("no more values\n");
                break;
            }
        }
    }

#ifdef DEBUG
    if (count_v != NUMBER_OF_LETTERS*NUMBER_OF_LETTERS) 
    {
        printf("%d values appear in the text file (expected %d values, (padding with \"0\"))\n", 
        count_v, NUMBER_OF_LETTERS*NUMBER_OF_LETTERS);
    }
#endif

    fclose(fp);
    return 0;
}


/*-----------------------------------------------------------------------
@brief  ~~> void print_score_table() 
    Function to Print the score table.

@ return: void
-----------------------------------------------------------------------*/
void print_score_table() 
{
    for (int i = 0; i < NUMBER_OF_LETTERS; i++) 
    {
        for (int j = 0; j < NUMBER_OF_LETTERS; j++)
        {
            printf("%4d", init_pack.score_table[i][j]);
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