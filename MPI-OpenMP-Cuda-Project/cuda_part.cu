#include <stdio.h>
#include <string.h>

#define BLOCK_DIM 1024 // number of threads in a block
#define TABLE_SIZE 26
//#define DEBUG

typedef struct
{
    int score;
    int offset;
    int k;
    int seq_len;
    char seq[2000];
} Seq_Info;





__global__ void work_cuda(char* d_seq2, int seq2_len,char* d_seq1,int seq1_len,int score_table[][26],int table_size,int* arr_max_result);

int computeOnGPU(char* seq1, char* seq2 ,int seq1_length,int seq2_length, int score_table[][26],int table_size,int* max_arr);
__device__ void scan_plus(int *array, int size);
__device__ void scan_maxl(int *array, int size);
__device__ char next_char(char c);


__global__ void work_cuda(char* d_seq2, int seq2_len,char* d_seq1,int seq1_len,int score_table[][26],int table_size,int* arr_max_result)
{

    int tid = threadIdx.x;
    __shared__ int scores[2* BLOCK_DIM];

        if(threadIdx.x == 0 ){
            scores[0] = 0;
            for(int i = 0 ; i<seq2_len ; i++){
                scores[0] += score_table[d_seq1[i+blockIdx.x]-'A'][next_char(d_seq2[i])-'A'];
                 }
            }
      

            while(tid < 2*BLOCK_DIM - 1){
                if(tid< seq2_len){
                int Score = score_table[d_seq1[tid+blockIdx.x]-'A'][d_seq2[tid]-'A'];
                int ScoreM = score_table[d_seq1[tid+blockIdx.x]-'A'][next_char(d_seq2[tid])-'A'];
                scores[tid+1] = Score-ScoreM;
                }else{
                     scores[tid+1] = 0;
                }
                tid += BLOCK_DIM;
                
            }


    __syncthreads();
    scan_plus(scores, BLOCK_DIM*2);
    __syncthreads();
    scan_maxl(scores, BLOCK_DIM*2);

    __syncthreads();

    tid = threadIdx.x ;

    while(tid < 2*BLOCK_DIM ){
        if (tid == 0 && scores[0] == scores[2*BLOCK_DIM - 1]){
            arr_max_result[2*blockIdx.x] = scores[tid]; //score
            arr_max_result[2*blockIdx.x+1] = tid;//k
        }else if(tid != 0 && scores[tid-1] != scores[tid] && scores[tid] == scores[2*BLOCK_DIM - 1]){
            arr_max_result[2*blockIdx.x] = scores[tid]; //score
            arr_max_result[2*blockIdx.x+1] = tid;//k
        }
    tid += BLOCK_DIM;
    }


}



__device__ void scan_maxl(int *array, int size)
{
   for (unsigned int stride=1; stride <= size/2; stride *= 2) {
        int v,v1;
        int tid =threadIdx.x + blockDim.x ;

        if (threadIdx.x >= stride) {
            v = array[threadIdx.x - stride];
        }
        if (tid >= stride) {
            v1 = array[ tid- stride];
        }
        __syncthreads(); // wait untill all threads get to this line

        if (threadIdx.x >= stride && array[threadIdx.x] < v )
            array[threadIdx.x] = v;

        if (tid >= stride && array[tid] < v1 )
            array[tid] = v1;

        __syncthreads(); // wait untill all threads get to this line
     }
     
}


__device__ void scan_plus(int *array, int size)
{
   for (unsigned int stride=1; stride <= size/2; stride *= 2) {
        int v,v1;
        int tid =threadIdx.x + blockDim.x;

        if (threadIdx.x >= stride) {
            v = array[threadIdx.x - stride];
        }

        if (tid >= stride) {
            v1 = array[ tid- stride];
        }
        
        __syncthreads(); // wait untill all threads get to this line

        if (threadIdx.x >= stride)
            array[threadIdx.x] += v;

        if (tid >= stride)
            array[tid] += v1;   

        __syncthreads(); // wait untill all threads get to this line
     }
     
} 

int computeOnGPU(char* seq1, char* seq2 ,int seq1_length,int seq2_length, int score_table[][26],int table_size,int* max_arr)
{
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
   
    // Allocate memory on GPU to copy the data from the host
  char* d_seq1;
    err = cudaMalloc((void **)&d_seq1, seq1_length*sizeof(char));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        return 1;
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_seq1, seq1, seq1_length*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        return(1);
    }
    
    char *d_seq2;
    err = cudaMalloc((void **)&d_seq2, (seq2_length)*sizeof(char));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(-1);
    }
      err = cudaMemcpy(d_seq2, seq2, (seq2_length)*sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "(4)Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(-1);
    }
    
    // Allocate memory on GPU to copy the TABLE from the host
    int *d_table;
    err = cudaMalloc((void **)&d_table, table_size*table_size*sizeof(int));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        return 1;
    }

    // Copy TABLE from host to the GPU memory
    err = cudaMemcpy(d_table, score_table, table_size*table_size*sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        return(1);
    }

    dim3 blockDim;
    blockDim.x = BLOCK_DIM;
    int blocks = seq1_length - seq2_length + 1;

    int* arr_max_result;
          err = cudaMalloc((void **)&arr_max_result, 2*blocks*sizeof(int));
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(-1);
    }
   



    work_cuda<<<blocks, blockDim.x>>>(d_seq2,seq2_length,d_seq1, seq1_length,(int(*)[TABLE_SIZE])d_table,table_size,arr_max_result);
    /* note: next lines may be executed before the kernel is done */
    err = cudaGetLastError();
    if (err != cudaSuccess) {
         fprintf(stderr, "Failed to launch incrementByOne kernel -  %s\n", cudaGetErrorString(err));
        return(1);
        }


  
    // Copy the  result from GPU to the host memory.
    err = cudaMemcpy(max_arr, arr_max_result, 2*blocks*sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        return(1);
    }

#ifdef DEBUG
    
    for(int i = 0 ; i < blocks; i++){
        printf("offset  = %d, k = %d  max score:%d\n",i,max_arr[(i*2)+1],max_arr[i*2] );
    }

#endif


    // Free allocated memory on GPU
    if (cudaFree(arr_max_result) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        return(1);
    }
    if (cudaFree(d_table) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        return(1);
    }

    return 0;



}

__device__ char next_char(char c)
{
    c = (c == 'Z') ? 'A' : c+1;
    return c;
}