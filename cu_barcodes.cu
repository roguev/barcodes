/* cu_barcodes
 * computes a set of DNA sequences where each sequence in the set
 * differs from every other sequence by a minimum of given number
 * of nucleotides. THis is an interesting problem as it becomes more
 * difficult the larger the set becomes. This is the CUDA C implementation
 * Written by Assen Roguev, 2018
 */

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <ctime>

// compile with -DCUDA_ERROR_CHECK toturn on error checking
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )

inline void __cudaCheckError( const char *file, const int line ) {
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", cudaGetErrorString( err ) );
        exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}

inline void __cudaSafeCall( cudaError err, const char *file, const int line ) {
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err ) {
        fprintf( stderr, "%s:%i : %s\n", file, line, cudaGetErrorString( err ) );
        
        exit( -1 );
    }
#endif
    return;
}

/* gen_rand_sequence generates a random DNA sequence
 * of given length
 * Arguments:
 *  char* to hold the generated sequence
 *  int sequence length
 */ 
void gen_rand_sequence( char* b, int size ) {
    const char charset[] = "ATGC";
    
    // attempting to generate very random sequences
    // so seeding the random generator every time the
    // function is called
    clock_t t = clock();
    srand(t);
    
    if (size) {
        for (int n = 0; n < size; n++) {
            int key = rand() % (int) (sizeof charset - 1);
            b[n] = charset[key];
        }
    }
}

/* check_sequence is a CUDA kernel that checks a
 * sequence against a pool of sequences
 * Arguments:
 *  device char* to a sequence to be chacked
 *  device char* to a pool of sequences to check against
 *  device int* to number of pool sequences passed
 *  int to number of sequences in the pool
 *  int sequence length
 *  int minimum desired mismatches
 */
__global__
void check_sequence( char* d_seq, char* d_pool, int* d_pass, int pool_size, int LEN, int diffs ) {
    
    // to save some atomicAdd time on the slower global memory
    __shared__ int sh_pass;
    
    // number of differences
    int n_diff = 0;
    
    // position within the valid sequence pool
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    
    // set the share variable to 0
    if (threadIdx.x == 0)
        sh_pass = 0;
    __syncthreads();
    
    // overflow control
    if (ix < pool_size) {
        // check strings letter by letter and indicate differences
        for (int j = 0; j < LEN; j++) {
            if (d_seq[j] != d_pool[ix*LEN + j]) {
                if (++n_diff == diffs) {
                    atomicAdd(&sh_pass, 1); 
                    break;
                }
            }
        }
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        atomicAdd(d_pass, sh_pass);
}

/* populate_h_seq generates random sequences and puts
 * them in a single string
 * Argumenst:
 *  char* to a target string
 *  int desired sequence lenght
 *  int desired number of sequences
 */
void populate_h_seq( char* h_seq, int LEN, int N ) {
    // delete contents of the target string
    memset(h_seq, 0, N*LEN*sizeof(char));
    
    // generate random sequenes
    for (int i = 0; i < N; i++)
        gen_rand_sequence(h_seq+i*LEN, LEN);
}
    
/* main expects 3 command line argumenst
 *  int length of desired sequences
 *  int minimum number of mismatches
 *  int number of sequences to generate
 */
int main( int argc, char** argv ) {
    // change this according to available memory
    int CHUNK = 100000;
    
    if (argc < 3) {
        fprintf(stderr, "\n%s: Insufficient arguments\n\n", argv[0]);
        exit(1);
    }

// used for timing portions of the code
// compile with -DTIMEIT to turn this on
#ifdef TIMEIT
    clock_t t0,t1;
#endif
    // process command line
    int LEN = atoi(argv[1]);
    int min_diff = atoi(argv[2]);
    int N = atoi(argv[3]);  // up to 1M, if more change the kernel parameters below
    
    // pool of valid sequences on device
    char* d_pool;    // free later
    CudaSafeCall( cudaMalloc((void**)&d_pool, N*LEN*sizeof(char)) );
    CudaSafeCall( cudaMemset(d_pool, 0, N*LEN*sizeof(char)) );
    
    // current number of valid sequences in pool
    int pool_size = 0;
    
    // pass criteria variables
    int h_pass = 0; // place-holder for d_pass on host
    int* d_pass;    // free later
    CudaSafeCall( cudaMalloc((void**)&d_pass, sizeof(int)) );
    
    // barcodes pool
    char* h_seq = (char*)calloc(CHUNK*LEN, sizeof(char));    // free later
    char* d_seq; // free later
    CudaSafeCall( cudaMalloc((void**)&d_seq, CHUNK*LEN*sizeof(char)) );
    
    // temporary sequence for printing
    char tmp_seq[LEN+1];
    memset(tmp_seq,0,(LEN+1)*sizeof(char));
    
    bool done = false;
    int chunk_counter = 0;
    while (!done) {
        populate_h_seq(h_seq, LEN, CHUNK);
        
        // copy sequences to device
        CudaSafeCall( cudaMemcpy(d_seq, h_seq, CHUNK*LEN*sizeof(char), cudaMemcpyHostToDevice) );
        
        // seed the data
        if (pool_size == 0) {
            pool_size = 1;
            CudaSafeCall( cudaMemcpy(d_pool, h_seq, LEN*sizeof(char), cudaMemcpyHostToDevice) );
        }
        
        // iterate over all the random sequences
        for (int i = 0; i < CHUNK; i++) {
            
            // reset d_pass
            CudaSafeCall( cudaMemset(d_pass, 0, sizeof(int)) );
            
            // launch kernel, initialize for a pool of 1024x1024 sequences
#ifdef TIMEIT
            t0 = clock();
#endif
            check_sequence<<<1024,1024>>>( d_seq+i*LEN, d_pool, d_pass, pool_size, LEN, min_diff );
#ifdef TIMEIT
            t1 = clock();
#endif
            CudaCheckError();
            
            // check d_pass put it in h_pass
            CudaSafeCall( cudaMemcpy(&h_pass, d_pass, sizeof(int), cudaMemcpyDeviceToHost) );
            
            // valid sequence?
            if (h_pass == pool_size) {
                // output
                memcpy(tmp_seq, h_seq+i*LEN, LEN*sizeof(char));
                fprintf(stdout,"%i\t%i\t%s", i + chunk_counter, pool_size, tmp_seq);
#ifdef TIMEIT
                fprintf(stdout, "\t%f", (float)(t1-t0)/CLOCKS_PER_SEC);
#endif
                fprintf(stdout, "\n");
                
                // reached desired nimber of sequences?
                if (pool_size < N) {
                    
                    // copy sequence to the end of the d_pool string and increment pool_size
                    CudaSafeCall( cudaMemcpy(d_pool+pool_size*LEN, h_seq+i*LEN, LEN*sizeof(char), cudaMemcpyHostToDevice) );
                    pool_size++;
                } else {
                    done = true;
                    break;
                }
            }
        }
        chunk_counter += CHUNK;
    }

    // cleanup
    cudaFree(d_pool);
    cudaFree(d_seq);
    cudaFree(d_pass);
    free(h_seq);
    
    return 0;
}