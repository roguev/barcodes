/* ser_barcodes.cc
 * computes a set of DNA sequences where each sequence in the set
 * differs from every other sequence by a minimum of a given number
 * of nucleotides. This is an interesting problem as it becomes more
 * difficult the larger the set becomes. This is the serial implementation
 * in C.
 * compile:
 *  g++ -O3 -std=c++11 -DTIMEIT ser_barcodes.cc -o ser_barcodes
 * Written by Assen Roguev, 2018
 */

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <random>

// uncomment or compile with -DTIMEIT to enable benchmarking
//#define TIMEIT

/* gen_rand_sequence generates a random DNA sequence
 * of given length
 * Arguments:
 *  size_t sequence length
 *  std::uniform_int_distribution<int>& random number distrobution to sample from
 *  std::default_random_engine& random number generator
 * Returns:
 *  char* to the generated sequence
 */ 
inline char* gen_rand_sequence( uint16_t LEN, std::uniform_int_distribution<int>& dist, std::default_random_engine& gen ) {
    const char charset[] = "ATGC";
    
    // allocate output on the heap, \0 terminated for printing
    char* seq = (char*)calloc(LEN+1, sizeof(char));
    uint8_t key;
    
    for (uint16_t n = 0; n < LEN; n++) {
        key = (uint8_t)dist(gen);
        seq[n] = charset[key];
    }
    return seq;
}

/* check_sequence checks a given sequence against a pool of sequences
 * Arguments:
 *  char* to the sequence to be checked
 *  char** to an array of sequences in the pool
 *  size_t sequence length
 *  size_t number of sequences in the pool
 *  uint8_t minimum desired mismatches
 * Rerurns:
 * true if the sequence in question passes the test
 * false otherwise
 */
inline bool check_barcode( char* b, char** seqs, uint16_t LEN, uint32_t pool_size, uint16_t diffs ) {
    uint16_t d = 0;
    size_t passed = 0;
    
    if (pool_size == 0)
        return true;
    
    for (uint32_t i = 0; i < pool_size; i++) {
        d = 0;
        for (uint16_t j = 0; j < LEN; j++) {
            if (b[j] != seqs[i][j]) {
                if (++d >= diffs) {
                    passed++;
                    break;
                }
            }
        }
    }
    
    if (passed == pool_size)
        return true;
    
    return false;
}

/* main takes 3 command line arguments
 *  int desired length
 *  int minimum desired mismatches
 *  int desired number of sequences
 */
int main( int argc, char** argv ) {
    #ifdef TIMEIT
    clock_t t0, t1, t2;
    #endif

    char* tmp;
    uint32_t pool_size = 0;
    
    if (argc < 3) {
        fprintf(stderr, "\n%s: Insufficient arguments\n\n", argv[0]);
        exit(1);
    }
    
    uint16_t LEN = (uint16_t)atoi(argv[1]);
    uint16_t min_diff = (uint16_t)atoi(argv[2]);
    uint32_t N = (uint32_t)atoi(argv[3]);
    
    char** seqs = (char**)calloc(N,sizeof(char*));
    
    uint32_t counter = 0;
    
    // for random sequence generation
    std::default_random_engine rng(clock());
    std::uniform_int_distribution<int> rndist(0,3);
    
    #ifdef TIMEIT
    t0 = clock();
    #endif  
    
    while (pool_size < N) {
        tmp = gen_rand_sequence(LEN, rndist, rng);
        counter++;
        
        #ifdef TIMEIT
        t1 = clock();
        #endif
        if (check_barcode(tmp, seqs, LEN, pool_size, min_diff)) {
            #ifdef TIMEIT
            t2 = clock();
            #endif
            seqs[pool_size] = tmp;
            pool_size++;
            
            // output total number of sequences tested, current number of
            // sequences identified and the actual sequence
            #ifdef TIMEIT
            fprintf(stdout, "%i\t%i\t%f\t%f\n", counter, pool_size, (float)(t2-t1)/CLOCKS_PER_SEC, (float)(t2-t0)/CLOCKS_PER_SEC);
            #endif
            
            // output
            fprintf(stdout,"%s\n", tmp);
        }
    }

    // cleanup
    for (uint32_t i = 0; i < N; i++)
        free(seqs[i]);
    free(seqs);
    
    return 0;
}
