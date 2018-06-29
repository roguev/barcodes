/* ser_barcodes
 * computes a set of DNA sequences where each sequence in the set
 * differs from every other sequence by a minimum of given number
 * of nucleotides. THis is an interesting problem as it becomes more
 * difficult the larger the set becomes. This is the serial implementation
 * in C.
 * Written by Assen Roguev, 2018
 */

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <ctime>

// uncomment or compile with -DTIMEIT to enable benchmarking
//#define TIMEIT

/* gen_rand_sequence generates a random DNA sequence
 * of given length
 * Arguments:
 *  size_t sequence length
 * Returns:
 *  char* to the generated sequence
 */ 
char* gen_rand_sequence( size_t size ) {
    const char charset[] = "ATGC";
    
    // allocate output on the heap
    char* str = (char*)calloc(size, sizeof(char));
    
    // attempting to generate very random sequences
    // so seeding the random generator every time the
    // function is called
    clock_t t = clock();
    srand(t);
    
    if (size) {
        for (size_t n = 0; n < size; n++) {
            int key = rand() % (int) (sizeof charset - 1);
            str[n] = charset[key];
        }
    }
    return str;
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
bool check_barcode( char* b, char** seqs, size_t LEN, size_t pool_size, uint8_t diffs ) {
    uint8_t d = 0;
    size_t passed = 0;
    
    if (pool_size == 0)
        return true;
    
    for (size_t i = 0; i < pool_size; i++) {
        d = 0;
        for (size_t j = 0; j < LEN; j++) {
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
    clock_t t0, t1;
#endif
    
    char* tmp;
    size_t pool_size = 0;
    
    if (argc < 3) {
        fprintf(stderr, "\n%s: Insufficient arguments\n\n", argv[0]);
        exit(1);
    }
    
    size_t LEN = (size_t)atoi(argv[1]);
    uint8_t min_diff = (uint8_t)atoi(argv[2]);
    size_t N = (size_t)atoi(argv[3]);
    
    char** seqs = (char**)calloc(N,sizeof(char*));
    
    int counter = 0;
    while (pool_size < N) {
        tmp = gen_rand_sequence(LEN);
        counter++;
#ifdef TIMEIT
        t0 = clock();
#endif
        if (check_barcode(tmp, seqs, LEN, pool_size, min_diff)) {
#ifdef TIMEIT
            t1 = clock();
#endif
            seqs[pool_size] = tmp;
            pool_size++;
            
            // output total number of sequences tested, current number of
            // sequences identified and the actual sequence
            fprintf(stdout,"%i\t%lu\t%s", counter, pool_size, tmp);
#ifdef TIMEIT
            fprintf(stdout, "\t%f", (float)(t1-t0)/CLOCKS_PER_SEC);
#endif
            fprintf(stdout,"\n");
        }
    }

    // cleanup
    for (size_t i = 0; i < N; i++)
        free(seqs[i]);
    free(seqs);
    
    return 0;
}
