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
bool check_barcode( char* b, char** seqs, size_t LEN, size_t N_seqs, uint8_t diffs ) {
    uint8_t d = 0;
    size_t passed = 0;
    
    if (N_seqs == 0)
        return true;
    
    for (size_t i = 0; i < N_seqs; i++) {
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
    
    if (passed == N_seqs)
        return true;
    
    return false;
}

/* main takes 3 command line arguments
 *  int desired length
 *  int minimum desired mismatches
 *  int desired number of sequences
 */
int main( int argc, char** argv ) {
    char* tmp;
    size_t current_N = 0;
    
    if (argc < 3) {
        fprintf(stderr, "\n%s: Insufficient arguments\n\n", argv[0]);
        exit(1);
    }
    
    size_t LEN = (size_t)atoi(argv[1]);
    uint8_t min_diff = (uint8_t)atoi(argv[2]);
    size_t N = (size_t)atoi(argv[3]);
    
    char** seqs = (char**)calloc(N,sizeof(char*));
    
    int counter = 0;
    while (current_N < N) {
        tmp = gen_rand_sequence(LEN);
        counter++;
        if (check_barcode(tmp, seqs, LEN, current_N, min_diff)) {
            seqs[current_N] = tmp;
            current_N++;
            
            // output total number of sequences tested, current number of
            // sequences identified and the actual sequence
            fprintf(stdout,"%i\t%lu\t%s\n", counter, current_N, tmp);
        }
    }

    // cleanup
    for (size_t i = 0; i < N; i++)
        free(seqs[i]);
    free(seqs);
    
    return 0;
}
