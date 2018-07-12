 Three programs to compute a set of DNA sequences where each sequence in the set
 differs from every other sequence by a minimum of a given number
 of nucleotides. This is an interesting problem as it becomes more
 difficult the larger the set becomes.
 
 Three implementations are available:
 * ser_barcodes.cc is the serial implementation in C
 * cu_barcodes.cu and cu_barcodes_curand.cu are CUDA C implementations differing by
 the way the random sequences are generated.
 
 The CUDA C code is between 2x and more than 30x faster than the serial C code depending on the
 difficulty of the problem. For example a set of 15,000 sequences of 100 nucleotides differing
 by 60 nucleotides takes about 220 seconds to complete using the CUDA C code (slightly faster using the
 cuRAND variation) and about 1900 seconds on an AMD Phenom(tm) 9750 Quad-Core Processor machine equiped with an 
 NVIDIA GTX 580 graphics card, a performance gain of 8.6x.
