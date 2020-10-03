const char* dgemm_desc = "My awesome dgemm.";

#define GCC_ALN(var,alignment)\
	__builtin_assume_aligned(var,alignment)
#define MEMBD ((int) 64)
void square_dgemm(const int M, const double* restrict A, const double* restrict B, double* restrict C)
//void square_dgemm(const int M, const double *A, const double *B, double *C)
{
    int i, j, k;
   A= GCC_ALN(A,MEMBD);
   B= GCC_ALN(B,MEMBD);
   C= GCC_ALN(C,MEMBD);
    for (j = 0; j < M; ++j) {
        for (k = 0; k < M; ++k) {
            for (i = 0; i < M; ++i){
                C[j*M+i] += A[k*M+i] * B[j*M+k];
	    }
        }
    }
}
