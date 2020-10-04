#include <stdlib.h>
const char* dgemm_desc = "My awesome dgemm.";

#define GCC_ALN(var,alignment)\
	__builtin_assume_aligned(var,alignment)
#define MEMBD ((int) 64)
void square_dgemm(const int M, const double* restrict A, const double* restrict B, double* restrict C)
//void square_dgemm(const int M, const double *A, const double *B, double *C)
{
	double* restrict newA =  (double*) aligned_alloc(MEMBD,M*M*sizeof(double));
	double* restrict newB =  (double*) aligned_alloc(MEMBD,M*M*sizeof(double));
	double* restrict newC =  (double*) aligned_alloc(MEMBD,M*M*sizeof(double));
	 newA= GCC_ALN(newA,MEMBD);
	 newB= GCC_ALN(newB,MEMBD);
   	 newC= GCC_ALN(newC,MEMBD);


    int i, j, k;
  
    for (j = 0; j < M; ++j) {
            for (i = 0; i < M; ++i){

		newA[j*M+i] = A[j*M+i];
		newB[j*M+i] = B[j*M+i];
                newC[j*M+i] = C[j*M+i];
	    }
    }
    for (j = 0; j < M; ++j) {
        for (k = 0; k < M; ++k) {
            for (i = 0; i < M; ++i){
                newC[j*M+i] += newA[k*M+i] * newB[j*M+k];
	    }
        }
    }
 
    for (j = 0; j < M; ++j) {
            for (i = 0; i < M; ++i){
                C[j*M+i] += newC[j*M+i];
	    }
    }
    free(newA);
    free(newB);
    free(newC);
}
