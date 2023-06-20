/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) 
{
	double* C = calloc(N * N, sizeof(double));
	double* AB = calloc(N * N, sizeof(double));
	double* ABAt = calloc(N * N, sizeof(double));
	double* Bt = calloc(N * N, sizeof(double));
	//A * B
	cblas_dcopy(N * N, B, 1, AB, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, N, N, 1, A, N, AB, N);
	//A * B * At
	cblas_dcopy(N * N, A, 1, ABAt, 1);
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
		CblasNonUnit, N, N, 1, ABAt, N, AB, N);
	//Bt * Bt
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, 1, B, N, 
		B, N, 0, Bt, N);
	cblas_dcopy(N * N, Bt, 1, C, 1);
	cblas_daxpy(N * N, 1, AB, 1, C, 1);
	
	free(AB);
	free(ABAt);
	free(Bt);
	return C;
}
