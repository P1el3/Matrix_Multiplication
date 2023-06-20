/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* matrix_multiplyer(int N, double *m1, double* m2)
{
	double* m = calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
		{
			for(int k = 0; k < N; k++)
			{
				m[i * N + j] += m1[i* N + k] * m2[k * N + j];
			}
	}
	return m;
}
double* matrix_multiplyer_upper(int N, double *m1, double* m2)
{
	double* m = calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
		{
			for(int k = i; k < N; k++)
			{
				m[i * N + j] += m1[i* N + k] * m2[k * N + j];
			}
	}
	return m;
}

double* matrix_multiplyer_lower(int N, double *m1, double* m2)
{
	double* m = calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
		{
			for(int k = j; k < N; k++)
			{
				m[i * N + j] += m1[i* N + k] * m2[k * N + j];
			}
	}
	return m;
}


double* matrix_transpose(int N, double *m)
{
	double* aux = calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)	
			aux[j * N + i] = m[i * N + j];	
	return aux;
}

double* matrix_addition(int N, double* m1, double* m2)
{
	double* m = calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			m[i * N + j] = m1[i * N + j] + m2[i * N + j];
	return m;
}

double* my_solver(int N, double *A, double* B) {
	//C=A×B×At+Bt×Bt
	
	//AxB
	double *AB;
	AB = matrix_multiplyer_upper(N, A, B);
	//At
	double *At;
	At = matrix_transpose(N, A);
	//ABxAt
	double *ABAt;
	ABAt = matrix_multiplyer_lower(N, AB, At);
	//Bt
	double *Bt;
	Bt = matrix_transpose(N, B);
	//BtBt
	double *BtBt;
	BtBt = matrix_multiplyer(N, Bt, Bt);
	//final ecuation
	double *C;
	C = matrix_addition(N, ABAt, BtBt);

	free(AB);
	free(At);
	free(ABAt);
	free(Bt);
	free(BtBt);
	return C;
}
