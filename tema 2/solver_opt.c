/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */

double* my_solver(int N, double *A, double* B) 
{//C + A x B x At + Bt x Bt
	double* AB = calloc(N * N, sizeof(double));
	double* At = calloc(N * N, sizeof(double));
	double* ABAt = calloc(N * N, sizeof(double));
	double* BtBt = calloc(N * N, sizeof(double));
	double* C = calloc(N * N, sizeof(double));
	//pt inmultiri folosesc funtia din laborator
	for(register int i = 0; i < N; i++)
	{//A x B
  		register double *orig_pa = &A[i * N + i];
  		for(register int j = 0; j < N; j++)
		{
    		register double *pa = orig_pa;
    		register double *pb = &B[i * N + j];
    		register double suma_AB = 0.0;
    		for(register int k = i; k < N; k++)
			{	
      			suma_AB += *pa * *pb;
      			pa++;
      			pb += N;
    		}
    		AB[i * N + j] = suma_AB;
  		}
	}
	for(register int i = 0; i < N; i++)
	{//Bt x Bt
		register double *orig_Bt = &B[i * N];
		for(register int j = 0; j < N; j++)
		{
			register double *bt1 = orig_Bt;
			register double *bt2 = &B[j];
			register double suma_BtBt = 0;	
			for(register int k = 0; k < N; k++)
			{	
				suma_BtBt += *bt1 * *bt2;
				bt1++;
				bt2 += N;
			}
			BtBt[j * N + i] = suma_BtBt;			
		}

	}
	//At
	
	for(register int i = 0; i < N; i++)
	{
		for(register int j = 0; j < N; j++)
		{
			At[i* N + j] = A[j * N + i];
		}
	}

	//AB * At 
	for(register int i = 0; i < N; i++)
	{
  		register double *orig_pa = &AB[i * N];
  		for(register int j = 0; j < N; j++)
		{
    	register double *pa = orig_pa + j;
    	register double *pb = &At[j * N + j];
    	register double suma = 0.0;
		for(register int k = j; k < N; k++)
		{
			suma += *pa * *pb;
			pa++;
			pb += N;
		}	
    	ABAt[i * N + j] = suma;
  		}
	}
	for(register int i = 0; i < N; i++)
	{
		for(register int j = 0; j < N; j++)
		{
			C[i * N + j] = ABAt[i * N + j] + BtBt[i * N + j];
		}
	}
	free(AB);
	free(At);
	free(ABAt);
	free(BtBt);
	return C;
}
