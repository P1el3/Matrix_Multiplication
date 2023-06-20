#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

double* sup_mat(int N)
{
    double* m = malloc(N * N * sizeof(double));
    int no = 1;
    for(int i = 0; i < N; i++)
    {
        for(int j = i; j < N; j++)
        {
            m[i * N + j] = no;  // fill with random values
            no ++;
        }
        
    }
    
    return m;
}

double* generate_matrix(int N)
{
    double* m = malloc(N * N * sizeof(double));
    int no = 1;
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
        {
            m[i * N + j] = no;
            no ++;
        }
    return m;
}

void print_matrix(int N, double* m)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf("%f " ,m[i * N + j]);
        }
        printf("\n");
    }
}




/*
	C = A * B * At + Bt * Bt 
*/

double* my_solver(int N, double *A, double *B) 
{
	double* C = calloc(N * N, sizeof(double));
	
	//copiez A in A_aux
	double* A_aux = calloc(N * N, sizeof(double));
	cblas_dcopy(N * N, A, 1, A_aux, 1);
	
	//copiez B in B_aux
	double* B_aux = calloc(N * N, sizeof(double));
	cblas_dcopy(N * N, A, 1, B_aux, 1);
	
	//calculez A x B
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
		CblasNonUnit, N, N, 1.0, B, N, A_aux, N);
	
	//calculez A x B x At 
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
		N, N, 1.0, A, N, A_aux, N);
	
	//calculez Bt x Bt
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, 1.0, B, N, 
		B, N, 0.0, B_aux, N);
	
	//copiez A_aux in C
	cblas_dcopy(N * N, A_aux, 1, C, 1);
	
	//calculez A * B * At + Bt * Bt
	cblas_daxpy(N * N, 1, B_aux, 1, C, 1);
	
	free_all(A_aux, B_aux);
	return C;
}


int main()
{
    int i,j,k;
    int N = 3;
    double* m = malloc(N * N * sizeof(double));
    double* m1 = malloc(N * N * sizeof(double));
    double* m2 = malloc(N * N * sizeof(double));
    m1 = sup_mat(N);
    print_matrix(N, m1);
    printf("\n");
    m2 = generate_matrix(N);
    print_matrix(N, m2);
    printf("\n");

    m = my_solver(N, m1, m2);
    print_matrix(N, m);
        
    return 0;
}