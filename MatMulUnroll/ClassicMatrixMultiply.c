// classic.cpp : "Textbook" implementation of matrix multiply

// Author:  Paul J. Drongowski
// Address: Boston Design Center
//	  Advanced Micro Devices, Inc.
//	  Boxborough, MA 01719
// Date:    20 October 2005
//
// Copyright (c) 2005 Advanced Micro Devices, Inc.

// The purpose of this program is to demonstrate measurement
// and analysis of program performance using AMD CodeAnalyst(tm).
// All engineers are familiar with simple matrix multiplication,
// so this example should be easy to understand.
//
// This implementation of matrix multiplication is a direct
// translation of the "classic" textbook formula for matrix multiply.
// Performance of the classic implementation is affected by an
// inefficient data access pattern, which we should be able to
// identify using CodeAnalyst(TM).

#include <stdlib.h>
#include <stdio.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ROWS     2040     // Number of rows in each matrix
#define COLUMNS  2040  // Number of columns in each matrix

double wallclock_time(void);
static float matrix_a[ROWS][COLUMNS] ;    // Left matrix operand
//static float nep[17] ;
static float matrix_b[ROWS][COLUMNS] ;    // Right matrix operand
//static float nep[17] ; 
static float matrix_r[ROWS][COLUMNS] ;    // Matrix result
//static float nep[17] ; 
static int mask[ROWS][COLUMNS] ;    // mask array

void sgemm_(char *transA, char *transb, int *M, int *N, int *K, float *alpha, float *A, int *lda, float *B, int *ldb, float *beta, float *C, int *ldc);

FILE *result_file ;

void initialize_matrices()
{
	int i,j;

	// Define initial contents of the matrices
	for (i = 0 ; i < ROWS ; i++) {
		for (j = 0 ; j < COLUMNS ; j++) {
			matrix_a[i][j] = (float) rand() / RAND_MAX ;
			matrix_b[i][j] = (float) rand() / RAND_MAX ;
			matrix_r[i][j] = 0.0 ;
			mask[i][j] = NINT(matrix_b[i][j]);
		}
	}
}

void print_result()
{
	int i,j;

	// Print the result matrix
	for (i = 0 ; i < ROWS ; i++) {
		for (j = 0 ; j < COLUMNS ; j++) {
			fprintf(result_file, "%6.4f ", matrix_r[i][j]) ;
		}
		fprintf(result_file, "\n") ;
	}
}


void multiply_matrices2() 
{ 
	int i,j,k;

	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			for (j = 0 ; j < COLUMNS ; j++) { 
				matrix_r[i][j] = matrix_r[i][j] + matrix_a[i][k] * matrix_b[k][j] ; 
			} 
		} 
	} 
} 

void multiply_matricesUnroll1() 
{ 
	int i,j,k;
	float tmp;
	
	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			tmp = matrix_a[i][k];
			for (j = 0 ; j < COLUMNS ; j+=4) { 
				matrix_r[i][j+0] = matrix_r[i][j+0] + tmp*matrix_b[k][j+0];
				matrix_r[i][j+1] = matrix_r[i][j+1] + tmp*matrix_b[k][j+1];
				matrix_r[i][j+2] = matrix_r[i][j+2] + tmp*matrix_b[k][j+2];
				matrix_r[i][j+3] = matrix_r[i][j+3] + tmp*matrix_b[k][j+3];
			} 
		} 
	} 
} 

void multiply_matricesUnroll2() 
{ 
	int i,j,k;
	float tmp0, tmp1, tmp2, tmp3;

	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k+=4) { 
			tmp0 = matrix_a[i][k+0]; 
            tmp1 = matrix_a[i][k+1];
			tmp2 = matrix_a[i][k+2]; 
            tmp3 = matrix_a[i][k+3];
			for (j = 0 ; j < COLUMNS ; j++) { 
				matrix_r[i][j] = matrix_r[i][j] + tmp0*matrix_b[k+0][j]
				+ tmp1*matrix_b[k+1][j]
				+ tmp2*matrix_b[k+2][j]
				+ tmp3*matrix_b[k+3][j];
			} 
		} 
	} 
} 
 
void multiply_matricesUnroll3() 
{ 
	int i,j,k;
	float tmp0, tmp1, tmp2, tmp3;
	
	for (i = 0 ; i < ROWS ; i+=4) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			tmp0 = matrix_a[i+0][k]; 
            tmp1 = matrix_a[i+1][k];
			tmp2 = matrix_a[i+2][k]; 
            tmp3 = matrix_a[i+3][k];
			for (j = 0 ; j < COLUMNS ; j++) { 
				matrix_r[i+0][j] = matrix_r[i+0][j] + tmp0*matrix_b[k][j];
				matrix_r[i+1][j] = matrix_r[i+1][j] + tmp1*matrix_b[k][j];
				matrix_r[i+2][j] = matrix_r[i+2][j] + tmp2*matrix_b[k][j];
				matrix_r[i+3][j] = matrix_r[i+3][j] + tmp3*matrix_b[k][j];
			} 
		} 
	} 
} 


void multiply_matrices3() 
{ 
	char *transa, *transb;
	float beta, alpha;
	int row, col;

    transa = "N";
    transb = "N";
    alpha = 1.0;
    beta = 0.0;
	row = ROWS;
	col = COLUMNS;

	sgemm_(transa, transb, &row, &col, &col, &alpha,
		&matrix_a[0][0], &col, &matrix_b[0][0], &col, &beta, &matrix_r[0][0], &col);

}



int main(int argc, char* argv[])
{
	double t0, t1;
	if ((result_file = fopen("classic.txt", "w")) == NULL) {
		fprintf(stderr, "Couldn't open result file\n") ;
		perror("classic") ;
		return( EXIT_FAILURE ) ;
	}

	fprintf(result_file, "Classic matrix multiplication\n") ;

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices2() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 2 No-unroll          =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matricesUnroll1() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices Unroll 1 Inner-loop  =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matricesUnroll2() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices Unroll 2 Middle-loop =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matricesUnroll3() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices Unroll 3 Outer-loop  =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices3() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 3 SGEMM              =%lf\n", t1-t0);

	fclose(result_file) ;

	return( 0 ) ;
}
