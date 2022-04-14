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
#define ROWS     1020     // Number of rows in each matrix
#define COLUMNS  1020  // Number of columns in each matrix

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

void multiply_matrices0()
{
	float sum;
	int i,j,k;

	// Multiply the two matrices
	for (i = 0 ; i < ROWS ; i++) {
		for (j = 0 ; j < COLUMNS ; j++) {
			sum = 0.0;
			for (k = 0 ; k < COLUMNS ; k++) {
				if (mask[i][k]==1) {
					sum += matrix_a[i][k] * matrix_b[k][j] ;
				}
			}
			matrix_r[i][j] = sum ;
		}
	}

}

void multiply_matrices1()
{
	float sum = 0.0 ;
	int i,j,k;

	// Multiply the two matrices
	for (i = 0 ; i < ROWS ; i++) {
		for (j = 0 ; j < COLUMNS ; j++) {
			sum = 0.0;
			for (k = 0 ; k < COLUMNS ; k++) {
				sum = sum + matrix_a[i][k] * matrix_b[k][j] ;
//				matrix_r[i][j] = matrix_r[i][j] + matrix_a[i][k] * matrix_b[k][j] ;
			}
			matrix_r[i][j] = sum ;
		}
	}

}

void multiply_matrices2() 
{ 
	int i,j,k;
	// Multiply the two matrices 
	// 
	// Please note that the nesting of the innermost 
	// loops has been changed. The index variables j 
	// and k change the most frequently and the access 
	// pattern through the operand matrices is sequential 
	// using a small stride (one.) This changes improves 
	// access to memory data through the data cache. Data 
	// translation lookaside buffer (DTLB) behavior is 
	// also improved. 

	// Multiply the two matrices 
	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			for (j = 0 ; j < COLUMNS ; j++) { 
				matrix_r[i][j] = matrix_r[i][j] + matrix_a[i][k] * matrix_b[k][j] ; 
			} 
		} 
	} 
} 

void multiply_matrices2b() 
{ 
	int i,j,k;

	// Multiply the two matrices 
	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			for (j = 0 ; j < COLUMNS ; j++) { 
				if (mask[i][k]==1) {
					matrix_r[i][j] = matrix_r[i][j] + matrix_a[i][k] * matrix_b[k][j] ; 
				}
			} 
		} 
	} 
} 

void multiply_matrices2c() 
{ 
	int i,j,k;

	// Multiply the two matrices 
	for (i = 0 ; i < ROWS ; i++) { 
		for (k = 0 ; k < COLUMNS ; k++) { 
			for (j = 0 ; j < COLUMNS ; j++) { 
				matrix_r[i][j] = mask[i][k]*(matrix_r[i][j] + matrix_a[i][k] * matrix_b[k][j]) ; 
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

while (1) {
#pragma omp parallel 
{
	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices0() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 0 =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices1() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 1 =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices2() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 2 =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices2b() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 2b =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices2c() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 2c =%lf\n", t1-t0);

	initialize_matrices() ;
	t0=wallclock_time();
	multiply_matrices3() ;
	t1=wallclock_time();
	fprintf(stderr,"Time multiply_matrices 3 =%lf\n", t1-t0);
}
}

	fclose(result_file) ;

	return( 0 ) ;
}
