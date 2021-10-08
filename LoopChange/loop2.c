#include<stdio.h>
#include<stdlib.h>

double wallclock_time();
#define N 1000

int main(int argc, char *argv[]) 
{
	static float  A[N][N];
	static float  B[N][N];
	static float  C[N][N];
	double t0, t1;
	int i, j, k;

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			A[i][j] = 0.0;
			B[i][j] = (float)(i*0.0025+j*0.000125);
			C[i][j] = (float)(i*0.00125+j*0.00025);
		}
	}

    t0 = wallclock_time();
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			for (k=0; k<N; k++) {
               A[i][j] = A[i][j]+B[i][k]*C[k][j];
			}
		}
	}
	t1 = wallclock_time();

	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
	fprintf(stderr,"i j k = %f\n", t1-t0);
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			A[i][j] = 0.0;
		}
	}

    t0 = wallclock_time();
    for (i=0; i<N; i++) {
        for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
                A[i][j] = A[i][j]+B[i][k]*C[k][j];
            }
        }
    }
    t1 = wallclock_time();
            
	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
    fprintf(stderr,"i k j = %f\n", t1-t0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0.0;
        }
    }
   
    t0 = wallclock_time();
        for (j=0; j<N; j++) {
    for (i=0; i<N; i++) {
        for (k=0; k<N; k++) {
                A[i][j] = A[i][j]+B[i][k]*C[k][j];
            }
        }
    }
    t1 = wallclock_time();
            
	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
    fprintf(stderr,"j i k = %f\n", t1-t0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0.0;
        }
    }
   
    t0 = wallclock_time();
    for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
    for (i=0; i<N; i++) {
                A[i][j] = A[i][j]+B[i][k]*C[k][j];
            }
        }
    }
    t1 = wallclock_time();
            
	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
    fprintf(stderr,"j k i = %f\n", t1-t0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0.0;
        }
    }
   
    t0 = wallclock_time();
    for (k=0; k<N; k++) {
    for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
                A[i][j] = A[i][j]+B[i][k]*C[k][j];
            }
        }
    }
    t1 = wallclock_time();
            
	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
    fprintf(stderr,"k i j = %f\n", t1-t0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0.0;
        }
    }
   
    t0 = wallclock_time();
    for (k=0; k<N; k++) {
    for (j=0; j<N; j++) {
    for (i=0; i<N; i++) {
                A[i][j] = A[i][j]+B[i][k]*C[k][j];
            }
        }
    }
    t1 = wallclock_time();
            
	fprintf(stderr,"A=%f %f %f\n",A[0][0], A[0][N-1], A[N-1][0]);
    fprintf(stderr,"k j i = %f\n", t1-t0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0.0;
        }
    }
   
	return 0;
}

