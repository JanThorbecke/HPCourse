#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double wallclock_time();
#define NMAX 409600
#define NLOOP 3000

int main(int argc, char *argv[])
{
        float a[NMAX],  b[NMAX],  c[NMAX],  d[NMAX];
        int i,j;
        double t0, t1, t2;

        for (i=0; i<NMAX; i++) {
                a[i] = 0.0;
                b[i] = i;
                c[i] = 2*i;
                d[i] = -0.5*i;
        }

        t0 = wallclock_time();
        #pragma omp parallel private(i,j) shared(a,b,c)
        for (j=0; j<NLOOP; j++) {
        #pragma omp for schedule(static,1)
        for (i=0; i<NMAX; i++) {
                a[i] += b[i]+c[i]*d[i];
        }
        }
        t1 = wallclock_time();
        fprintf(stderr,"shared cache loop = %e seconds\n", t1-t0);
        fprintf(stderr,"a[NMAX/2] = %f\n",a[NMAX/2]);

        for (i=0; i<NMAX; i++) {
                a[i] = 0.0;
                b[i] = i;
                c[i] = 2*i;
                d[i] = -0.5*i;
        }

        t0 = wallclock_time();
        #pragma omp parallel private(i,j) shared(a,b,c)
        for (j=0; j<NLOOP; j++) {
        #pragma omp for schedule(static)
        for (i=0; i<NMAX; i++) {
                a[i] += b[i]+c[i]*d[i];
        }
        }
        t1 = wallclock_time();
        fprintf(stderr,"non-shared loop = %e seconds\n", t1-t0);
        fprintf(stderr,"a[NMAX/2] = %f\n",a[NMAX/2]);

        return 0;
}

double wallclock_time()
{
        struct timeval s_val;
        static struct timeval b_val;
        double time;
        static int base=0;

        gettimeofday(&s_val,0);

        if (!base) {
                b_val = s_val;
                base = 1;
                return 0.0;
        }

        time = (double)(s_val.tv_sec-b_val.tv_sec) + 
                   (double)(1e-6*(s_val.tv_usec-b_val.tv_usec));

        return (double)time;
}


