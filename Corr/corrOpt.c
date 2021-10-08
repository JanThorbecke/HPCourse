#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "par.h"
#include "segy.h"
double wallclock_time(void);

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define TRCBYTES 240

/*
extern  int __intel_cpu_indicator;
extern  void __intel_cpu_indicator_init(void)  ;
extern  void __intel_new_proc_init_P(void)  ;
extern  void __intel_new_proc_init(void)  ;
__intel_cpu_indicator_init();
*/


/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" CORROPT - Correlation tests in time domain",
"  ",
" CORROPT file_out= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_out= ................ output file",
"  ",
" Optional parameters:",
" ",
"   dt=0.004 ................. time step",
"   ntout=1024 ............... number of output samples",
"   np=16 .................... maximum power of 2 for nt",
"   nrepeat=1 ................ repeat the calculation ",
"   cache=2.0 ................ cache size in MB",
" ",
"      Jan Thorbecke 2006",
"      TU Delft",
"      E-mail: janth@xs4all.nl ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *file_out;
	size_t  nread,  bytes, size, trace_sz, itrace, ntmax, nt;
	int     verbose, i, j, k, l, m, *its, ns, nT,np, seed;
	int 	swap, nb, Nb, nrepeat, ntout, cachesize;
	float   *A, *B, *C, *D, *E, dt, *amp, am, cache;
	float   *A0, *A1, *B0, *B1, *tmp, flops, cr;
	double  t0, t1, t2, tbias;
	char    *filename_out;
	segy    *hdr;


/* Read in parameters */

	initargs(argc,argv);
	requestdoc(0);
	t0 = wallclock_time();
	tbias = wallclock_time()-t0;
	fprintf(stderr,"tbias = %e\n",tbias);
	t0=t1=t2=0.0;

	if(!getparstring("file_out", &filename_out)) filename_out=NULL;
	if(!getparint("ntout", &ntout)) ntout = 1024;
	if(!getparint("np", &np)) np = 16;
	if(!getparint("nrepeat", &nrepeat)) nrepeat = 20;
	if(!getparfloat("cache", &cache)) cache = 2.0;
	cachesize = (int)(cache*1024*1024/sizeof(float));
	fprintf(stderr,"Cache can contain %d floats\n", cachesize);

	if(!getparint("verbose", &verbose)) verbose = 0;

	ntmax = (size_t)MAX(pow(2.,(float)np),cachesize) + ntout;
	fprintf(stderr,"ntmax = %ld\n",ntmax);

	A  = (float *)calloc(ntmax+1,sizeof(float));
	B  = (float *)calloc(ntmax+1,sizeof(float));
	C  = (float *)calloc(ntmax+1,sizeof(float));
	fprintf(stderr,"Memory allocation =%f MB\n",(ntmax*3.0/(1024*1024))*sizeof(float));

	file_out = fopen(filename_out,"w+");

	seed = 10;
	srand(seed);
	for (i=0; i<ntmax; i++) {
		A[i] = -0.5+(rand()/((float)RAND_MAX));
		B[i] = -0.5+(rand()/((float)RAND_MAX));
	}

/* Compute full correlation */

	for (i=1; i<np; i++) {
		nt= (int)pow(2.,(float)i);
		t2=0.0;
		for (k=0; k<nrepeat; k++) {
#pragma ivdep
			for (j=0; j<cachesize; j++) {
				A[j] = 0.0;
				B[j] = 0.0;
			}
			memset(C,0,sizeof(float)*nt);
			t0 = wallclock_time();
#pragma ivdep
			for (l=0; l<ntout; l++) {
#pragma ivdep
/*
                     cr = 0.0;
                     for (j=0; j<nt; j++) {
                         cr += A[j]*B[j+l];
                     }
                    C[l] += cr;

*/
#pragma ivdep
				for (j=0; j<nt; j++) {
					C[l] += A[j]*B[j+l];
				}
			}
			t2 += (wallclock_time() - t0);
		}
		t2 -= tbias;
		flops = ((ntout*2.0*nrepeat)*(nt/(1024.*1024.)))/(t2);
		fprintf(file_out," %ld  %f \n",nt,flops);
		fprintf(stderr,"N=%ld compute time %f with %f Mflop/s\n",nt, t2, flops);
	}
	fprintf(file_out,"\n");


	free(A);
	free(B);
	free(C);
	fclose(file_out);

	return 0;
}
