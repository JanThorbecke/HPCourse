#include "par.h"
#include "segyhdr.h"
#include <genfft.h>
#include <assert.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void extr3D(complex *wavelet, float df, float fmin, float fmax, int nfreq, int nx, float dx, int ny, float dy, int oplx, int oply, float dz, int nz, float alpha, float cp, float wfacto, float wfacts, float *data3d, int order, int McC, int opt, int d1op, int oper, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" PULS3D - 3D puls generation.",
" 								",
" PULS3D oper= [optional parameters]",
" 							        ",
" Required parameters:",
" ",
"   oper= .................... Type of 3D-extrapolation ",
"   file_out= ................ Output file with the puls response",
" 							        ",
" Optional parameters:",
" ",
" MEDIUM DEFINITION ",
"   cp=1000 .................. Wave velocity (homogeneous medium)",
"   nx=111 ................... number of x samples",
"   ny=nx .................... number of y samples",
"   dx=10.0 .................. stepsize in spatial direction",
"   dy=dx .................... stepsize in spatial direction",
" SAMPLING AND SOURCE DEFINITION ",
"   file_src=spike ........... source wavelet (overrules dt)",
"   nt=256 ................... number of samples",
"   dt=0.004 ................. stepsize in time-direction ",
"   fmin=5 ................... minimum frequency",
"   fmax=45 .................. maximum frequency",
" EXTRAPOLATION DEFINITION ",
"   dz=dx .................... extrapolation step in W-operators",
"   oplx=19 .................. length of the x-operator",
"   oply=oplx ................ length of the y-operator",
"   alpha=65 ................. maximum angle of interest (degrees)",
"   wfacto=5e-5 .............. weight factor in operator calculation",
"   wfacts=1e-2 .............. weight factor in Series coefficients",
"   perc=5 ................... percentage of smoothing",
"   ndepth=55 ................ number of depth steps",
"   order=13 ................. order in McClellan and Series expansion",
"   McC=1 .................... order of McClellan filter for cos(kr)",
"   opt=0 .................... optimization of the McClellan filter",
"   verbose= ................. 1: verbose option",
"",
"   Options for oper:",
"         - 1  = 2D convolution simple",
"         - 4  = 2D convolution Q4",
"         - 8  = 2D convolution Q8",
"         - 0  = 2D convolution Vector",
" ",
"Jan Thorbecke: janth@xs4all.nl ", 
" ",
NULL
};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	FILE    *fp;
	int     oplx, oply, verbose, ndepth;
	int     j, nx, ny, nt, nfreq, oper, order, nread, nwrite, k;
	int 	sign, McC, opt, d1op, i, n1;
	float   dx, dy, dt, dz, fmin, fmax, df, d1;
	double  t0, t1, t2;
	float   alpha, wfacto, wfacts, cp, perc;
	float 	*data3d, *wavelet, *tmpdata;
	complex *cwave;
    char    *file_src, *file_out;
    segyhdr tr;

/* Reading parameters */

	t0 = wallclock_time();

	initargs(argc, argv);
	requestdoc(1);

	if(!getparint("oper", &oper)) oper=0;
    if(!getparint("verbose", &verbose)) verbose = 0;
    if(!getparstring("file_out", &file_out)){
        if (verbose) fprintf(stderr,"parameter file_out not found, assume pipe\n");
        file_out = NULL;
    }
    if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparfloat("cp", &cp)) cp = 1000;
	if(!getparint("nx", &nx)) nx = 111;
	if(!getparfloat("dx", &dx)) dx = 10.0;
	if(!getparint("ny", &ny)) ny = nx;
	if(!getparfloat("dy", &dy)) dy = dx;
	if(!getparint("nt", &nt)) nt = 256;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getparfloat("dz", &dz)) dz = dx;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparint("oplx", &oplx)) oplx = 19;
	if(!getparint("oply", &oply)) oply = oplx;
	if(!getparfloat("wfacto", &wfacto)) wfacto = 5e-5;
	if(!getparfloat("wfacts", &wfacts)) wfacts = 1e-2;
	if(!getparfloat("perc", &perc)) perc = 5.0;
	if(!getparfloat("fmin", &fmin)) fmin = 5.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45.0;
	if(!getparint("ndepth", &ndepth)) ndepth = 55;
	if(!getparint("opt", &opt)) opt = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;
	d1op = 1;
	order = 1;
	McC = 1;
	opt = 1;

/* check parameters */

	if (verbose) {
		fprintf(stderr, "Minimum frequency to treat ............... : ");
		fprintf(stderr, "%f\n", fmin);
		fprintf(stderr, "Maximum frequency to treat ............... : ");
		fprintf(stderr, "%f\n", fmax);
		fprintf(stderr, "Maximum unaliased frequency .............. : ");
		fprintf(stderr, "%f\n", cp/(2.0*dx*sin(alpha*PI/180.0)));
		fprintf(stderr, "Number of x samples ...................... : ");
		fprintf(stderr, "%d\n", nx);
		fprintf(stderr, "Stepsize in spatial(x) domain ............ : ");
		fprintf(stderr, "%f\n", dx);
		fprintf(stderr, "Number of y samples ...................... : ");
		fprintf(stderr, "%d\n", ny);
		fprintf(stderr, "Stepsize in spatial(y) domain ............ : ");
		fprintf(stderr, "%f\n", dy);
		fprintf(stderr, "Number of time samples ................... : ");
		fprintf(stderr, "%d\n", nt);
		fprintf(stderr, "Stepsize in time domain .................. : ");
		fprintf(stderr, "%f\n", dt);
		fprintf(stderr, "Extrapolation step ....................... : ");
		fprintf(stderr, "%f\n", dz);
		fprintf(stderr, "Operator length in x (odd) ............... : ");
		fprintf(stderr, "%d\n", oplx);
		fprintf(stderr, "Operator length in y (odd) ............... : ");
		fprintf(stderr, "%d\n", oply);
		fprintf(stderr, "Maximum angle of interest ................ : ");
		fprintf(stderr, "%f\n", alpha);
		fprintf(stderr, "P-wave velocity  ......................... : ");
		fprintf(stderr, "%f\n", cp);
		fprintf(stderr, "Number of depth steps .................... : ");
		fprintf(stderr, "%d\n", ndepth);
	}

/* ========================= Opening wavelet file ====================== */

    if (file_src == NULL){
        if(!getparfloat("dt", &dt)) dt = 0.004;
        wavelet = (float *)calloc(nt, sizeof(float));
        wavelet[nt/2] = 1.0;
    }
    else {
        if (verbose) 
			fprintf(stderr,"Reading wavelet from file %s.\n", file_src);

		fp = fopen(file_src, "r");
		if (fp == NULL) fprintf(stderr,"error in opening file %s.\n", file_src);

		nread = fread(&tr, 1, TRCBYTES, fp);
		assert(nread == TRCBYTES);
		n1 = tr.ns;
		d1 = (float)tr.dt/1000000.0;
		assert ( (int)(d1*1000) == (int)(dt*1000) );

		wavelet = (float *)calloc(n1, sizeof(float));
		nread = fread(&wavelet[0], sizeof(float), n1, fp);
		assert(nread == n1);
		fclose(fp);

        if (n1 <= nt) for (i = n1; i < nt; i++) wavelet[i] = 0.0;
        else fprintf(stderr,"file_src has more samples than output\n");
    }

	nfreq = nt/2+1;
	df = 1.0/(nt*dt);
	cwave = (complex *)malloc(nfreq*sizeof(complex));
	sign = -1;
	rc1fft(wavelet, cwave, nt, sign);
	free(wavelet);

	if (verbose) {
		fprintf(stderr, "Frequency step ........................... : ");
		fprintf(stderr, "%f\n", df);
		t1 = wallclock_time();
		fprintf(stderr, "\nComputation Time (seconds)\n");
		fprintf(stderr, "- Reading wavelet forward FFT ... = %f s.\n",t1-t0);
	}

/* 3D extrapolation and imaging */

	data3d = (float *)calloc(nx*ny*ndepth,sizeof(float));

	extr3D(cwave, df, fmin, fmax, nfreq, nx, dx, ny, dy, 
		oplx, oply, dz, ndepth, alpha, cp, wfacto, wfacts, data3d,
		order, McC, opt, d1op, oper, verbose);

	if (verbose) {
		t2 = wallclock_time();
		fprintf(stderr, "- Total Migration time .......... = %f s.\n",t2-t1);
	}

/******** Write data to output file ********/

	free(cwave);
	tmpdata = (float *)malloc(ndepth*sizeof(float));
	memset(&tr, 0, TRCBYTES);
	tr.f1 = -(nx-1)*0.5*dx;
	tr.f2 = -(ny-1)*0.5*dy;
	tr.d1 = dz;
	tr.d2 = dx;
	tr.dt = (int)dz*1000;
	tr.ns = ndepth;
	tr.scalco = -1000;
	tr.fldr = 1;
	tr.trwf = nx;

	fp = fopen(file_out, "w");
    if (fp == NULL) fprintf(stderr,"error in opening file %s.\n", file_out);

	for (j = 0; j < ny; j++) {
		tr.gy = 1000*(tr.f2 + j*dy);
		for (i = 0; i < nx; i++) {
            tr.gx = NINT((tr.f1 + i*dx)*1000);
			fwrite(&tr, 1, TRCBYTES, fp);
			for (k = 0; k < ndepth; k++) tmpdata[k] = data3d[k*nx*ny+j*nx+i];
        	nwrite = fwrite(tmpdata, sizeof(float), ndepth, fp);
			if (nwrite != ndepth) fprintf(stderr,"Error in writing samples to file %s.\n", file_out);
        }
    }
	fclose(fp);

	if (verbose) {
		t1 = wallclock_time();
		fprintf(stderr, "- Writing 3D data ............... = %f s.\n",t1-t2);
		fprintf(stderr, "- Total time .................... = %f s.\n",t1-t0);
	}

	free(data3d);
	free(tmpdata);
	return 0;
}
