#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "par.h"
#include "SUsegy.h"
#include "segy.h"

#define TRCBYTES 240

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define ISODD(n) ((n) & 01)
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2)
{
	FILE    *file_out;
	size_t  nwrite, itrace;
	segy    *hdr;
//	char    *ptr;

/* Read in parameters */

//    ptr = strstr(filename, " ");
//    *ptr = '\0';

	file_out = fopen( filename, "w+" );
	assert( file_out );

	hdr = (segy *)calloc(1,TRCBYTES);
	hdr->ns = n1;
	hdr->dt = NINT(1000000*(d1));
	hdr->d1 = d1;
	hdr->d2 = d2;
	hdr->f1 = f1;
	hdr->f2 = f2;
	hdr->fldr = 1;
	hdr->trwf = n2;

	for (itrace=0; itrace<n2; itrace++) {
		hdr->tracl = itrace+1;
		nwrite = fwrite( hdr, 1, TRCBYTES, file_out );
		assert (nwrite == TRCBYTES);
		nwrite = fwrite( &data[itrace*n1], sizeof(float), n1, file_out );
		assert (nwrite == n1);
	} 
	fclose(file_out);
	free(hdr);

	return 0;
}


int writeSUfile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2)
{
	FILE    *file_out;
	size_t  nwrite, itrace;
	SUsegy  *SUhdr;
	char    *ptr;

/* Read in parameters */

    ptr = strstr(filename, " ");
    *ptr = '\0';

	file_out = fopen( filename, "w+" );
	assert( file_out );

	SUhdr = (SUsegy *)calloc(1,TRCBYTES);
	SUhdr->ns = n1;
	SUhdr->dt = NINT(1000000*(d1));
	SUhdr->d1 = d1;
	SUhdr->d2 = d2;
	SUhdr->f1 = f1;
	SUhdr->f2 = f2;
	SUhdr->fldr = 1;
	SUhdr->trwf = n2;

	for (itrace=0; itrace<n2; itrace++) {
		SUhdr->tracl = itrace+1;
		nwrite = fwrite( SUhdr, 1, TRCBYTES, file_out );
		assert (nwrite == TRCBYTES);
		nwrite = fwrite( &data[itrace*n1], sizeof(float), n1, file_out );
		assert (nwrite == n1);
	} 
	fclose(file_out);
	free(SUhdr);

	return 0;
}

