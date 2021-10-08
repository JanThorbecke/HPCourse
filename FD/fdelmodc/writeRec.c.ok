#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "segyhdr.h"
#include "par.h"
#include "SUsegy.h"
#include "fdelmod.h"

extern FILE *fopen64 (__const char *__restrict __filename,
                      __const char *__restrict __modes);

FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segyhdr *hdr, float *data, int n, FILE *fp) ;
void name_ext(char *filename, char *extension);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeRec(recPar rec, modPar mod, int ixsrc, int izsrc, int nsam, int ishot, int fileno, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, float *rec_p, float *rec_pp, float *rec_ss, int verbose)
{
    FILE    *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss;
    int irec;
	int append;
	char number[16], filename[1024];
    segyhdr hdr;
    
	if (!rec.n) return 0;
	if (ishot) append=1;
	else append=0;

	/* check if number of samples fits in hdr.ns (max=16384) */
	strcpy(filename, rec.file_rcv);
	if (fileno) {
		sprintf(number,"_%03d\0",fileno);
		name_ext(filename, number);
	}
	if (verbose>2) vmess("Writing receiver data to file %s", filename);
	if (nsam != rec.nt && verbose) vmess("Number of samples written to last file = %d",nsam);

	memset(&hdr,0,TRCBYTES);
	hdr.dt     = 1000000*(mod.dt*rec.skipdt);
	hdr.scalco = 1000;
	hdr.scalel = 1000;
	hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
	hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
	hdr.fldr   = ishot+1;
	hdr.trid   = 1;
	hdr.ns     = nsam;
	hdr.trwf   = rec.n;
	hdr.ntr    = (ishot+1)*rec.n;
	hdr.f1     = 0.0;
	hdr.d1     = mod.dt*rec.skipdt;
	hdr.d2     = (rec.x[1]-rec.x[0])*mod.dx;
	hdr.f2     = mod.x0+rec.x[0]*mod.dx;

	if (rec.type.vx)  fpvx  = fileOpen(filename, "_rvx", append);
	if (rec.type.vz)  fpvz  = fileOpen(filename, "_rvz", append);
	if (rec.type.p)   fpp   = fileOpen(filename, "_rp", append);
	if (rec.type.txx) fptxx = fileOpen(filename, "_rtxx", append);
	if (rec.type.tzz) fptzz = fileOpen(filename, "_rtzz", append);
	if (rec.type.txz) fptxz = fileOpen(filename, "_rtxz", append);
	if (rec.type.pp)  fppp  = fileOpen(filename, "_rpp", append);
	if (rec.type.ss)  fpss  = fileOpen(filename, "_rss", append);

	for (irec=0; irec<rec.n; irec++) {
		hdr.tracf  = irec+1;
		hdr.tracl  = ishot*rec.n+irec+1;
		hdr.gx     = 1000*(mod.x0+rec.x[irec]*mod.dx);
		hdr.offset = (rec.x[irec]-ixsrc)*mod.dx;

		if (rec.type.vx) {
			traceWrite( &hdr, &rec_vx[irec*rec.nt], nsam, fpvx) ;
		}
		if (rec.type.vz) {
			traceWrite( &hdr, &rec_vz[irec*rec.nt], nsam, fpvz) ;
		}
		if (rec.type.p) {
			traceWrite( &hdr, &rec_p[irec*rec.nt], nsam, fpp) ;
		}
		if (rec.type.txx) {
			traceWrite( &hdr, &rec_txx[irec*rec.nt], nsam, fptxx) ;
		}
		if (rec.type.tzz) {
			traceWrite( &hdr, &rec_tzz[irec*rec.nt], nsam, fptzz) ;
		}
		if (rec.type.txz) {
			traceWrite( &hdr, &rec_txz[irec*rec.nt], nsam, fptxz) ;
		}
		if (rec.type.pp) {
			traceWrite( &hdr, &rec_pp[irec*rec.nt], nsam, fppp) ;
		}
		if (rec.type.ss) {
			traceWrite( &hdr, &rec_ss[irec*rec.nt], nsam, fpss) ;
		}
	}

	if (rec.type.vx) fclose(fpvx);
	if (rec.type.vz) fclose(fpvz);
	if (rec.type.p) fclose(fpp);
	if (rec.type.txx) fclose(fptxx);
	if (rec.type.tzz) fclose(fptzz);
	if (rec.type.txz) fclose(fptxz);
	if (rec.type.pp) fclose(fppp);
	if (rec.type.ss) fclose(fpss);

    return 0;
}

