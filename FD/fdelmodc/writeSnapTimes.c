#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "par.h"
#include "segyhdr.h"
#include "fdelmod.h"

extern FILE *fopen64 (__const char *__restrict __filename,
                      __const char *__restrict __modes);

FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segyhdr *hdr, float *data, int n, FILE *fp);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeSnapTimes(modPar mod, snaPar sna, int ixsrc, int izsrc, int itime, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose)
{
	FILE    *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss;
	int append, isnap;
	int n1, ibnd, ixs, izs, ize, i, j;
	int ix, iz, ix2, iz2;
	float *snap, sdx;
	segyhdr hdr;

    ibnd = mod.iorder/2-1;
	n1   = mod.naz;
	sdx  = 1.0/mod.dx;

	if (sna.nsnap==0) return 0;

	if ( (((itime-sna.delay) % sna.skipdt)==0) && 
		  (itime >= sna.delay) &&
		  (itime <= sna.delay+(sna.nsnap-1)*sna.skipdt) ) {

		isnap = NINT((itime-sna.delay)/sna.skipdt);
		if (verbose) vmess("Writing snapshot(%d) at time=%.3f", isnap+1, itime*mod.dt);
	
		if (isnap) append=1;
		else append=0;

		if (sna.type.vx)  fpvx  = fileOpen(sna.file_snap, "_svx", append);
		if (sna.type.vz)  fpvz  = fileOpen(sna.file_snap, "_svz", append);
		if (sna.type.p)   fpp   = fileOpen(sna.file_snap, "_sp", append);
		if (sna.type.txx) fptxx = fileOpen(sna.file_snap, "_stxx", append);
		if (sna.type.tzz) fptzz = fileOpen(sna.file_snap, "_stzz", append);
		if (sna.type.txz) fptxz = fileOpen(sna.file_snap, "_stxz", append);
		if (sna.type.pp)  fppp  = fileOpen(sna.file_snap, "_spp", append);
		if (sna.type.ss)  fpss  = fileOpen(sna.file_snap, "_sss", append);
	
		memset(&hdr,0,TRCBYTES);
		hdr.dt     = 1000000*(mod.dt);
		hdr.scalco = 1000;
		hdr.scalel = 1000;
		hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
		hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
		hdr.fldr   = isnap+1;
		hdr.trid   = 1;
		hdr.ns     = sna.nz;
		hdr.trwf   = sna.nx;
		hdr.ntr    = (isnap+1)*sna.nx;
		hdr.f1     = sna.z1*mod.dz+mod.z0;
		hdr.f2     = sna.x1*mod.dx+mod.x0;
		hdr.d1     = mod.dz*sna.skipdz;
		hdr.d2     = mod.dx*sna.skipdx;

/***********************************************************************
* vx velocities have one sample less in x-direction
* vz velocities have one sample less in z-direction
* txz stresses have one sample less in z-direction and x-direction
***********************************************************************/

		snap = (float *)malloc(sna.nz*sizeof(float));

		for (ixs=sna.x1, i=0; ixs<=sna.x2; ixs+=sna.skipdx, i++) {
			hdr.tracf  = i+1;
			hdr.tracl  = isnap*sna.nx+i+1;
			hdr.gx     = 1000*(mod.x0+ixs*mod.dx);
			ix = ixs+ibnd;
			ix2 = ix+1;

			izs = sna.z1+ibnd;
			ize = sna.z2+ibnd;
			if (sna.type.vx) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = vx[ix2*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fpvx);
			}
			if (sna.type.vz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = vz[ix*n1+iz+1];
				}
				traceWrite(&hdr, snap, sna.nz, fpvz);
			}
			if (sna.type.p) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = tzz[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fpp);
			}
			if (sna.type.tzz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = tzz[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fptzz);
			}
			if (sna.type.txx) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = txx[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fptxx);
			}
			if (sna.type.txz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = txz[ix2*n1+iz+1];
				}
				traceWrite(&hdr, snap, sna.nz, fptxz);
			}
			/* calculate divergence of velocity field */
			if (sna.type.pp) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					iz2 = iz+1;
					snap[j] = sdx*((vx[ix2*n1+iz]-vx[ix*n1+iz])+
									(vz[ix*n1+iz2]-vz[ix*n1+iz]));
				}
				traceWrite(&hdr, snap, sna.nz, fppp);
			}
			/* calculate rotation of velocity field */
			if (sna.type.ss) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					iz2 = iz+1;
					snap[j] = sdx*((vx[ix2*n1+iz2]-vx[ix2*n1+iz])-
									(vz[ix2*n1+iz2]-vz[ix*n1+iz2]));
				}
				traceWrite(&hdr, snap, sna.nz, fpss);
			}

		}

		if (sna.type.vx) fclose(fpvx);
		if (sna.type.vz) fclose(fpvz);
		if (sna.type.p) fclose(fpp);
		if (sna.type.txx) fclose(fptxx);
		if (sna.type.tzz) fclose(fptzz);
		if (sna.type.txz) fclose(fptxz);
		if (sna.type.pp) fclose(fppp);
		if (sna.type.ss) fclose(fpss);

		free(snap);
	}

	return 0;
}

