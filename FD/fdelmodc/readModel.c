#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segyhdr.h"
#include "par.h"
#include "fdelmod.h"

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int readModel(modPar mod, float *rox, float *roz, float *l2m, float *lam, float *muu, float *tss, float *tes, float *tep)
{
    FILE    *fpcp, *fpcs, *fpro;
	FILE    *fpqp=NULL, *fpqs=NULL;
    size_t  nread;
    int i, tracesToDo;
	int n1, ix, iz, nz, nx;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;
	float cp2, cs2, cs11, cs12, cs21, cs22, mul, mu, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, cpx, cpz, bx, bz, fac;
	float *cp, *cs, *ro, *qp, *qs;
	float a, b;
    segyhdr hdr;
    

	/* grid size and start positions for the components */
	nz = mod.nz;
	nx = mod.nx;
	n1 = mod.naz;
	fac = mod.dt/mod.dx;

	/* Vx: rox */
	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
	/* Vz: roz */
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;
	/* Txz: muu */
	ioTx=mod.iorder/2;
	ioTz=ioTx;

/* open files and read first header */

	cp = (float *)malloc(nz*nx*sizeof(float));
   	fpcp = fopen( mod.file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

	ro = (float *)malloc(nz*nx*sizeof(float));
   	fpro = fopen( mod.file_ro, "r" );
   	assert( fpro != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpro);
   	assert(nread == TRCBYTES);

	if (mod.ischeme>2) {
		cs = (float *)malloc(nz*nx*sizeof(float));
		fpcs = fopen( mod.file_cs, "r" );
   		assert( fpcs != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpcs);
   		assert(nread == TRCBYTES);
	}

/* for visco acoustic/elastic media open Q file(s) if given as parameter */

	if (mod.file_qp != NULL && (mod.ischeme==2 || mod.ischeme==4)) {
		qp = (float *)malloc(nz*sizeof(float));
		fpqp = fopen( mod.file_qp, "r" );
   		assert( fpqp != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpqp);
   		assert(nread == TRCBYTES);
	}
	if (mod.file_qs != NULL && mod.ischeme==4) {
		qs = (float *)malloc(nz*sizeof(float));
		fpqs = fopen( mod.file_qs, "r" );
   		assert( fpqs != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpqs);
   		assert(nread == TRCBYTES);
	}


/* read all traces */

	tracesToDo = mod.nx;
	i = 0;
	while (tracesToDo) {
       	nread = fread(&cp[i*nz], sizeof(float), hdr.ns, fpcp);
       	assert (nread == hdr.ns);
       	nread = fread(&ro[i*nz], sizeof(float), hdr.ns, fpro);
       	assert (nread == hdr.ns);
		if (mod.ischeme>2) {
       		nread = fread(&cs[i*nz], sizeof(float), hdr.ns, fpcs);
       		assert (nread == hdr.ns);
		}

/*************************************************************

	Converts the Qp,Qs-value to tau-epsilon and tau-sigma

      tau-sigma    = (sqrt(1.0+(1.0/Qp**2))-(1.0/Qp))/w
      tau-epsilonP = 1.0/(w*w*tau-sigma)
      tau-epsilonS = (1.0+(w*Qs*tau-sigma))/(w*Qs-(w*w*tau-sigma));

*************************************************************/

		/* visco-acoustic */
		if (mod.ischeme==2 || mod.ischeme==4) {
			if (mod.file_qp != NULL) {
       			nread = fread(&qp[0], sizeof(float), nz, fpqp);
       			assert (nread == hdr.ns);
				for (iz=0; iz<nz; iz++) {
					a = sqrt(1.0+(1.0/(qp[iz]*qp[iz]))-(1.0/qp[iz]))/mod.fw;
					b = 1.0/(mod.fw*mod.fw*a);
					tss[(i+ioPx)*n1+iz+ioPz] = 1.0/a;
					tep[(i+ioPx)*n1+iz+ioPz] = b;
				}
			}
			else {
				for (iz=0; iz<nz; iz++) {
					a = sqrt(1.0+(1.0/(mod.Qp*mod.Qp))-(1.0/mod.Qp))/mod.fw;
					b = 1.0/(mod.fw*mod.fw*a);
					tss[(i+ioPx)*n1+iz+ioPz] = 1.0/a;
					tep[(i+ioPx)*n1+iz+ioPz] = b;
				}
			}
		}

		/* visco-elastic */
		if (mod.ischeme==4) {
			if (mod.file_qs != NULL) {
       			nread = fread(&qs[0], sizeof(float), hdr.ns, fpqs);
       			assert (nread == hdr.ns);
				for (iz=0; iz<nz; iz++) {
					a = 1.0/tss[(i+ioPx)*n1+iz+ioPz];
					tes[(i+ioPx)*n1+iz+ioPz] = (1.0+(mod.fw*qs[iz]*a))/(mod.fw*qs[iz]-(mod.fw*mod.fw*a));
				}
			}
			else {
				for (iz=0; iz<nz; iz++) {
					a = 1.0/tss[(i+ioPx)*n1+iz+ioPz];
					tes[(i+ioPx)*n1+iz+ioPz] = (1.0+(mod.fw*mod.Qs*a))/(mod.fw*mod.Qs-(mod.fw*mod.fw*a));
				}
			}
		}

       	nread = fread(&hdr, 1, TRCBYTES, fpcp);
       	if (nread==0) break;
       	nread = fread(&hdr, 1, TRCBYTES, fpro);
       	if (nread==0) break;
		if (mod.ischeme>2) {
       		nread = fread(&hdr, 1, TRCBYTES, fpcs);
       		if (nread==0) break;
		}
		if (mod.file_qp != NULL && (mod.ischeme==2 || mod.ischeme==4)) {
       		nread = fread(&hdr, 1, TRCBYTES, fpqp);
       		if (nread==0) break;
		}
		if (mod.file_qs != NULL && mod.ischeme==4) {
       		nread = fread(&hdr, 1, TRCBYTES, fpqs);
       		if (nread==0) break;
		}
		i++;
	}
   	fclose(fpcp);
   	fclose(fpro);
   	if (mod.ischeme>2) fclose(fpcs);
	if (fpqp != NULL) fclose(fpqp);
	if (fpqs != NULL) fclose(fpqs);

/* check for zero densities */

	for (i=0;i<nz*nx;i++) {
		if (ro[i]==0.0) {
			vwarn("Zero density for trace=%d sample=%d", i/nz, i%nz);
			verr("ERROR zero density is not a valid value, program exit");
		}
	}

/* calculate the grids needed for the FD scheme */

/* the edges of the model */

	if (mod.ischeme>2) { /* Elastic Scheme */
		iz = nz-1;
		for (ix=0;ix<nx-1;ix++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
			cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
			cs2a = cs[(ix+1)*nz+iz]*cs[(ix+1)*nz+iz];
			cs11 = cs2*ro[ix*nz+iz];
			cs12 = cs2*ro[ix*nz+iz];
			cs21 = cs2a*ro[(ix+1)*nz+iz];
			cs22 = cs2a*ro[(ix+1)*nz+iz];
//			cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//			cpz  = cp[ix*nz+iz];

			if (cs11 > 0.0) {
				mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
			}
			else {
				mul = 0.0;
			}
			mu   = cs2*ro[ix*nz+iz];
			lamda2mu = cp2*ro[ix*nz+iz];
			lamda    = lamda2mu - 2*mu;

			bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
			bz = ro[ix*nz+iz];
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
			muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
		}

		ix = nx-1;
		for (iz=0;iz<nz-1;iz++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
			cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
			cs2b = cs[ix*nz+iz+1]*cs[ix*nz+iz+1];
			cs11 = cs2*ro[ix*nz+iz];
			cs12 = cs2b*ro[ix*nz+iz+1];
			cs21 = cs2*ro[ix*nz+iz];
			cs22 = cs2b*ro[ix*nz+iz+1];
//			cpx  = cp[ix*nz+iz];
//			cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1]);

			if (cs11 > 0.0) {
				mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
			}
			else {
				mul = 0.0;
			}
			mu   = cs2*ro[ix*nz+iz];
			lamda2mu = cp2*ro[ix*nz+iz];
			lamda    = lamda2mu - 2*mu;

			bx = ro[ix*nz+iz];
			bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
			muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
		}
		ix=nx-1;
		iz=nz-1;
		cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
		cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
		mu   = cs2*ro[ix*nz+iz];
		lamda2mu = cp2*ro[ix*nz+iz];
		lamda    = lamda2mu - 2*mu;
		bx = ro[ix*nz+iz];
		bz = ro[ix*nz+iz];
		rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
		roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
		l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
		muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;

		for (ix=0;ix<nx-1;ix++) {
			for (iz=0;iz<nz-1;iz++) {
				cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
				cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
				cs2a = cs[(ix+1)*nz+iz]*cs[(ix+1)*nz+iz];
				cs2b = cs[ix*nz+iz+1]*cs[ix*nz+iz+1];
				cs2c = cs[(ix+1)*nz+iz+1]*cs[(ix+1)*nz+iz+1];

/*
Compute harmonic average of mul for accurate and stable fluid-solid interface
see Finite-difference modeling of wave propagation in a fluid-solid configuration 
Robbert van Vossen, Johan O. A. Robertsson, and Chris H. Chapman
*/

				cs11 = cs2*ro[ix*nz+iz];
				cs12 = cs2b*ro[ix*nz+iz+1];
				cs21 = cs2a*ro[ix*nz+iz];
				cs22 = cs2c*ro[ix*nz+iz+1];
//				cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//				cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

				if (cs11 > 0.0) {
					mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
				}
				else {
					mul = 0.0;
				}
				mu   = cs2*ro[ix*nz+iz];
				lamda2mu = cp2*ro[ix*nz+iz];
				lamda    = lamda2mu - 2*mul;
	
				bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
				bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
				rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
				roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
				l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
				lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
				muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
			}
		}
	}
	else { /* Acoustic Scheme */
		iz = nz-1;
		for (ix=0;ix<nx-1;ix++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//			cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//			cpz  = cp[ix*nz+iz];

			lamda2mu = cp2*ro[ix*nz+iz];

			bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
			bz = ro[ix*nz+iz];
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		}

		ix = nx-1;
		for (iz=0;iz<nz-1;iz++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//			cpx  = cp[ix*nz+iz];
//			cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

			lamda2mu = cp2*ro[ix*nz+iz];

			bx = ro[ix*nz+iz];
			bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		}
		ix=nx-1;
		iz=nz-1;
		cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
		lamda2mu = cp2*ro[ix*nz+iz];
		bx = ro[ix*nz+iz];
		bz = ro[ix*nz+iz];
		rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
		roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
		l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;

		for (ix=0;ix<nx-1;ix++) {
			for (iz=0;iz<nz-1;iz++) {
				cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//				cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//				cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

				lamda2mu = cp2*ro[ix*nz+iz];
	
				bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
				bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
				rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
				roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
				l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			}
		}
	}

	/* For topography free surface check for zero-velocity and set rox and roz also to zero */
	for (ix=0;ix<nx;ix++) {
		for (iz=0;iz<nz;iz++) {
			if (l2m[(ix+ioPx)*n1+iz+ioPz]==0.0) {
				rox[(ix+ioXx)*n1+iz+ioXz]=0.0;
				roz[(ix+ioZx)*n1+iz+ioZz]=0.0;
			}
		}
	}

	free(cp);
	free(ro);
   	if (mod.ischeme>2) free(cs);

    return 0;
}


