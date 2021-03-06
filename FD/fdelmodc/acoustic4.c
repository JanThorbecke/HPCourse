#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmod.h"

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *src_nwav, int verbose);

int acoustic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float *vz, float *p, float *rox, float *roz, float *l2m, int verbose)
{
/*********************************************************************
       COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 

  The captial symbols T (=Txx,Tzz) Txz,Vx,Vz represent the actual grid
  The indices ix,iz are related to the T grid, so the capital T 
  symbols represent the actual modelled grid.

  one cel (iz,ix)
       |
       V                              extra column of vx,txz
                                                      |
    -------                                           V
   | txz vz| txz vz  txz vz  txz vz  txz vz  txz vz txz
   |       |      
   | vx  t | vx  t   vx  t   vx  t   vx  t   vx  t  vx
    -------
     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   | 
     txz vz  txz Vz--Txz-Vz--Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  t   vx  t   vx  t   vx  t   vx  t  vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz  <--| 
                                                             |
                                         extra row of txz/vz |

***********************************************************************/

	float c1, c2;
	float tmps;
	int   ix, iz, ixs, izs, ibnd, store;
	int   nx, nz, n1;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;

	ibnd = mod.iorder/2-1;

	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;

	/* calculate vx for all grid points except on the virtual boundary*/
	for (ix=ioXx; ix<nx+1; ix++) {
		for (iz=ioXz; iz<nz+1; iz++) {
			vx[ix*n1+iz] += rox[ix*n1+iz]*(
						c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
						c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]));
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
	for (ix=ioZx; ix<nx+1; ix++) {
		for (iz=ioZz; iz<nz+1; iz++) {
			vz[ix*n1+iz] += roz[ix*n1+iz]*(
						c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
						c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, src_nwav, verbose);
	}

	/* rigid boundary condition clears velocities on boundaries */
	if (bnd.rig[0]) { /* rigid surface at top */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			vx[ix*n1+ibnd] = 0.0;
			vz[ix*n1+ibnd] = -vz[ix*n1+ibnd+1];
			if (mod.iorder == 4) vz[ix*n1+0] = -vz[ix*n1+3];
		}
	}
	if (bnd.rig[1]) { /* rigid surface at right */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			vx[(nx+ibnd)*n1+iz]   = -vx[(nx+ibnd-1)*n1+iz];
			vz[(nx+ibnd-1)*n1+iz] = 0.0;
			if (mod.iorder == 4) vx[(nx+2)*n1+iz] = -vx[(nx-1)*n1+iz];
		}
	}
	if (bnd.rig[2]) { /* rigid surface at bottom */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			vx[ix*n1+nz+ibnd-1] = 0.0;
			vz[ix*n1+nz+ibnd]   = -vz[ix*n1+nz+ibnd-1];
			if (mod.iorder == 4) vz[ix*n1+nz+2] = -vz[ix*n1+nz-1];
		}
	}
	if (bnd.rig[3]) { /* rigid surface at left */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			vx[ibnd*n1+iz] = -vx[(ibnd-1)*n1+iz];
			vz[ibnd*n1+iz] = 0.0;
			if (mod.iorder == 4) vx[1*n1+iz] = -vx[3*n1+iz];
		}
	}


	/* calculate p/tzz for all grid points except on the virtual boundary */
	for (ix=ioPx; ix<nx+1; ix++) {
		for (iz=ioPz; iz<nz+1; iz++) {
			p[ix*n1+iz] += l2m[ix*n1+iz]*(
						c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
						c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]) +
						c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
						c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]));
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, src_nwav, verbose);
	}

/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
	store=0;
	if (src.type==1 || src.type==6) {
		ixs = ixsrc + ibnd;
		izs = izsrc + ibnd;
		if (ixs == ibnd) store=1;
		if (ixs == nx+ibnd-1) store=1;
		if (izs == ibnd) store=1;
		if (izs == nz+ibnd-1) store=1;
		if (store) {
			if (src.type==1) tmps = p[ixs*n1+izs];
			else tmps = vx[ixs*n1+izs];
		}
	}

	if (bnd.free[0]) { /* free surface at top */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			iz = bnd.surface[ix-ibnd];
			p[ix*n1+iz] = 0.0;
		}
	}
	if (bnd.free[1]) { /* free surface at right */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			p[(ibnd+nx-1)*n1+iz] = 0.0;
		}
	}
	if (bnd.free[2]) { /* free surface at bottom */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			p[ix*n1+nz+ibnd-1] = 0.0;
		}
	}
	if (bnd.free[3]) { /* free surface at left */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			p[ibnd*n1+iz] = 0.0;
		}
	}

	/* restore source positions on the edge */
	if (store) {
		if (src.type==1) p[ixs*n1+izs] = tmps;
		else vx[ixs*n1+izs] = tmps;
	}

	return 0;
}
