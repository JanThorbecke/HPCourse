#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmod.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *src_nwav, int verbose);

int elastic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int verbose)
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
	float tmpx,tmpz, dp, dvx, dvz;
	float *dvvx, *dvvz;
	int   ix, iz, ixs, izs, izp, ibnd, store;
	int   nx, nz, n1;
	int   ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;
	dvvx = (float *)malloc(n1*sizeof(float));
	dvvz = (float *)malloc(n1*sizeof(float));

	ibnd = mod.iorder/2-1;

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

	/* calculate vx for all grid points except on the virtual boundary*/
	for (ix=ioXx; ix<nx+1; ix++) {
		for (iz=ioXz; iz<nz+1; iz++) {
			vx[ix*n1+iz] += rox[ix*n1+iz]*(
						c1*(txx[ix*n1+iz]     - txx[(ix-1)*n1+iz] +
							txz[ix*n1+iz+1]   - txz[ix*n1+iz])    +
						c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
							txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
		}
	}


	/* calculate vz for all grid points except on the virtual boundary */
	for (ix=ioZx; ix<nx+1; ix++) {
		for (iz=ioZz; iz<nz+1; iz++) {
			vz[ix*n1+iz] += roz[ix*n1+iz]*(
						c1*(tzz[ix*n1+iz]     - tzz[ix*n1+iz-1] +
							txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
						c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
							txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, tzz, txx, txz, rox, roz, src_nwav, verbose);
	}

	/* rigid boundary condition clears velocities on boundaries */
	if (bnd.rig[0]) { /* rigid surface at top */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			vx[ix*n1+ibnd] = 0.0;
			vz[ix*n1+ibnd] = -vz[ix*n1+ibnd+1];
			vz[ix*n1+0]    = -vz[ix*n1+3];
		}
	}
	if (bnd.rig[1]) { /* rigid surface at right */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			vz[(nx+ibnd-1)*n1+iz] = 0.0;
			vx[(nx+ibnd)*n1+iz]   = -vx[(nx+ibnd-1)*n1+iz];
			vx[(nx+2)*n1+iz]      = -vx[(nx-1)*n1+iz];
		}
	}
	if (bnd.rig[2]) { /* rigid surface at bottom */
		for (ix=ibnd; ix<=ibnd+nx-1; ix++) {
			vx[ix*n1+nz+ibnd-1] = 0.0;
			vz[ix*n1+nz+ibnd]   = -vz[ix*n1+nz+ibnd-1];
			vz[ix*n1+nz+2]      = -vz[ix*n1+nz-1];
		}
	}
	if (bnd.rig[3]) { /* rigid surface at left */
		for (iz=ibnd; iz<=ibnd+nz-1; iz++) {
			vz[ibnd*n1+iz] = 0.0;
			vx[ibnd*n1+iz] = -vx[(ibnd-1)*n1+iz];
			vx[1*n1+iz]    = -vx[3*n1+iz];
		}
	}

	/* calculate Txx/tzz for all grid points except on the virtual boundary */
	for (ix=ioPx; ix<nx+1; ix++) {
		for (iz=ioPz; iz<nz+1; iz++) {
			dvvx[iz] = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
					   c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
		}
		for (iz=ioPz; iz<nz+1; iz++) {
			dvvz[iz] = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
					   c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
		}
		for (iz=ioPz; iz<nz+1; iz++) {
			txx[ix*n1+iz] += l2m[ix*n1+iz]*dvvx[iz] + lam[ix*n1+iz]*dvvz[iz];
			tzz[ix*n1+iz] += l2m[ix*n1+iz]*dvvz[iz] + lam[ix*n1+iz]*dvvx[iz];
		}
	}

/*
	for (ix=ioPx; ix<nx+1; ix++) {
		for (iz=ioPz; iz<nz+1; iz++) {
			dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
					   c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
			dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
					   c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
			txx[ix*n1+iz] += l2m[ix*n1+iz]*dvx + lam[ix*n1+iz]*dvz;
			tzz[ix*n1+iz] += l2m[ix*n1+iz]*dvz + lam[ix*n1+iz]*dvx;
		}
	}
*/
/*
	for (ix=ixs-5; ix<ixs+5; ix++) {
		for (iz=izs-1; iz<izs+5; iz++) {
	fprintf(stderr,"ix=%d iz=%d txx=%e tzz=%e\n", ix, iz, txx[ix*n1+iz], tzz[ix*n1+iz]);
		}
	}
*/
	/* calculate Txz for all grid points except on the virtual boundary */
	for (ix=ioTx; ix<nx+1; ix++) {
		for (iz=ioTz; iz<nz+1; iz++) {
			txz[ix*n1+iz] += mul[ix*n1+iz]*(
					c1*(vx[ix*n1+iz]     - vx[ix*n1+iz-1] +
						vz[ix*n1+iz]     - vz[(ix-1)*n1+iz]) +
					c2*(vx[ix*n1+iz+1]   - vx[ix*n1+iz-2] +
						vz[(ix+1)*n1+iz] - vz[(ix-2)*n1+iz]) );
		}
	}

/*
	for (ix=ixs-5; ix<ixs+5; ix++) {
		for (iz=izs-1; iz<izs+5; iz++) {
	fprintf(stderr,"ix=%d iz=%d txz=%e\n", ix, iz, txz[ix*n1+iz]);
		}
	}
*/
	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, tzz, txx, txz, rox, roz, src_nwav, verbose);
	}

/* Free surface: calculate free surface conditions for stresses 
 *     Conditions (for upper boundary):
 *     1. Tzz = 0
 *     2. Txz = 0
 *     3. Txx: remove term with dVz/dz, computed in e2/e4 routines
 *             and add extra term with dVx/dx,
 *             corresponding to free-surface condition for Txx.
 *             In this way, dVz/dz is not needed in computing Txx
 *             on the upper stress free boundary. Other boundaries
 *             are treated similar.
 *             For the 4th order schemes, the whole virtual boundary
 *             must be taken into account in the removal terms, 
 *             because FndCalcFree4 sets
 *             velocities on this boundary!
 *
 *    Compute the velocities on the virtual boundary to make interpolation
 *    possible for receivers. */


	/* check if there are sources placed on the free surface */
	store=0;
	if (src.type==1 || src.type==3 || src.type==4 || src.type==6) {
		ixs = ixsrc + ibnd;
		izs = izsrc + ibnd;
		if (ixs == ibnd) store=1;
		if (ixs == nx+ibnd-1) store=1;
		if (izs == ibnd) store=1;
		if (izs == nz+ibnd-1) store=1;
		if (store) {
			if (src.type==1) {
				tmpx = txx[ixs*n1+izs];
				tmpz = tzz[ixs*n1+izs];
			}
			else if (src.type==3) {
				tmpz = tzz[ixs*n1+izs];
			}
			else if (src.type==4) {
				tmpx = txx[ixs*n1+izs];
			}
			else if (src.type==6) {
				tmpx = vx[ixs*n1+izs];
			}
		}
	}

	if (bnd.free[0]) { /* free surface at top */
		izp = bnd.surface[0];
		for (ix=1; ix<nx+1; ix++) {
			iz = bnd.surface[ix-1];
			if ( izp==iz ) {
				/* clear normal pressure */
				tzz[ix*n1+iz] = 0.0;
				/* assure that txz=0 on boundary by filling virtual boundary */
				txz[ix*n1+iz] = -txz[ix*n1+iz+1] ;
/* This update to Vz is unstable */
//				vz[ix*n1+iz] = vz[ix*n1+iz+1]+ (vx[(ix+1)*n1+iz]-vx[ix*n1+iz])*
//								lam[ix*n1+iz]/l2m[ix*n1+iz];
			}
			izp=iz;
		}
		/* extra line of txz has to be copied */
		izp = bnd.surface[0];
		for (ix=1; ix<nx+1; ix++) {
			iz = bnd.surface[ix-1];
			if ( izp==iz ) txz[ix*n1+iz-1] = -txz[ix*n1+iz+2];
			izp=iz;
		}
		/* calculate txx on top stress-free boundary */
		izp = bnd.surface[0];
		for (ix=1; ix<nx+1; ix++) {
			iz = bnd.surface[ix-1];
			if ( izp==iz ) {
				dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
				dvx = c1*(vx[(ix+1)*n1+iz] - vx[(ix)*n1+iz]) +
				  	c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
				txx[ix*n1+iz] = dvx*dp;
			}
			izp=iz;
		}
		/* if surface has also left or right edges */
		izp = bnd.surface[0];
		for (ix=1; ix<nx+1; ix++) {
			iz = bnd.surface[ix-1];
			if ( izp < iz ) { /* right boundary */
				/* clear normal pressure */
				txx[ix*n1+iz] = 0.0;
				if ( (iz-izp) >= 2 ) { /* VR point */
					/* assure that txz=0 on boundary */
					txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
					txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
					/* calculate tzz on right stress-free boundary */
					dvz = c1*(vz[ix*n1+iz+1] - vz[ix*n1+iz]) +
						c2*(vz[ix*n1+iz+2] - vz[ix*n1+iz-1]);
					dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
					tzz[ix*n1+iz] = dvz*dp;
				}
				else {
//					if (izp) { /* IR point */	
//						txz[ix*n1+iz] = -txz[ix*n1+iz+1] ;
//						txz[ix*n1+iz-1] = -txz[ix*n1+iz+2];
						txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
						txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
						tzz[ix*n1+iz] = 0.0;
//					}
//					else { /* OR point */
						txz[(ix-1)*n1+iz] = 0.0;
						txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
						txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
				vz[ix*n1+iz] = vz[ix*n1+iz+1]+ (vx[(ix+1)*n1+iz]-vx[ix*n1+iz])*
								lam[ix*n1+iz]/l2m[ix*n1+iz];
//					}
				}
				
			} /* end if right */
			if ( izp > iz ) { /* left boundary */
				/* clear normal pressure */
				txx[ix*n1+iz] = 0.0;
				/* assure that txz=0 on boundary */
				txz[(ix-1)*n1+iz] = -txz[ix*n1+iz];
				/* extra line of txz has to be copied */
				txz[(ix-2)*n1+iz] = -txz[(ix+1)*n1+iz] ;
				/* calculate tzz on left stress-free boundary */
				dvz = c1*(vz[ix*n1+iz+1] - vz[ix*n1+iz]) +
					c2*(vz[ix*n1+iz+2] - vz[ix*n1+iz-1]);
				dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
				tzz[ix*n1+iz] = dvz*dp;
			} /* end if left */
			izp=iz;
//			izp=bnd.surface[MAX(ix-2,0)];;
		} /* end ix loop */
	}
	if (bnd.free[1]) { /* free surface at right */
		for (iz=1; iz<nz+1; iz++) {
			/* clear normal pressure */
			txx[nx*n1+iz] = 0.0;
			/* assure that txz=0 on boundary by filling virtual boundary */
			txz[(nx+1)*n1+iz] = -txz[nx*n1+iz];
		}
		/* extra line of txz has to be copied */
		for (iz=2; iz<nz+1; iz++) {
			txz[(nx+2)*n1+iz] = -txz[(nx-1)*n1+iz] ;
		}
		/* calculate tzz on right stress-free boundary */
		ix = nx;
		for (iz=1; iz<nz+1; iz++) {
			dvz = c1*(vz[nx*n1+iz+1] - vz[nx*n1+iz]) +
				  c2*(vz[nx*n1+iz+2] - vz[nx*n1+iz-1]);
			dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
			tzz[nx*n1+iz] = dvz*dp;
		}
	}
	if (bnd.free[2]) { /* free surface at bottom */
		for (ix=1; ix<nx+1; ix++) {
			/* clear normal pressure */
			tzz[ix*n1+nz] = 0.0;
			/* assure that txz=0 on boundary by filling virtual boundary */
			txz[ix*n1+nz+1] = -txz[ix*n1+nz] ;
		}
		/* extra line of txz has to be copied */
		for (ix=2; ix<nx+1; ix++) {
			txz[ix*n1+nz+2] = -txz[ix*n1+nz-1];
		}
		/* calculate txx on bottom stress-free boundary */
		iz = nz;
		for (ix=1; ix<nx+1; ix++) {
			dvx = c1*(vx[(ix+1)*n1+nz-1] - vx[ix*n1+nz-1]) +
				  c2*(vx[(ix+2)*n1+nz-1] - vx[(ix-1)*n1+nz-1]);
			dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
			txx[ix*n1+nz] = dvx*dp;
		}
	}
	if (bnd.free[3]) { /* free surface at left */
		for (iz=1; iz<nz+1; iz++) {
			/* clear normal pressure */
			txx[n1+iz] = 0.0;
			/* assure that txz=0 on boundary by filling virtual boundary */
			txz[n1+iz] = -txz[2*n1+iz];
		}
		/* extra line of txz has to be copied */
		for (iz=2; iz<nz+1; iz++) {
			txz[iz] = -txz[3*n1+iz] ;
		}
		/* calculate tzz on left stress-free boundary */
		for (iz=1; iz<nz+1; iz++) {
			dvz = c1*(vz[n1+iz+1] - vz[n1+iz]) +
				  c2*(vz[n1+iz+2] - vz[n1+iz-1]);
			dp = l2m[ioPx*n1+iz]-lam[ioPx*n1+iz]*lam[ioPx*n1+iz]/l2m[ioPx*n1+iz];
			tzz[n1+iz] = dvz*dp;
		}
	}

	/* restore source positions on the edge */
	if (store) {
		if (src.type==1) {
			txx[ixs*n1+izs]+=tmpx;
			tzz[ixs*n1+izs]+=tmpz;
		}
		else if (src.type==3) {
			tzz[ixs*n1+izs]=tmpz;
		}
		else if (src.type==4) {
			txx[ixs*n1+izs]=tmpx;
		}
		else if (src.type==4) {
			vx[ixs*n1+izs]=tmpx;
		}
	}

	free(dvvx);
	free(dvvz);
      return 0;
}
