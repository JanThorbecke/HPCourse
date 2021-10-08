#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmod.h"

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *src_nwav, int verbose)
{

/*********************************************************************
* 
* description: implements the source on the grid.
* 
* For the acoustic schemes, the source-type must not be txx tzz or txz.
*
**********************************************************************/
	int is0, ibnd;
	int isrc, ix, iz, n1;
	int id1, id2;
	float src_ampl, time, scl, dt;
	static int first=1;

    ibnd = mod.iorder/2-1;
	n1   = mod.naz;
	dt   = mod.dt;

	/* special txz source activated? */

	if (first) {
		iz = izsrc + ibnd;
		if ((iz==0) && (bnd.free[0]) && (src.type==2)) {
			fprintf(stderr," - Using special txz source at top free surface\n");
            if (src.orient != 1) {
				fprintf(stderr,"Only monopole type allowed at surface. Resetted to monopole\n");
				src.orient=1;
			}
		}
		first = 0;
	}
             
/*
* for multiple sources per shot the sources are placed 
* around the central shot position 
* the first source position has an offset in x of is0
*
* itime = 0 corresponds with time=0
* itime = 1 corresponds with time=dt
* src[0] (the first sample) corresponds with time = 0
*/

	is0 = -1*floor((src.n-1)/2);
	for (isrc=0; isrc<src.n; isrc++) {
		src_ampl=0.0;
		/* calculate the source position */
		if (src.random) {
			ix = src.x[isrc] + ibnd;
			iz = src.z[isrc] + ibnd;
		}
		else {
            ix = ixsrc + ibnd + is0 + isrc;
            iz = izsrc + ibnd;
		}
		time = itime*dt - src.tbeg[isrc];
		id1 = floor(time/dt);
		id2 = id1+1;
		/* delay not reached or no samples left in source wavelet? */
		if ( (time < 0.0) || ( (itime*dt) >= src.tend[isrc]) ) continue;

		if (!src.multiwav) { /* only one wavelet */
			src_ampl = src_nwav[id1]*(id2-time/dt) + src_nwav[id2]*(time/dt-id1);
//			src_ampl = src_nwav[id1];
		}
		else { // multi-wavelet sources.
			src_ampl = src_nwav[isrc*wav.nt+id1]*(id2-time/dt) + src_nwav[isrc*wav.nt+id2]*(time/dt-id1);
//			src_ampl = src_nwav[isrc*wav.nt+id1];
		}

//		if (isrc==0) fprintf(stderr,"source %d at %d %d at time %d time=%f tbeg=%f wavampl=%e\n", isrc, ix, iz, id1, time, src.tbeg[isrc],src_nwav[id1]);

		if (src_ampl==0.0) continue;
		if ( ((ix-ibnd)<0) || ((ix-ibnd)>mod.nx) ) continue; /* source outside grid */

		if (verbose>4) {
			fprintf(stderr,"Source %d at grid [ix=%d,iz=%d] at itime %d has value %e\n",isrc, ix,iz, itime, src_ampl);
		}

/* cosine squared windowing to reduce edge effects on one shot arrays */

		if ( (src.n>1) && src.window) {
			if (isrc < src.window) {
				scl = cos(0.5*M_PI*(src.window - isrc)/src.window);
			}
			else if (isrc > src.n-src.window+1) {
				scl = cos(0.5*M_PI*(src.window - (src.n-isrc+1))/src.window);
			}
			src_ampl *= scl*scl;
		}

/* Force source */

		if (src.type == 6) {
			vx[ix*n1+iz] += src_ampl*mod.dx*rox[ix*n1+iz];
			return 0;
		}
		else if (src.type == 7) {
			vz[ix*n1+iz] += src_ampl*mod.dx*roz[ix*n1+iz];
			return 0;
		}

/* Stress source */

		if (mod.ischeme <= 2) { // Acoustic scheme
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient != 1) src_ampl=src_ampl/mod.dx;

				if (src.orient==1) { // monopole
					tzz[ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { // dipole +/-
					tzz[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { // dipole - +
					tzz[ix*n1+iz] += src_ampl;
					tzz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { // dipole +/0/-
					if (iz > ibnd) 
						tzz[ix*n1+iz-1]+= 0.5*src_ampl;
					if (iz < mod.nz+ibnd-1) 
						tzz[ix*n1+iz+1] -= 0.5*src_ampl;
				}
				else if (src.orient==5) { // dipole + -
					tzz[ix*n1+iz] += src_ampl;
					tzz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
		}
		else { // Elastic scheme
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient==1) { // monopole
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { // dipole +/-
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[ix*n1+iz+1] -= src_ampl;
					tzz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { // dipole - +
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[(ix-1)*n1+iz] -= src_ampl;
					tzz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { // dipole +/0/-
					if (iz > ibnd) {
						txx[ix*n1+iz-1]+= 0.5*src_ampl;
						tzz[ix*n1+iz-1]+= 0.5*src_ampl;
					}
					if (iz < mod.nz+ibnd-1) {
						txx[ix*n1+iz+1] -= 0.5*src_ampl;
						tzz[ix*n1+iz+1] -= 0.5*src_ampl;
					}
				}
				else if (src.orient==5) { // dipole + -
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[(ix+1)*n1+iz] -= src_ampl;
					tzz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
			else if (src.type == 2) {
				/* Txz source */
				if ((iz == ibnd) && bnd.free[0]) {
/* offset is not corrected to allow for surface sources at free surface!! */
					txz[ix*n1+iz] += src_ampl;
					txz[(ix+1)*n1+iz] += src_ampl;
				}
				else {
					ix=ix+1; // correct for offset from x boundary
					iz=iz+1; // correct for offset from z boundary
					txz[ix*n1+iz] += src_ampl;
				}
/* possible dipole orientations for a txz source */
				if (src.orient == 2) { // dipole +/-
					txz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient == 3) { // dipole - +
					txz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient == 4) { //  dipole +/O/-
/* correction: subtrace previous value to prevent z-1 values. */
					txz[ix*n1+iz] -= 2.0*src_ampl;
					txz[ix*n1+iz+1] += src_ampl;
				}
				else if (src.orient == 5) { // dipole + -
					txz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
			/* Tzz source */
			else if(src.type == 3) {
				tzz[ix*n1+iz] += src_ampl;
			} 
			/* Txx source */
			else if(src.type == 4) {
				txx[ix*n1+iz] += src_ampl;
			} 

/***********************************************************************
* pure shear source (experimental)
***********************************************************************/
			else if(src.type == 5) {
				ix=ix+1; // correct for offset from x boundary
				iz=iz+1; // correct for offset from z boundary
				if (src.orient == 3 || src.orient == 4) src_ampl = -src_ampl;
				vx[ix*n1+iz] -= src_ampl;
				vz[ix*n1+iz] += src_ampl;
				vx[ix*n1+iz-1] += src_ampl;
				vz[(ix-1)*n1+iz] -= src_ampl;
/* determine second position of dipole */
				if (src.orient == 2) { // dipole +/-
					iz += 1;
				}
				else if (src.orient == 3) { // dipole - +
					ix += 1;
				}
				else if (src.orient == 4) { // dipole +/0/-
					iz += 1;
				}
				else if (src.orient == 5) { // dipole + -
					ix += 1;
				}
				if (src.orient != 1) {
					src_ampl = -src_ampl;
					vx[ix*n1+iz] -= src_ampl;
					vz[ix*n1+iz] += src_ampl;
					vx[ix*n1+iz-1] += src_ampl;
					vz[(ix-1)*n1+iz] -= src_ampl;
				}
			} // src.type
		} // ischeme
	} // loop over isrc

	return 0;
}
