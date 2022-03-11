#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmod.h"

int getRecTimes(modPar mod, recPar rec, int itime, int isam, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, float *rec_p, float *rec_pp, float *rec_ss, int verbose)
{
	int n1, ibnd;
	int irec, iz, ix2, iz2, ix;

    ibnd = mod.iorder/2-1;
	n1   = mod.naz;

	if (!rec.n) return 0;

/***********************************************************************
* velocity or txz or potential registrations issues:
* rec_x and rec_z are related to actual txx/tzz/p positions.
* offsets from virtual boundaries must be taken into account.
*
* vx velocities have one sample less in x-direction
* vz velocities have one sample less in z-direction
* txz stresses have one sample less in z-direction and x-direction
***********************************************************************/

	for (irec=0; irec<rec.n; irec++) {
		iz = rec.z[irec]+ibnd;
		ix = rec.x[irec]+ibnd;
		iz2 = iz+1;
		ix2 = ix+1;
//			fprintf(stderr,"irec=%d with ix=%d iz=%d\n", irec, ix, iz);
		if (rec.type.p)   rec_p[irec*rec.nt+isam]   = tzz[ix*n1+iz];
		if (rec.type.txx) rec_txx[irec*rec.nt+isam] = txx[ix*n1+iz];
		if (rec.type.tzz) rec_tzz[irec*rec.nt+isam] = tzz[ix*n1+iz];
		if (rec.type.txz) rec_txz[irec*rec.nt+isam] = txz[ix2*n1+iz2];
		if (rec.type.pp) {
				rec_pp[irec*rec.nt+isam] = (vx[ix2*n1+iz]-vx[ix*n1+iz] +
											vz[ix*n1+iz2]-vz[ix*n1+iz])/mod.dx;
		}
		if (rec.type.ss) {
				rec_ss[irec*rec.nt+isam] = (vx[ix2*n1+iz2]-vx[ix2*n1+iz] -
										   (vz[ix2*n1+iz2]-vz[ix*n1+iz2]))/mod.dx;
		}
		if (rec.type.vz) {
/* interpolate vz to vx position to the right and above of vz */
			if (rec.int_vz == 1) {
				rec_vz[irec*rec.nt+isam] = 0.25*(
				vz[ix*n1+iz2]+vz[ix2*n1+iz2]+vz[ix*n1+iz]+vz[ix2*n1+iz]);
			}
/* interpolate vz to Txx/Tzz position by taking the mean of 2 values */
			else if (rec.int_vz == 2) {
				rec_vz[irec*rec.nt+isam] = 0.5*(
				vz[ix*n1+iz2]+vz[ix*n1+iz]);
			}
			else {
				rec_vz[irec*rec.nt+isam] = vz[ix*n1+iz2];
			}
		}
		if (rec.type.vx) {
/* interpolate vx to vz position to the left and below of vx */
			if (rec.int_vx == 1) {
				rec_vx[irec*rec.nt+isam] = 0.25*(
				vx[ix2*n1+iz]+vx[ix2*n1+iz2]+vx[ix*n1+iz]+vx[ix*n1+iz2]);
			}
/* interpolate vx to Txx/Tzz position by taking the mean of 2 values */
			else if (rec.int_vx == 2) {
				rec_vx[irec*rec.nt+isam] = 0.5*(
				vx[ix2*n1+iz]+vx[ix*n1+iz]);
			}
			else {
				rec_vx[irec*rec.nt+isam] = vx[ix2*n1+iz];
			}
		}

	} /* end of irec loop */

	return 0;
}
