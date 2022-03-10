#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"par.h"
#include"fdelmod.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

int getParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose);

int readModel(modPar mod, float *rox, float *roz, float *l2m, float *lam, float *muu, float *tss, float *tes, float *tep);

int defineSource(wavPar wav, srcPar src, float *src_nwav, int verbose);

int acoustic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float
*vz, float *p, float *rox, float *roz, float *l2m, int verbose);

int viscoacoustic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float
*vz, float *p, float *rox, float *roz, float *l2m, float *tss, float *tep, float *q, int verbose);

int elastic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int verbose);

int viscoelastic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *src_nwav, float *vx, float
*vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, float *ts, float *tep, float
*tes, float *r, float *q, float *p, int verbose);

int getRecTimes(modPar mod, recPar rec, int itime, int rsam, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, float *rec_p, float *rec_pp, float *rec_ss, int verbose);

int writeRec(recPar rec, modPar mod, int ixsrc, int izsrc, int nsam, int ishot, int fileno, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, float *rec_p, float *rec_pp, float *rec_ss, int verbose);

int writeSnapTimes(modPar mod, snaPar sna, int ixsrc, int izsrc, int itime, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose);

/* Self documentation */
char *sdoc[] = {
" ",
" fdelmodc - elastic acoustic finite difference wavefield modeling ",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity file",
"   file_cs= .......... S (cs) velocity file",
"   file_den= ......... density (ro) file",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   file_snap=snap.su . base name for snapshot files",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   dt= ............... read from file_src: if dt==0 then dt= can be used to set it",
"" ,
" OPTIONAL PARAMETERS:",
"   ischeme=3 ......... 1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic",
"   tmod=(nt-1)*dt .... total registration time (nt from file_src)",
"   boundary=1 ........ 1=free, 2=absorbing, 3=rigid, 4=taper",
"   ntaper=0 .......... length of taper at edges of model",
"   tapleft=0 ......... =1: taper left edge of model",
"   tapright=0 ........ =1: taper right edge of model",
"   taptop=0 .......... =1: taper top edge of model",
"   tapbottom=0 ....... =1: taper bottom edge of model",
"   cfree=0 ........... 1=free surface",
//"   grid_dir=f ........ direction of time modeling (r=reverse time)",
"   Qp=15 ............. global Q-value for P-waves in visco-elastic (ischeme=2,4)",
"   file_qp= .......... model file Qp values as function of depth",
"   Qs=Qp ............. global Q-value for S-waves in visco-elastic (ischeme=4)",
"   file_qs= .......... model file Qs values as function of depth",
"   fw=0.5*fmax ....... central frequency for which the Q's are used",
"   sinkdepth=0 ....... grid point below topography (defined bij cp=0.0)",
"   nxmax=512 ......... maximum number of x-positions in model",
"   nzmax=512 ......... maximum number of z-positions in model",
"   ntmax=2048 ........ maximum number of samples in source",
"   verbose=0 ......... silent mode; =1: display info",
"",
" SOURCE DEFINITION:",
"   src_type=1 ........ 1=P 2=Txz 3=Tzz 4=Txx 5=S 6=Fx 7=Fz",
"   src_orient=1 ...... orientation of the source",
"                     - 1=monopole",
"                     - 2=dipole +/-",
"                     - 3=dipole - +",
"                     - 4=dipole +/0/-",
"                     - 5=dipole + -",
"   xsrc=middle ....... x-position of (first) source element",
"   zsrc=zmin ......... z-position of (first) source element",
"   dxsrc=dx .......... x-shift of each next source element",
"   dzsrc=0.0 ......... z-shift of each next source element",
"   nshot=1 ........... number of shots to model",
"   dxshot=dx ......... if nshot > 1: x-shift in shot locations",
"   dzshot=0 .......... if nshot > 1: z-shift in shot locations",
"   fmax=from_src ..... maximum frequency in wavelet",
"   src_multiwav=0 .... use traces in file_src as areal source",
"" ,
" PLANE WAVE SOURCE DEFINITION:",
"   plane_wave=0 ...... model plane wave with nsrc= sources",
"   nsrc=1 ............ number of sources to simulate the plane wave",
"   src_angle=0 ....... angle of plane source array",
"   src_velo=1500 ..... velocity to use in src_angle definition",
"",
" RANDOM SOURCE DEFINITION FOR SEISMIC INTERFEROMTERY:",
"   src_random=0 ...... 1 enables nsrc random sources positions in one modeling",
"   wav_random=1 ...... 1 generates (band limited by fmax) noise signatures ",
"   nsrc=1 ............ number of sources to simulate for one shot",
"   xsrc1=0 ........... left bound for x-position of sources",
"   xsrc2=0 ........... right bound for x-position of sources",
"   zsrc1=0 ........... left bound for z-position of sources",
"   zsrc2=0 ........... right bound for z-position of sources",
"   tsrc1=0.5 ......... begin time interval for random sources",
"   tsrc2=tmod ........ end time interval for random sources",
"   tlength=tsrc2-tsrc1 average duration of random source signal",
"   amplitude=0 ....... distribution of source amplitudes",
"   distribution=0 .... for amplitude and tlenght 0=flat 1=Gaussian ",
"   seed=10 ........... seed for start of random sequence ",
"" ,
" SNAP SHOT SELECTION:",
"   tsnap1=0.1 ........ first snapshot time (s)",
"   tsnap2=0.0 ........ last snapshot time (s)",
"   dtsnap=0.1 ........ snapshot time interval (s)",
"   dxsnap=dx ......... sampling in snapshot in x-direction",
"   dzsnap=dz ......... sampling in snapshot in z-direction",
"   sna_type_p=1 ...... p registration _sp",
"   sna_type_vz=1 ..... Vz registration _svz",
"   sna_type_vx=0 ..... Vx registration _svx",
"   sna_type_txx=0 .... Txx registration _stxx",
"   sna_type_tzz=0 .... Tzz registration _stzz",
"   sna_type_txz=0 .... Txz registration _stxz",
"   sna_type_pp=0 ..... P registration _sP",
"   sna_type_ss=0 ..... S registration _sS",
"   sna_vxvztime=0 .... registration of vx/vx times",
"                       The fd scheme is also staggered in time.",
"                       Time at which vx/vz snapshots are written:",
"                     - 0=previous vx/vz relative to txx/tzz/txz at time t",
"                     - 1=next     vx/vz relative to txx/tzz/txz at time t",
"" ,
" RECEIVER SELECTION:",
"   xrcv1=xmin ........ first x-position of receiver array(s)",
"   xrcv2=xmax ........ last x-position of receiver array(s)",
"   dxrcv=dx .......... x-position increment of receivers",
"   zrcv1=zmin ........ first z-position of receiver array(s)",
"   zrcv2=zrcv1 ....... last z-position of receiver array(s)",
"   dzrcv=0.0 ......... z-position increment of receivers",
"   dtrcv=.004 ........ desired sampling in receiver data (s)",
"   max_nrec=10000 .... maximum number of receivers",
"   largeSUfile=0 ..... writing large SU file (nt > 16000)",
"   rec_ntsam=nt ...... desired number of time samples",
"   rec_delay=0 ....... sample to start recording",
"   dxspread=0 ........ if nshot > 1: x-shift of rcv spread",
"   dzspread=0 ........ if nshot > 1: z-shift of rcv spread",
"   rec_type_p=1 ...... p registration _rp",
"   rec_type_vz=1 ..... Vz registration _rvz",
"   rec_type_vx=0 ..... Vx registration _rvx",
"   rec_type_txx=0 .... Txx registration _rtxx",
"   rec_type_tzz=0 .... Tzz registration _rtzz",
"   rec_type_txz=0 .... Txz registration _rtxz",
"   rec_type_pp=0 ..... P registration _rP",
"   rec_type_ss=0 ..... S registration _rS",
"   rec_int_vx=0  ..... interpolation of Vx receivers",
"                     - 0=Vx->Vx (no interpolation)",
"                     - 1=Vx->Vz",
"                     - 2=Vx->Txx/Tzz",
"   rec_int_vz=0 ...... interpolation of Vz receivers",
"                     - 0=Vz->Vz (no interpolation)",
"                     - 1=Vz->Vx",
"                     - 2=Vz->Txx/Tzz",
"" ,
//"  average=0      - average for 12 points for p-receivers",
"NOTES: For viscoelastic media dispersion and stability are not always",
" quaranteed by the calculated criteria. Specially for Q values smaller than 13",
"",
"      Jan Thorbecke 2008",
"      TU Delft",
"      E-mail: janth@xs4all.nl ",
"  ",
NULL};


int main(int argc, char **argv)
{
	modPar mod;
	recPar rec;
	snaPar sna;
	wavPar wav;
	srcPar src;
	bndPar bnd;
	shotPar shot;
	float *src_nwav;
	float *rox, *roz, *l2m, *lam, *mul;
	float *tss, *tes, *tep, *p, *q, *r;
	float *vx, *vz, *tzz, *txz, *txx;
	float *rec_vx, *rec_vz, *rec_p;
	float *rec_txx, *rec_tzz, *rec_txz;
	float *rec_pp, *rec_ss;
	double t0, t1, t2, t3, tt, tinit;
	size_t size, sizem;
	int n1, n2, nx, nz, nsrc, ix, iz, ir, ishot;
	int ioPx, ioPz, ioTx, ioTz;
	int it0, it1, its, it, fileno, isam;
	int ixsrc, izsrc;
	int verbose;

	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	getParameters(&mod, &rec, &sna, &wav, &src, &shot, &bnd, verbose);

	/* allocate arrays for model parameters */

	n1 = mod.nz+mod.iorder-1;
	n2 = mod.nx+mod.iorder-1;
	mod.naz = n1;
	mod.nax = n2;
	sizem=n1*n2;
	sizem=n2;

	rox = (float *)calloc(sizem,sizeof(float));
	roz = (float *)calloc(sizem,sizeof(float));
	l2m = (float *)calloc(sizem,sizeof(float));
	if (mod.ischeme==2) {
		tss = (float *)calloc(sizem,sizeof(float));
		tep = (float *)calloc(sizem,sizeof(float));
		q = (float *)calloc(sizem,sizeof(float));
	}
	if (mod.ischeme>2) {
		lam = (float *)calloc(sizem,sizeof(float));
		mul = (float *)calloc(sizem,sizeof(float));
	}
	if (mod.ischeme==4) {
		tss = (float *)calloc(sizem,sizeof(float));
		tes = (float *)calloc(sizem,sizeof(float));
		tep = (float *)calloc(sizem,sizeof(float));
		r = (float *)calloc(sizem,sizeof(float));
		p = (float *)calloc(sizem,sizeof(float));
		q = (float *)calloc(sizem,sizeof(float));
	}

	/* read velocity and density files */

	readModel(mod, rox, roz, l2m, lam, mul, tss, tes, tep);

	/* read and/or define source wavelet(s) */

	src_nwav = (float *)calloc(wav.nt*wav.nx,sizeof(float));
	assert(src_nwav != NULL);

	defineSource(wav, src, src_nwav, verbose);

	/* allocate arrays for wavefield and receiver arrays */

	vx  = (float *)malloc(sizem*sizeof(float));
	vz  = (float *)malloc(sizem*sizeof(float));
	tzz = (float *)malloc(sizem*sizeof(float)); /* =P field for acoustic */
	if (mod.ischeme>2) {
		txz = (float *)malloc(sizem*sizeof(float));
		txx = (float *)malloc(sizem*sizeof(float));
	}
	
	size = rec.n*rec.nt;
	if (rec.type.vz)  rec_vz  = (float *)calloc(size,sizeof(float));
	if (rec.type.vx)  rec_vx  = (float *)calloc(size,sizeof(float));
	if (rec.type.p)   rec_p   = (float *)calloc(size,sizeof(float));
	if (rec.type.txx) rec_txx = (float *)calloc(size,sizeof(float));
	if (rec.type.tzz) rec_tzz = (float *)calloc(size,sizeof(float));
	if (rec.type.txz) rec_txz = (float *)calloc(size,sizeof(float));
	if (rec.type.pp)  rec_pp  = (float *)calloc(size,sizeof(float));
	if (rec.type.ss)  rec_ss  = (float *)calloc(size,sizeof(float));

	t1= wallclock_time();
	if (verbose) {
		tinit = t1-t0;
		vmess("*******************************************");
		vmess("************* runtime info ****************");
		vmess("*******************************************");
		vmess("CPU time for intializing arrays and model = %f", tinit);
	}

	/* Sinking source and receiver arrays. 
	   If P-velocity is zero place source and receiver deeper 
	   and store surface position in bnd.surface[ix] */

    ioPx=mod.iorder/2-1;
    ioPz=ioPx;
	ioTx=mod.iorder/2;
	ioTz=ioTx;
	for (ir=0; ir<rec.n; ir++) {
		iz = rec.z[ir];
		ix = rec.x[ir];
		while(l2m[(ix+ioPx)*n1+iz+ioPz] == 0.0) iz++;
		rec.z[ir]=iz+rec.sinkdepth;
		if (verbose>3) vmess("receiver %d at x=%d z=%d", ir, ix, rec.z[ir]);
	}
	for (ishot=0; ishot<shot.n; ishot++) {
		iz = shot.z[ishot];
		ix = shot.x[ishot];
		while(l2m[(ix+ioPx)*n1+iz+ioPz] == 0.0) iz++;
//		if(mul[(ix+ioTx)*n1+iz+ioTz] == 0.0) iz++;
		shot.z[ishot]=iz;
	}
	for (ix=0; ix<mod.nx; ix++) {
		iz = ioPz;
		while(l2m[(ix+ioPx)*n1+iz] == 0.0) iz++;
		bnd.surface[ix] = iz;
		if (verbose>4) vmess("Topgraphy surface x=%.2f z=%.2f", mod.x0+mod.dx*ix, mod.z0+mod.dz*iz);
	}

	/* Outer loop over number of shots */

	for (ishot=0; ishot<shot.n; ishot++) {

		izsrc = shot.z[ishot];
		ixsrc = shot.x[ishot];
		fileno= 0

		memset(vx,0,sizem*sizeof(float));
		memset(vz,0,sizem*sizeof(float));
		memset(tzz,0,sizem*sizeof(float));
		if (mod.ischeme==2) {
			memset(q,0,sizem*sizeof(float));
		}
		if (mod.ischeme>2) {
			memset(txz,0,sizem*sizeof(float));
			memset(txx,0,sizem*sizeof(float));
		}
		if (mod.ischeme==4) {
			memset(r,0,sizem*sizeof(float));
			memset(p,0,sizem*sizeof(float));
			memset(q,0,sizem*sizeof(float));
		}
		if (verbose) {
			vmess("Modeling source at gridpoints ix=%d iz=%d", shot.x[ishot], shot.z[ishot]);
			vmess(" which are actual positions x=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ishot], mod.z0+mod.dz*shot.z[ishot]);
			vmess("Receivers at gridpoint range ix=%d - %d", rec.x[0], rec.x[rec.n-1]);
			vmess(" which are actual positions x=%.2f - %.2f", mod.x0+mod.dx*rec.x[0], mod.x0+mod.dx*rec.x[rec.n-1]);
			vmess("Receivers at gridpoint range iz=%d - %d", rec.z[0], rec.z[rec.n-1]);
			vmess(" which are actual positions z=%.2f - %.2f", mod.z0+mod.dz*rec.z[0], mod.z0+mod.dz*rec.z[rec.n-1]);
		}

		if (mod.grid_dir) {
			it1=0;
			it0=mod.nt;
			its=-1;
		}
		else {
			it0=0;
			it1=mod.nt;
			its=1;
		}

		/* loop over the number of time steps */

		for (it=it0; it<it1; it+=its) {

			switch ( mod.ischeme ) {
				case 1 : /* Acoustic */
					acoustic4(mod, src, wav, bnd, it, ixsrc, izsrc, src_nwav, 
							vx, vz, tzz, rox, roz, l2m, verbose);
					break;
				case 2 : /* Visco-Acoustic */
					viscoacoustic4(mod, src, wav, bnd, it, ixsrc, izsrc, src_nwav, 
							vx, vz, tzz, rox, roz, l2m, tss, tep, q, verbose);
					break;
				case 3 : /* Elastic */
					elastic4(mod, src, wav, bnd, it, ixsrc, izsrc, src_nwav, 
							vx, vz, tzz, txx, txz, rox, roz, l2m, lam, mul, verbose);
					break;
				case 4 : /* Visco-Elastic */
					viscoelastic4(mod, src, wav, bnd, it, ixsrc, izsrc, src_nwav, 
						vx, vz, tzz, txx, txz, rox, roz, l2m, lam, mul, 
						tss, tep, tes, r, q, p, verbose);
					break;
			}

			/* write samples to file if rec.nt samples are calculated */

			if ( (((it-rec.delay) % rec.skipdt)==0) && (it >= rec.delay) ) {
				int writeToFile, itwritten;

				writeToFile = ! ( (((it-rec.delay)/rec.skipdt)+1)%rec.nt );
				itwritten   = fileno*(rec.nt)*rec.skipdt;
				isam        = (it-rec.delay-itwritten)/rec.skipdt;

//				fprintf(stderr,"ToFile = %d isam=%d fileno=%d it=%d itwritten=%d\n", writeToFile,isam, fileno, it, itwritten);

				/* store time at receiver positions */
				getRecTimes(mod, rec, it, isam, &vx, &vz, &tzz, &txx, &txz, 
					rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
					rec_p, rec_pp, rec_ss, verbose);
			
				if (writeToFile && (it != it1-1) ) {
//					fprintf(stderr,"ToFile = %d isam=%d fileno=%d it=%d itwritten=%d\n", writeToFile,isam, fileno, it, itwritten);
					fileno = ( ((it-rec.delay)/rec.skipdt)+1)/rec.nt;
					writeRec(rec, mod, ixsrc, izsrc, isam+1, ishot, fileno,
						rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
						rec_p, rec_pp, rec_ss, verbose);
				}
			}

			/* write snapshots to file */
			writeSnapTimes(mod, sna, ixsrc, izsrc, it, vx, vz, tzz, txx, txz, verbose);

			/* beams */
					
			/* taper the edges of the model */
			taperEdges(mod, bnd, vx, vz, verbose);

			if (verbose) {
				if (it==it0) t2=wallclock_time();
				if (it==it0+500*its) {
					t3=wallclock_time();
					tt=(t3-t2)*((it1-it0)*its)/500.0;
					vmess("Estimated compute time = %.2f s. per shot.",tt);
					vmess("Estimated total compute time = %.2f s.",tinit+shot.n*tt);
				}
			}

		} /* end of loop over time steps */
		}

		/* write output files: receivers or beams */
		if (fileno) fileno++;
		 writeRec(rec, mod, ixsrc, izsrc, isam+1, ishot, fileno,
			rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
			rec_p, rec_pp, rec_ss, verbose);

	} /* end of loop over number of shots */


	t1= wallclock_time();
	if (verbose) {
		vmess("Total compute time FD modelling = %.2f s.", t1-t0);
	}

	/* free arrays */
	
	free(rox);
	free(roz);
	free(l2m);
	free(src_nwav);
	free(vx);
	free(vz);
	free(tzz);
	if (rec.type.vz)  free(rec_vz);
	if (rec.type.vx)  free(rec_vx);
	if (rec.type.p)   free(rec_p);
	if (rec.type.txx) free(rec_txx);
	if (rec.type.tzz) free(rec_tzz);
	if (rec.type.txz) free(rec_txz);
	if (rec.type.pp)  free(rec_pp);
	if (rec.type.ss)  free(rec_ss);

	if (mod.ischeme==2) {
		free(tss);
		free(tep);
		free(q);
	}
	if (mod.ischeme>2) {
		free(lam);
		free(mul);
		free(txz);
		free(txx);
	}
	if (mod.ischeme==4) {
		free(tss);
		free(tes);
		free(tep);
		free(r);
		free(p);
		free(q);
	}


	return 0;
}
