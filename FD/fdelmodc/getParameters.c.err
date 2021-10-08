#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"fdelmod.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

float gaussGen();

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose);

int recvPar(int *rec_x, int *rec_z, int *nrcv, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz);

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);

int getParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose)
{
	int isnapmax1, isnapmax2, isnapmax, sna_nrsna;
	int n1, n2, nx, nz, nsrc, ix, axis;
	int idzsrc, idxsrc;
	int src_ix0, src_iz0, src_ix1, src_iz1;
	float cp_min, cp_max, cs_min, cs_max, ro_min, ro_max;
	float stabfactor,dispfactor, cmin, dt, fmax, scl, wfct, tapfact;
	float xstart,d1,d2,f1,f2,sub_x0,sub_z0;
	float srcendx, srcendz, dx, dz;
	float xsrc, zsrc, dxsrc, dzsrc, dxshot, dzshot, dtshot;
	float dxrcv,dzrcv,dxspread,dzspread;
	float tsnap1, tsnap2, dtsnap, dxsnap, dzsnap, dtrcv;
	float xsnap1, xsnap2, zsnap1, zsnap2, xmax, zmax;
	float xsrc1, xsrc2, zsrc1, zsrc2, tsrc1, tsrc2, tlength;
	float src_angle, src_velo, p, grad2rad;
	int i, max_nrec;
	int boundary, ibnd, cfree;
	int ntaper,tapleft,tapright,taptop,tapbottom;
	int largeSUfile;
	int seed,is,ntraces;
	float rand;
	char *name, tmpname[1024];

	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparint("iorder",&mod->iorder)) mod->iorder=4;
	if (!getparint("ischeme",&mod->ischeme)) mod->ischeme=3;
	if (!getparstring("file_cp",&mod->file_cp)) {
		verr("parameter file_cp required!");
	}
	if (!getparstring("file_den",&mod->file_ro)) {
		verr("parameter file_den required!");
	}
	if (mod->ischeme>2) {
		if (!getparstring("file_cs",&mod->file_cs)) {
			verr("parameter file_cs required!");
		}
	}
	if (!getparstring("file_src",&wav->file_src)) wav->file_src=NULL;
	if (!getparstring("file_snap",&sna->file_snap)) sna->file_snap="snap.su";
	if (!getparstring("file_beam",&rec->file_beam)) rec->file_beam="beam.su";
	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparint("grid_dir",&mod->grid_dir)) mod->grid_dir=0;
	if (!getparint("src_random",&src->random)) src->random=0;
	if (src->random) {
		if (!getparint("wav_random",&wav->random)) wav->random=1;
	}
	else {
		if (!getparint("wav_random",&wav->random)) wav->random=0;
	}
	if (!wav->random) assert (wav->file_src != NULL);
	
	getModelInfo(mod->file_cp, &nz, &nx, &dz, &dx, &sub_z0, &sub_x0, &cp_min, &cp_max, &axis, 1, verbose);
	getModelInfo(mod->file_ro, &n1, &n2, &d1, &d2, &xstart, &xstart, &ro_min, &ro_max, &axis, 0, verbose);
	assert(ro_min != 0.0);
	if (NINT(100*(dx/d2)) != 100) 
		vwarn("dx differs for file_cp and file_den!");
	if (NINT(100*(dz/d1)) != 100) 
		vwarn("dz differs for file_cp and file_den!");
	if (nx != n2) 
		vwarn("nx differs for file_cp and file_den!");
	if (nz != n1) 
		vwarn("nz differs for file_cp and file_den!");

	if (mod->ischeme>2) {
		getModelInfo(mod->file_cs, &n1, &n2, &d1, &d2, &xstart, &xstart, &cs_min, &cs_max, &axis, 1, verbose);
		if (NINT(100*(dx/d2)) != 100) 
			vwarn("dx differs for file_cp and file_cs!");
		if (NINT(100*(dz/d1)) != 100) 
			vwarn("dz differs for file_cp and file_cs!");
		if (nx != n2) 
			vwarn("nx differs for file_cp and file_cs!");
		if (nz != n1) 
			vwarn("nz differs for file_cp and file_cs!");
	}
	mod->dz = dz;
	mod->dx = dx;
	mod->nz = nz;
	mod->nx = nx;
	
	/* define wavelet(s) and model time */

	if (wav->file_src!=NULL) {
		getWaveletInfo(wav->file_src, &wav->nt, &wav->nx, &wav->dt, &d2, &f1, &f2, &fmax, &ntraces, verbose);
		if (wav->dt <= 0.0) {
			vwarn("dt in wavelet (file_src) equal to 0.0 or negative.");
			vwarn("Use parameter dt= to overule dt from file_src.");
		}
		if(!getparfloat("tmod",&mod->tmod)) mod->tmod = (wav->nt-1)*wav->dt;
		if(!getparfloat("dt",&mod->dt)) mod->dt=wav->dt;
		if(!getparfloat("fmax",&wav->fmax)) wav->fmax=fmax;
	}
	else {
		fmax = 50;
		if(!getparfloat("dt",&mod->dt)) verr("dt must be given or use file_src=");
		if(!getparfloat("tmod",&mod->tmod)) verr("tmod must be given");
		if(!getparfloat("fmax",&wav->fmax)) wav->fmax=fmax;
	}
	assert(mod->dt!=0.0);
	mod->nt = NINT(mod->tmod/mod->dt)+1;
	dt = mod->dt;
	if (src->random) {
		if (wav->random) {
			wav->nt=mod->nt;
			wav->dt=mod->dt;
			wav->nx=1;
		}
	}

	if (!getparint("src_type",&src->type)) src->type=1;
	if (!getparint("src_orient",&src->orient)) src->orient=1;
	if (mod->ischeme<=2) {
		if (src->type>1 && src->type<7)
			verr("Invalid src_type for acoustic scheme!");
	}
	if (mod->ischeme==2 || mod->ischeme==4) {
		if (!getparstring("file_qp",&mod->file_qp)) mod->file_qp=NULL;
		if (!getparstring("file_qs",&mod->file_qs)) mod->file_qs=NULL;
		if (!getparfloat("Qp",&mod->Qp)) mod->Qp=1;
		if (!getparfloat("Qs",&mod->Qs)) mod->Qs=mod->Qp;
		if (!getparfloat("fw",&mod->fw)) mod->fw=0.5*wav->fmax;
	}
	assert(src->type > 0);

/* dispersion factor to 10 points per wavelength (2nd order)
   or 5 points per wavelength (4th order) */

	if (mod->iorder == 2) {
		dispfactor=10;
		stabfactor = 1.0/sqrt(2.0);
	}
	else {
		dispfactor = 5;
		stabfactor = 0.606; /* courant number */
	}

	mod->x0 = sub_x0;
	mod->z0 = sub_z0;
	xmax = sub_x0+(nx-1)*dx;
	zmax = sub_z0+(nz-1)*dz;

	if (verbose) {
		vmess("*******************************************");
		vmess("************** general info ***************");
		vmess("*******************************************");
		vmess("tmod    = %f",mod->tmod);
		vmess("ntsam   = %d   dt      = %f",mod->nt, mod->dt);
		if (mod->ischeme == 1) vmess("Acoustic staggered grid, pressure/velocity");
		if (mod->ischeme == 2) vmess("Visco-Acoustic staggered grid, pressure/velocity");
		if (mod->ischeme == 3) vmess("Elastic staggered grid, stress/velocity");
		if (mod->ischeme == 4) vmess("Visco-Elastic staggered grid, stress/velocity");
		if (mod->grid_dir) vmess("Time reversed modelling");
		else vmess("Forward modelling");
		vmess("*******************************************");
		vmess("*************** model info ****************");
		vmess("*******************************************");
		vmess("nz      = %8d   nx      = %8d", nz, nx);
		vmess("dz      = %8.4f   dx      = %8.4f", dz, dx);
		vmess("zmin    = %8.4f   zmax    = %8.4f", sub_z0, zmax);
		vmess("xmin    = %8.4f   xmax    = %8.4f", sub_x0, xmax);
		vmess("min(cp) = %9.3f  max(cp) = %9.3f", cp_min, cp_max);
		if (mod->ischeme>2) vmess("min(cs) = %9.3f  max(cs) = %9.3f", cs_min, cs_max);
		vmess("min(ro) = %9.3f  max(ro) = %9.3f", ro_min, ro_max);
		if (mod->ischeme==2 || mod->ischeme==4) {
			if (mod->file_qp!=NULL) vmess("Qp from file %s   ", mod->file_qp);
			else vmess("Qp      = %9.3f   ", mod->Qp);
			vmess("at freq = %5.3f", mod->fw);
		}
		if (mod->ischeme==4) {
			if (mod->file_qs!=NULL) vmess("Qs from file %s   ", mod->file_qs);
			else vmess("Qs      = %9.3f ", mod->Qs);
			vmess("at freq = %5.3f", mod->fw);
		}

		if (mod->ischeme <= 2) {
			cmin = cp_min;
		}
		else {
			cmin = cs_min; 
			if ( (cmin<1e-20) || (cp_min<cs_min) ) cmin=cp_min;
		}
		vmess("*******************************************");
		vmess("************* stability info **************");
		vmess("*******************************************");
		vmess("Dispersion criterion is %3d points per wavelength: ", NINT(dispfactor));
		vmess(" ====> wavelength > %f m [dx*disp]", dx*dispfactor);
//		vmess("The minimum velocity in the model is %f",cmin);
//		vmess("Hence, for acceptable grid-dispersion the maximum");
		vmess("The maximum frequency in source wavelet must be:");
		vmess(" ====> frequency < %f Hz. [Cmin/dx*disp]", cmin/(dx*dispfactor));
		vmess("Stability criterion for current settings: ");
		vmess(" ====> Cp < %f m/s [dx*disp/dt]", dx*stabfactor/dt);
//		vmess("With dt = %f  maximum velocity = %f",dt, dx*stabfactor/dt);
		vmess("For current model:");
		vmess(" Maximum velocity = %f dtmax = %e", cp_max,dx*stabfactor/cp_max);
		if (wav->file_src != NULL) vmess(" For wavelet(s) in file_src fmax = %f", fmax);
	}

	/* Check stability and dispersion setting */

	if (cp_max > dx*stabfactor/dt) {
		vwarn("***************** !!! *********************");
		vwarn("According to the input file , the maximum P-wave velocity");
		vwarn("in the current model is %f !!", cp_max);
		vwarn("Hence, adjust either dx, dt or lower the maximum velocity.");
		vwarn("***************** !!! *********************");
		verr("********* leaving program *********");
	}
	if (wav->fmax > cmin/(dx*dispfactor)) {
		vwarn("***************** !!! *********************");
		vwarn("The maximum frequency in the source wavelet is");
		vwarn("%.3f for stable modeling fmax < %.3f ", wav->fmax, cmin/(dx*dispfactor));
		vwarn("Hence, adjust either dx, dt or lower the maximum frequency.");
		vwarn("***************** !!! *********************");
		verr("********* leaving program *********");
	}

	/* define boundaries */

	if (!getparint("boundary",&boundary)) boundary=1;
	if (!getparint("ntaper",&ntaper)) ntaper=0;
	if (!getparint("cfree",&cfree)) cfree=0;
	
	for (ibnd=0;ibnd<4;ibnd++) {
		if (boundary == 1) {
			bnd->free[ibnd]=1;
			bnd->rig[ibnd]=0;
			bnd->tap[ibnd]=0;
		}
		else if (boundary == 3) {
			bnd->free[ibnd]=0;
			bnd->rig[ibnd]=1;
			bnd->tap[ibnd]=0;
		}
		else if (boundary == 4) {
			bnd->free[ibnd]=0;
			bnd->rig[ibnd]=0;
			bnd->tap[ibnd]=ntaper;
		}
	}
	if (!getparint("tapleft",&tapleft)) tapleft=0;
	if (!getparint("tapright",&tapright)) tapright=0;
	if (!getparint("taptop",&taptop)) taptop=0;
	if (!getparint("tapbottom",&tapbottom)) tapbottom=0;

	if (tapleft) {
		bnd->free[3]=0;
		bnd->rig[3]=0;
		bnd->tap[3]=ntaper;
	}
	else {
		bnd->tap[3]=0;
		bnd->free[3]=1;
	}
	if (tapright) {
		bnd->free[1]=0;
		bnd->rig[1]=0;
		bnd->tap[1]=ntaper;
	}
	else {
		bnd->tap[1]=0;
		bnd->free[1]=1;
	}
	
	if (taptop) {
		bnd->free[0]=0;
		bnd->rig[0]=0;
		bnd->tap[0]=ntaper;
	}
	else {
		bnd->tap[0]=0;
		bnd->free[0]=1;
	}
	if (tapbottom) {
		bnd->free[2]=0;
		bnd->rig[2]=0;
		bnd->tap[2]=ntaper;
	}
	else {
		bnd->tap[2]=0;
		bnd->free[2]=1;
	}
	
	if (cfree) {
		bnd->free[0]=1;
		bnd->rig[0]=0;
		bnd->tap[0]=0;
	}
	if (ntaper) {
		bnd->tapx = (float *)malloc(ntaper*sizeof(float));
		bnd->tapz = (float *)malloc(ntaper*sizeof(float));
		tapfact = 0.15;
		scl = tapfact/((float)ntaper);
		for (i=0; i<ntaper; i++) {
			wfct = (scl*i);
			bnd->tapx[i] = exp(-(wfct*wfct));

			wfct = (scl*(i+0.5));
			bnd->tapz[i] = exp(-(wfct*wfct));
//			fprintf(stderr,"i=%d old=%f new=%f\n", i, bnd->tapz[i], bnd->tapx[i]);
		}
/*
		scl = 20.0/ntaper;
		for (i=0; i<ntaper; i++) {
			wfct = (0.015*scl*(i));
			bnd->tapx[i] = exp(-(wfct*wfct));
			wfct = (0.015*scl*(i+0.5));
			bnd->tapz[i] = exp(-(wfct*wfct));
		}
*/
	}

	/* check for topography surface */
	bnd->surface = (int *)malloc((nx+nz)*sizeof(int));
	for (ix=0; ix<nx; ix++) {
		bnd->surface[ix] = 1;
	}

	if (verbose) {
		vmess("*******************************************");
		vmess("************* boundary info ***************");
		vmess("*******************************************");
		for (ibnd=0;ibnd<4;ibnd++) {
			if (ibnd==0) name="Top";
			if (ibnd==1) name="Right";
			if (ibnd==2) name="Bottom";
			if (ibnd==3) name="Left";
			if (bnd->free[ibnd]) vmess("%6s boundary  : free",name);
			if (bnd->rig[ibnd]) vmess("%6s boundary  : rigid",name);
			if (bnd->tap[ibnd]) vmess("%6s boundary  : tapered %d points",name,ntaper);
		}
	}

	/* define number of sources to model */

	if (!getparfloat("xsrc",&xsrc)) xsrc=sub_x0+((nx-1)*dx)/2;
	if (!getparfloat("zsrc",&zsrc)) zsrc=sub_z0;
	if (!getparfloat("dxsrc",&dxsrc)) dxsrc=dx;
	if (!getparfloat("dzsrc",&dzsrc)) dzsrc=0.0;
	if (!getparint("nshot",&shot->n)) shot->n=1;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dzshot",&dzshot)) dzshot=0.0;

	if (shot->n>1) {
		idxsrc=MAX(0,NINT(dxsrc/dx));
		idzsrc=MAX(0,NINT(dzsrc/dz));
	}
	else {
		idxsrc=0.0;
		idzsrc=0.0;
	}
	src_ix0=MAX(0,NINT((xsrc-sub_x0)/dx));
	src_ix0=MIN(src_ix0,nx);
	src_iz0=MAX(0,NINT((zsrc-sub_z0)/dz));
	src_iz0=MIN(src_iz0,nz);
	srcendx=(shot->n-1)*dxshot+xsrc;
	srcendz=(shot->n-1)*dzshot+zsrc;
	src_ix1=MAX(0,NINT((srcendx-sub_x0)/dx));
	src_ix1=MIN(src_ix1,nx);
	src_iz1=MAX(0,NINT((srcendz-sub_z0)/dz));
	src_iz1=MIN(src_iz1,nz);

	shot->x = (int *)calloc(shot->n,sizeof(int));
	shot->z = (int *)calloc(shot->n,sizeof(int));
	for (is=0; is<shot->n; is++) {
		shot->x[is] = src_ix0+is*idxsrc;
		shot->z[is] = src_iz0+is*idzsrc;
		if (shot->x[is] > nx-1) shot->n = is-1;
		if (shot->z[is] > nz-1) shot->n = is-1;
	}

	/* define source array */

	if (!getparint("plane_wave",&src->plane)) src->plane=0;
	if (!getparint("src_window",&src->window)) src->window=0;
	if (!getparfloat("src_angle",&src_angle)) src_angle=0.;
	if (!getparfloat("src_velo",&src_velo)) src_velo=1500.;
	if (!getparint("distribution",&src->distribution)) src->distribution=0;
	if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (src->random) {
		if (!getparint("nsrc",&nsrc)) nsrc=1;
		if (!getparint("seed",&seed)) seed=10;
		if (!getparfloat("xsrc1", &xsrc1)) xsrc1=sub_x0;
		if (!getparfloat("xsrc2", &xsrc2)) xsrc2=xmax;
		if (!getparfloat("zsrc1", &zsrc1)) zsrc1=sub_z0;
		if (!getparfloat("zsrc2", &zsrc2)) zsrc2=zmax;
		if (!getparfloat("tsrc1", &tsrc1)) tsrc1=0.5;
		if (!getparfloat("tsrc2", &tsrc2)) tsrc2=mod->tmod;
		tsrc2  = MIN(tsrc2, mod->tmod);
		if (!getparfloat("tlength", &tlength)) tlength=tsrc2-tsrc1;
		dxshot = xsrc2-xsrc1;
		dzshot = zsrc2-zsrc1;
		dtshot = tsrc2-tsrc1;
		if (wav->random) src->multiwav = 1;
		if (wav->random) wav->nx = nsrc;
		if (wav->random) wav->nt = NINT(tlength/mod->dt)+1;
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		srand48(seed);
		for (is=0; is<nsrc; is++) {
			rand = (float)drand48();
			src->x[is] = NINT((xsrc1+rand*dxshot-sub_x0)/dx);
			rand = (float)drand48();
			src->z[is] = NINT((zsrc1+rand*dzshot-sub_z0)/dz);
			rand = (float)drand48();
			src->tbeg[is] = tsrc1+rand*(dtshot);
			if (wav->random) {
				if (src->distribution) rand = fabsf(tlength+gaussGen()*tlength);
				else rand = (float)drand48()*tlength;
				src->tend[is] = MIN(src->tbeg[is]+rand, tsrc2);
			}
			else {
				src->tend[is] = (wav->nt-1)*wav->dt;
			}
			if (verbose>3) {
				vmess("Random xsrc=%f zsrc=%f src_tbeg=%f src_tend=%f",src->x[is]*dx, src->z[is]*dz, src->tbeg[is], src->tend[is]);
			}
		}

/* write time and length os source signals */

		if (verbose>3) {
			float *dum;
			dum = (float *)calloc(mod->nt, sizeof(float));
			for (is=0; is<nsrc; is++) {
				dum[(int)floor(src->tbeg[is]/mod->dt)] = src->tend[is]-src->tbeg[is];
			}
			FILE *fp;
			sprintf(tmpname,"srcTimeLengthN=%d.bin\0",mod->nt);
			fp = fopen(tmpname, "w+");
			fwrite(dum, sizeof(float), mod->nt, fp);
			fclose(fp);
			free(dum);
		}

		/* write velocity field with positions of the sources */
		if (verbose>3) {
			float *dum;
			dum = (float *)calloc(nx*nz, sizeof(float));
			for (is=0; is<nsrc; is++) {
//				dum[src->x[is]*nz+src->z[is]] = src->tbeg[is];
				dum[src->x[is]*nz+src->z[is]] = 1.0;
			}
			writesufile("randomsrcPositions.su", dum, nz, nx, sub_z0, sub_x0, dz, dx);
			free(dum);
		}
	}
	else {
		if (src->plane) { if (!getparint("nsrc",&nsrc)) nsrc=1;}
		else nsrc=1;

		if (nsrc > nx) {
			vwarn("Number of sources used in plane wave is larger than ");
			vwarn("number of gridpoints in X. Plane wave will be clipped to the edges of the model");
		}

	/* for a source defined on mutliple gridpoint calculate p delay factor */

		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		grad2rad = 17.453292e-3;
		p = sin(src_angle*grad2rad)/src_velo;
		if (p < 0.0) {
			for (is=0; is<nsrc; is++) {
				src->tbeg[is] = fabs((nsrc-is-1)*dx*p);
			}
		}
		else {
			for (is=0; is<nsrc; is++) {
				src->tbeg[is] = is*dx*p;
			}
		}
		for (is=0; is<nsrc; is++) {
			src->tend[is] = (wav->nt-1)*wav->dt;
		}
	}

	if (src->multiwav) {
		if (wav->nx != nsrc) {
			vwarn("src_multiwav has been defined but number of traces in");
			vwarn("file_src = %d is not equal to nsrc = %d", wav->nx, nsrc);
			vwarn("last trace in file_src will be repeated.");
		}
		else {
			vmess("Using all traces in file_src for areal shot");
		}
	}
	src->n=nsrc;

	if (verbose) {
		vmess("*******************************************");
		vmess("************* wavelet info ****************");
		vmess("*******************************************");
		vmess("src_nt   = %6d   src_nx      = %d", wav->nt, wav->nx);
		vmess("src_type = %6d   src_orient  = %d", src->type, src->orient);
		vmess("fmax     = %8.2f", fmax);
		fprintf(stderr,"    %s: Source type         : ",xargv[0]);
		switch ( src->type ) {
			case 1 : fprintf(stderr,"P "); break;
			case 2 : fprintf(stderr,"Txz "); break;
			case 3 : fprintf(stderr,"Tzz "); break;
			case 4 : fprintf(stderr,"Txx "); break;
			case 5 : fprintf(stderr,"S "); break;
			case 6 : fprintf(stderr,"Fx "); break;
			case 7 : fprintf(stderr,"Fz "); break;
		}
		fprintf(stderr,"\n");
		if (wav->random) vmess("Wavelet has a random signature with fmax=%.2f", wav->fmax);
		if (nsrc>1) {
			vmess("*******************************************");
			vmess("*********** source array info *************");
			vmess("*******************************************");
			vmess("Areal source array is defined with %d sources.",nsrc);
			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(wav->nx*(wav->nt/(1024.0*1024.0))));
			if (src->plane) vmess("Computed p-value = %f.",p);
		}
		if (src->random) {
		vmess("Sources are placed at random locations withing the domain: ");
		vmess(" x[%.2f : %.2f]  z[%.2f : %.2f] ", xsrc1, xsrc2, zsrc1, zsrc2);
		vmess(" and time window  t[%.3f : %.3f]", tsrc1, tsrc2);
		}
	}

	/* define snapshots */

	if (!getparfloat("tsnap1", &tsnap1)) tsnap1=0.1;
	if (!getparfloat("tsnap2", &tsnap2)) tsnap2=0.0;
	if (!getparfloat("dtsnap", &dtsnap)) dtsnap=0.1;
	if (!getparfloat("dxsnap", &dxsnap)) dxsnap=dx;
	if (!getparfloat("dzsnap", &dzsnap)) dzsnap=dz;
	if (!getparfloat("xsnap1", &xsnap1)) xsnap1=sub_x0;
	if (!getparfloat("xsnap2", &xsnap2)) xsnap2=xmax;
	if (!getparfloat("zsnap1", &zsnap1)) zsnap1=sub_z0;
	if (!getparfloat("zsnap2", &zsnap2)) zsnap2=zmax;
	if (!getparint("sna_vxvztime", &sna->vxvztime)) sna->vxvztime=0;

	if (!getparint("sna_type_vz", &sna->type.vz)) sna->type.vz=1;
	if (!getparint("sna_type_vx", &sna->type.vx)) sna->type.vx=0;
	if (mod->ischeme>2) {
		sna->type.p=0;
		if (!getparint("sna_type_txx", &sna->type.txx)) sna->type.txx=0;
		if (!getparint("sna_type_tzz", &sna->type.tzz)) sna->type.tzz=0;
		if (!getparint("sna_type_txz", &sna->type.txz)) sna->type.txz=0;
		if (!getparint("sna_type_pp", &sna->type.pp)) sna->type.pp=0;
		if (!getparint("sna_type_ss", &sna->type.ss)) sna->type.ss=0;
	}
	else {
		if (!getparint("sna_type_p", &sna->type.p)) sna->type.p=1;
		sna->type.txx=0;
		sna->type.tzz=0;
		sna->type.txz=0;
		sna->type.pp=0;
		sna->type.ss=0;
	}

	sna->nsnap = 0;
	if (tsnap2 >= tsnap1) {
		sna_nrsna   = 1+NINT((tsnap2-tsnap1)/dtsnap);
		sna->skipdt = MAX(1,NINT(dtsnap/dt));
		sna->skipdx = MAX(1,NINT(dxsnap/dx));
		sna->skipdz = MAX(1,NINT(dzsnap/dz));
		sna->delay  = NINT(tsnap1/dt);
		isnapmax1   = (sna_nrsna-1)*sna->skipdt;
		isnapmax2   = floor( (mod->nt-(sna->delay + 1))/sna->skipdt) * sna->skipdt;
		isnapmax    = (sna->delay + 1) + MIN(isnapmax1,isnapmax2);
		sna->nsnap  = floor((isnapmax-(sna->delay + 1))/sna->skipdt) + 1;

		sna->x1=NINT((MIN(MAX(sub_x0,xsnap1),xmax)-sub_x0)/dx);
		sna->x2=NINT((MIN(MAX(sub_x0,xsnap2),xmax)-sub_x0)/dx);
		sna->z1=NINT((MIN(MAX(sub_z0,zsnap1),zmax)-sub_z0)/dz);
		sna->z2=NINT((MIN(MAX(sub_z0,zsnap2),zmax)-sub_z0)/dz);
		dxsnap=dx*sna->skipdx;
		dzsnap=dz*sna->skipdz;
		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);

		if (verbose) {
			vmess("*******************************************");
			vmess("************* snap shot info **************");
			vmess("*******************************************");
			vmess("tsnap1  = %f tsnap2  = %f ", tsnap1, tsnap2);
			vmess("dtsnap  = %f Nsnap   = %d ", dtsnap, sna->nsnap);
			vmess("nzsnap  = %d nxsnap  = %d ", sna->nz, sna->nx);
			vmess("dzsnap  = %f dxsnap  = %f ", dzsnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			if (sna->vxvztime) vmess("vx/vz snapshot time  : t+0.5*dt ");
			else vmess("vx/vz snapshot time  : t-0.5*dt ");
			fprintf(stderr,"    %s: Snapshot types        : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Vz ");
			if (sna->type.vx) fprintf(stderr,"Vx ");
			if (sna->type.p) fprintf(stderr,"p ");
			if (mod->ischeme>2) {
				if (sna->type.txx) fprintf(stderr,"Txx ");
				if (sna->type.tzz) fprintf(stderr,"Tzz ");
				if (sna->type.txz) fprintf(stderr,"Txz ");
				if (sna->type.pp) fprintf(stderr,"P ");
				if (sna->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
		}
	}
	else {
		sna->nsnap = 0;
		if (verbose) vmess("*************** no snapshots **************");
	}

	/* define receivers */

	if (!getparint("largeSUfile",&largeSUfile)) largeSUfile=0;
	if (!getparint("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparfloat("dtrcv",&dtrcv)) dtrcv=0.004;
	rec->skipdt=NINT(dtrcv/dt);
	dtrcv = mod->dt*rec->skipdt;
	if (!getparint("rec_ntsam",&rec->nt)) rec->nt=NINT(mod->tmod/dtrcv)+1;
	if (!getparint("rec_delay",&rec->delay)) rec->delay=0;
	if (!getparint("rec_int_vx",&rec->int_vx)) rec->int_vx=0;
	if (!getparint("rec_int_vz",&rec->int_vz)) rec->int_vz=0;
	if (!getparint("max_nrec",&max_nrec)) max_nrec=10000;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
	if(!getparfloat("dxrcv", &dxrcv)) dxrcv=dx;
	if(!getparfloat("dzrcv", &dzrcv)) dzrcv=0;
	rec->nt=MIN(rec->nt, NINT(mod->tmod/dtrcv)+1);

	rec->x = (int *)malloc(max_nrec*sizeof(int));
	rec->z = (int *)malloc(max_nrec*sizeof(int));
	recvPar(rec->x, rec->z, &rec->n, sub_x0, sub_z0, dx, dz, nx, nz);

	if (!getparint("rec_type_vz", &rec->type.vz)) rec->type.vz=1;
	if (!getparint("rec_type_vx", &rec->type.vx)) rec->type.vx=0;
	if (mod->ischeme>2) {
		rec->type.p=0;
		if (!getparint("rec_type_txx", &rec->type.txx)) rec->type.txx=0;
		if (!getparint("rec_type_tzz", &rec->type.tzz)) rec->type.tzz=0;
		if (!getparint("rec_type_txz", &rec->type.txz)) rec->type.txz=0;
		if (!getparint("rec_type_pp", &rec->type.pp)) rec->type.pp=0;
		if (!getparint("rec_type_ss", &rec->type.ss)) rec->type.ss=0;
	}
	else {
		if (!getparint("rec_type_p", &rec->type.p)) rec->type.p=1;
		rec->type.txx=0;
		rec->type.tzz=0;
		rec->type.txz=0;
		rec->type.pp=0;
		rec->type.ss=0;
	}

	if (!rec->type.vx) rec->int_vx=0;
	if (!rec->type.vz) rec->int_vz=0;
//	rec_movx = NINT(dxspread/dx);
//	rec_movz = NINT(dzspread/dz);

	if (verbose) {
		if (rec->n) {
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %d nrcv    = %d ", rec->nt, rec->n);
			vmess("dtrcv   = %f              ", dtrcv );
			vmess("dzrcv   = %f dxrcv   = %f ", dzrcv, dxrcv);
			vmess("receivers at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->z[0]*dz+sub_z0, rec->z[rec->n-1]*dz+sub_z0);
			vmess("xmin    = %f xmax    = %f ", rec->x[0]*dx+sub_x0, rec->x[rec->n-1]*dx+sub_x0);
			vmess("receivers at gridpoints: ");
			vmess("izmin   = %d izmax   = %d ", rec->z[0], rec->z[rec->n-1]);
			vmess("ixmin   = %d ixmax   = %d ", rec->x[0], rec->x[rec->n-1]);
			if (rec->type.vx) {
				fprintf(stderr,"    %s: Receiver interpolation for Vx:",xargv[0]);
				if(rec->int_vx==0) fprintf(stderr,"vx->vx\n");
				if(rec->int_vx==1) fprintf(stderr,"vx->vz\n");
				if(rec->int_vx==2) fprintf(stderr,"vx->txx/tzz\n");
			}
			if (rec->type.vz) {
				fprintf(stderr,"    %s: Receiver interpolation for Vz:",xargv[0]);
				if(rec->int_vz==0) fprintf(stderr,"vz->vz\n");
				if(rec->int_vz==1) fprintf(stderr,"vz->vx\n");
				if(rec->int_vz==2) fprintf(stderr,"vz->txx/tzz\n");
			}
			fprintf(stderr,"    %s: Receiver types        : ",xargv[0]);
			if (rec->type.vz) fprintf(stderr,"Vz ");
			if (rec->type.vx) fprintf(stderr,"Vx ");
			if (rec->type.p) fprintf(stderr,"p ");
			if (mod->ischeme>2) {
				if (rec->type.txx) fprintf(stderr,"Txx ");
				if (rec->type.tzz) fprintf(stderr,"Tzz ");
				if (rec->type.txz) fprintf(stderr,"Txz ");
				if (rec->type.pp) fprintf(stderr,"P ");
				if (rec->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
			if (mod->nt*mod->dt > rec->nt*dtrcv) {
				int nfiles = ceil((mod->nt*mod->dt)/(rec->nt*dtrcv));
				int lastn = floor((mod->nt)%(rec->nt*rec->skipdt)/rec->skipdt)+1;
				vmess("Receiver recordings will be written to %d files",nfiles);
				vmess("Last file will contain %d samples",lastn);
				
			}
		}
		else {
		 	vmess("*************** no receivers **************");
		}
	}

	return 0;
}

