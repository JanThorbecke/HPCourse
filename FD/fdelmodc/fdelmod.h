#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct _compType { /* Receiver Type */
	int vz;
	int vx;
	int p;
	int txx;
	int tzz;
	int txz;
	int pp;
	int ss;
} compType;

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	char *file_beam;
	compType type;
	int n;
	int nt;
	int delay;
	int skipdt;
	int *z;
	int *x;
	int int_vx;
	int int_vz;
	int sinkdepth;
} recPar;

typedef struct _snapshotPar { /* Snapshot Parameters */
	char *file_snap;
	compType type;
	int nsnap;
	int delay;
	int skipdt;
	int skipdz;
	int skipdx;
	int nz;
	int nx;
	int z1;
	int z2;
	int x1;
	int x2;
	int vxvztime;
} snaPar;

typedef struct _modelPar { /* Model Parameters */
	int iorder;
	int ischeme;
	int grid_dir;
	char *file_cp;
	char *file_ro;
	char *file_cs;
	char *file_qp;
	char *file_qs;
	float dz;
	float dx;
	float dt;
	float tmod;
	int nt;
	float z0;
	float x0;
	int nz;
	int nx;
	int naz;
	int nax;
	int ioz;
	int iox;
	float Qp;
	float Qs;
	float fw;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src;
	int nt;
	int nx;
	float dt;
	float fmax;
	int random;
} wavPar;

typedef struct _sourcePar { /* Source Array Parameters */
	int n;
	int type;
	int orient;
	int *z;
	int *x;
	int random;
	float *tbeg;
	float *tend;
	int multiwav;
	int plane;
	float angle;
	float velo;
	float amplitude;
	int distribution;
	int window;
} srcPar;

typedef struct _shotPar { /* Shot Parameters */
	int n;
	int *z;
	int *x;
} shotPar;

typedef struct _boundPar { /* Boundary Parameters */
	int free[4];
	int tap[4];
	int rig[4];
	float *tapz;
	float *tapx;
	int cfree;
	int *surface;
} bndPar;

