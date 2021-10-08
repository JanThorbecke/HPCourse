#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

void readtable2D(complex *oper, float k, int hoplx, int hoply, int mode);

void conv2D(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode)
{
	int 	ix, iy, j, hoplx, hoply, ix2, iy2, index3, index4;
	int 	lenx, leny, tmpsize, opersize, hoplx2, hoply2;
	float	c=0;
	register float dumr, dumi;
	complex *copy, *tmp1, *tmp2, *tmp3, *tmp4, *hopx, *dum;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hoplx2 = hoplx-1;
	hoply2 = hoply-1;
	lenx = nx+2*hoplx2;
	leny = ny+2*hoply2;
	tmpsize = leny*hoplx+nx;
	opersize = hoplx*hoply;

	copy = (complex *)calloc(lenx*leny, sizeof(complex));
	dum  = (complex *)calloc(lenx, sizeof(complex));
	tmp1 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp2 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp3 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp4 = (complex *)calloc(tmpsize, sizeof(complex));
	hopx = (complex *)calloc(opersize, sizeof(complex));

/* Copy data into another array with zero's added to the edges */

	for (iy = 0; iy < ny; iy++) {
		memcpy(&copy[(hoply2+iy)*lenx+hoplx2], &data[iy*nx],
			nx*sizeof(complex));
	}

	memset( (float *)&data[0].r, 0, 2*nx*ny*sizeof(float) );

/* fill temporary arrays */

	for (iy = 0; iy < leny; iy++) {
		#pragma ivdep
		#pragma vector always
		for (ix = 0; ix < hoplx; ix++) {
			tmp1[iy*hoplx+ix]    = copy[iy*lenx+hoplx2+ix];
			tmp2[iy*hoplx+ix+nx] = copy[iy*lenx+hoplx2-ix];
			tmp3[iy*hoplx+ix].r  = tmp1[iy*hoplx+ix].r + 
								   tmp2[iy*hoplx+ix+nx].r;
			tmp3[iy*hoplx+ix].i  = tmp1[iy*hoplx+ix].i + 
								   tmp2[iy*hoplx+ix+nx].i;
		}
	}

	for (iy = 0; iy < leny; iy++) {
		memcpy(&tmp4[iy*hoplx], &tmp3[(leny-iy-1)*hoplx],
			hoplx*sizeof(complex));
	}

/* The 2D-Convolution */

	for (ix = 0; ix < nx; ix++) {
		for (iy = 0; iy < ny; iy++) {

/* if velocity changes calculate new operator */

            if (velmod[iy*nx+ix] != c) {
                c = velmod[iy*nx+ix];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
				hopx[0].r *=0.25;
				hopx[0].i *=0.25;
				for (ix2 = 1; ix2 < hoplx; ix2++) {
					hopx[ix2].r *=0.5;
					hopx[ix2].i *=0.5;
				}
				for (iy2 = 1; iy2 < hoply; iy2++) {
					hopx[iy2*hoplx].r *= 0.5;
					hopx[iy2*hoplx].i *= 0.5;
				}
            }

			index3 = (hoply2+iy)*hoplx;
			index4 = (leny-hoply2-iy-1)*hoplx;
/*
			#pragma vector aligned
			#pragma ivdep
			for (j = 0; j < opersize; j++) {
//				dumr = tmp3[index3+j].r + tmp4[index4+j].r;
//				dumi = tmp3[index3+j].i + tmp4[index4+j].i;

				dum[j].r = tmp3[index3+j].r + tmp4[index4+j].r;
				dum[j].i = tmp3[index3+j].i + tmp4[index4+j].i;
			}
*/

			#pragma ivdep
			#pragma vector always
			for (j = 0; j < opersize; j++) {
				dumr = tmp3[index3+j].r + tmp4[index4+j].r;
				dumi = tmp3[index3+j].i + tmp4[index4+j].i;
				data[iy*nx+ix].r += dumr*hopx[j].r;
				data[iy*nx+ix].r += dumi*hopx[j].i;
				data[iy*nx+ix].i += dumi*hopx[j].r;
				data[iy*nx+ix].i -= dumr*hopx[j].i;
			}

			/*
			for (iy2 = 0; iy2 < hoply; iy2++) {
				#pragma ivdep
				for (ix2 = 0; ix2 < hoplx; ix2++) {
					dumr = tmp3[index3+iy2*hoplx+ix2].r + tmp3[(hoply2-iy2+iy)*hoplx+ix2].r;
					dumi = tmp3[index3+iy2*hoplx+ix2].i + tmp3[(hoply2-iy2+iy)*hoplx+ix2].i;

					data[iy*nx+ix].r += dumr*hopx[j].r;
					data[iy*nx+ix].r += dumi*hopx[j].i;
					data[iy*nx+ix].i += dumi*hopx[j].r;
					data[iy*nx+ix].i -= dumr*hopx[j].i;
				}
			}
			*/

		}

		for (iy2 = 0; iy2 < leny; iy2++) {
			tmp1[iy2*hoplx+hoplx+ix] = copy[iy2*lenx+2*hoplx2+ix+1];
			tmp2[iy2*hoplx+nx-ix-1]  = copy[iy2*lenx+hoplx+ix];
		}

		for (ix2 = 0; ix2 < leny*hoplx; ix2++) {
			tmp3[ix2].r = tmp1[ix2+ix+1].r + tmp2[ix2+nx-ix-1].r;
			tmp3[ix2].i = tmp1[ix2+ix+1].i + tmp2[ix2+nx-ix-1].i;
		}

		for (iy2 = 0; iy2 < leny; iy2++) {
			memcpy(&tmp4[iy2*hoplx], &tmp3[(leny-iy2-1)*hoplx],
				hoplx*sizeof(complex));
		}

	}

	free(copy);
	free(tmp1);
	free(tmp2);
	free(tmp3);
	free(tmp4);
	free(hopx);

	return;
}


void conv2D_q8(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo;
	float	c=0;
	complex *term1, *oct, dum, *hopx;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	oct   = (complex *)malloc(hoplx*hoply*sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		memcpy(&term1[(iy+hy2)*nxo+hx2], &data[iy*nx], nx*sizeof(complex));
	}

	for (iy = hy2; iy < nyo-hy2; iy++) {
		for (ix = hx2; ix < nxo-hx2; ix++) {

/* calculating the sum at the x-axis and the diagonal x=y */
/* First make use of symmetry in x and y axis */

			oct[0] = term1[iy*nxo+ix];
			#pragma ivdep
			#pragma vector always
			for (i = 1; i < hoply; i++) {
				oct[i*hoplx].r = (term1[(iy-i)*nxo+ix].r + 
								  term1[iy*nxo+ix-i].r + 
								  term1[(iy+i)*nxo+ix].r + 
								  term1[iy*nxo+ix+i].r);
				oct[i*hoplx].i = (term1[(iy-i)*nxo+ix].i + 
								  term1[iy*nxo+ix-i].i + 
								  term1[(iy+i)*nxo+ix].i + 
								  term1[iy*nxo+ix+i].i);
				oct[i*hoplx+i].r = (term1[(iy-i)*nxo+ix-i].r + 
								  term1[(iy+i)*nxo+ix-i].r + 
								  term1[(iy-i)*nxo+ix+i].r + 
								  term1[(iy+i)*nxo+ix+i].r);
				oct[i*hoplx+i].i = (term1[(iy-i)*nxo+ix-i].i + 
								  term1[(iy+i)*nxo+ix-i].i + 
								  term1[(iy-i)*nxo+ix+i].i + 
								  term1[(iy+i)*nxo+ix+i].i);
			}

/* Second make use of the diagonal symmetry (only if dx == dy) */

			for (i = 2; i < hoply; i++) {
				#pragma ivdep
				#pragma vector always
				for (j = 1; j < i; j++) {
					oct[i*hoplx+j].r = 
						term1[(iy+i)*nxo+ix+j].r + term1[(iy+i)*nxo+ix-j].r + 
						term1[(iy-i)*nxo+ix+j].r + term1[(iy-i)*nxo+ix-j].r +
						term1[(iy+j)*nxo+ix+i].r + term1[(iy+j)*nxo+ix-i].r + 
						term1[(iy-j)*nxo+ix-i].r + term1[(iy-j)*nxo+ix+i].r;
					oct[i*hoplx+j].i = 
						term1[(iy+i)*nxo+ix+j].i + term1[(iy+i)*nxo+ix-j].i + 
						term1[(iy-i)*nxo+ix+j].i + term1[(iy-i)*nxo+ix-j].i +
						term1[(iy+j)*nxo+ix+i].i + term1[(iy+j)*nxo+ix-i].i + 
						term1[(iy-j)*nxo+ix-i].i + term1[(iy-j)*nxo+ix+i].i;
				}
			}

/* if velocity changes calculate new operator */

            if (velmod[(iy-hy2)*nx+ix-hx2] != c) {
                c = velmod[(iy-hy2)*nx+ix-hx2];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
            }

/* convolution with the operator */

			dum.r = 0.0;
			dum.i = 0.0;
			for (i = 0; i < hoply; i++) {
				for (j = 0; j <= i; j++) {
					dum.r += oct[i*hoplx+j].r*hopx[i*hoplx+j].r;
					dum.r += oct[i*hoplx+j].i*hopx[i*hoplx+j].i;
					dum.i += oct[i*hoplx+j].i*hopx[i*hoplx+j].r;
					dum.i -= oct[i*hoplx+j].r*hopx[i*hoplx+j].i;
				}
			}
			data[(iy-hy2)*nx+ix-hx2] = dum;

			/*
			if (dum.r != 0) fprintf(stderr,"ix = %d iy =%d dum = %f, %f\n", ix-hx2, iy-hy2, dum.r, dum.i);
			*/
		}
	}

	free(term1);
	free(oct);
	free(hopx);

	return;
}



void conv2D_q4(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo;
	float	c=0;
	complex *term1, *oct, dum, *hopx;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	oct   = (complex *)malloc(hoplx*hoply*sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		memcpy(&term1[(iy+hy2)*nxo+hx2], &data[iy*nx], nx*sizeof(complex));
	}

	for (iy = hy2; iy < nyo-hy2; iy++) {
		for (ix = hx2; ix < nxo-hx2; ix++) {

/* First make use of symmetry in x and y axis */

			oct[0] = term1[iy*nxo+ix];
			for (i = 1; i < hoply; i++) {
				#pragma ivdep
				#pragma vector always
				for (j = 1; j < hoplx; j++) {
					oct[i*hoplx+j].r = (term1[(iy-i)*nxo+ix-j].r + 
								  term1[(iy+i)*nxo+ix-j].r + 
								  term1[(iy-i)*nxo+ix+j].r + 
								  term1[(iy+i)*nxo+ix+j].r);
					oct[i*hoplx+j].i = (term1[(iy-i)*nxo+ix-j].i + 
								  term1[(iy+i)*nxo+ix-j].i + 
								  term1[(iy-i)*nxo+ix+j].i + 
								  term1[(iy+i)*nxo+ix+j].i);
				}
			}
			for (j = 1; j < hoplx; j++) {
				oct[j].r = (term1[iy*nxo+ix-j].r + 
							term1[iy*nxo+ix+j].r );
				oct[j].i = (term1[iy*nxo+ix-j].i + 
							term1[iy*nxo+ix+j].i );
			}
			for (i = 1; i < hoply; i++) {
				oct[i*hoplx].r = (term1[(iy-i)*nxo+ix].r + 
								term1[(iy+i)*nxo+ix].r );
				oct[i*hoplx].i = (term1[(iy-i)*nxo+ix].i + 
								term1[(iy+i)*nxo+ix].i );
			}

/* if velocity changes calculate new operator */

            if (velmod[(iy-hy2)*nx+ix-hx2] != c) {
                c = velmod[(iy-hy2)*nx+ix-hx2];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
            }

/* convolution with the operator */

			dum.r = 0.0;
			dum.i = 0.0;
			for (i = 0; i < hoply; i++) {
				#pragma ivdep
				#pragma vector always
				for (j = 0; j < hoplx; j++) {
					dum.r += oct[i*hoplx+j].r*hopx[i*hoplx+j].r;
					dum.r += oct[i*hoplx+j].i*hopx[i*hoplx+j].i;
					dum.i += oct[i*hoplx+j].i*hopx[i*hoplx+j].r;
					dum.i -= oct[i*hoplx+j].r*hopx[i*hoplx+j].i;
				}
			}
			data[(iy-hy2)*nx+ix-hx2] = dum;
		}
	}

	free(term1);
	free(oct);
	free(hopx);

	return;
}


void conv2D_q1(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo, ix2, iy2;
	int     starty, endy, startx, endx, k, l;
	float	c=0;
	complex *convr, *opx, dum, *hopx;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	convr = (complex *)malloc(nxo*nyo*sizeof(complex));
	opx   = (complex *)calloc(oplx*oply,sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		starty = MAX(iy-hoply+1, 0);
		endy   = MIN(iy+hoply, ny);

		for (ix = 0; ix < nx; ix++) {
			startx = MAX(ix-hoplx+1, 0);
			endx   = MIN(ix+hoplx, nx);
            if (velmod[iy*nx+ix] != c) {
                c = velmod[iy*nx+ix];
				readtable2D(hopx, om/c, hoplx, hoply, mode);

				for (iy2 = oply/2; iy2 < oply; iy2++) {
					for (ix2 = oplx/2; ix2 < oplx; ix2++) {
						opx[iy2*oplx+ix2] = hopx[(iy2-hoply+1)*hoplx+ix2-hoplx+1];
					}
					for (ix2 = 0; ix2 < oplx/2; ix2++) {
						opx[iy2*oplx+ix2] = hopx[(iy2-hoply+1)*hoplx+hoplx-ix2-1];
					}
				}
				for (iy2 = 0; iy2 < oply/2; iy2++) {
					for (ix2 = 0; ix2 < oplx; ix2++) {
						opx[iy2*oplx+ix2] = opx[(oply-1-iy2)*oplx+ix2];
					}
				}

            }

			dum.r = dum.i = 0.0;
			k = MAX(hoply-1-iy, 0);
			for (i = starty; i < endy; i++) {
				l = MAX(hoplx-1-ix, 0);
				#pragma ivdep
				#pragma vector always
				for (j = startx; j < endx; j++) {
					dum.r += data[i*nx+j].r*opx[k*oplx+l].r;
					dum.r -= data[i*nx+j].i*opx[k*oplx+l].i;
					dum.i += data[i*nx+j].r*opx[k*oplx+l].i;
					dum.i += data[i*nx+j].i*opx[k*oplx+l].r;
					l++;
				}
				k++;
			}

			convr[iy*nx+ix] = dum;
			/*
			if (dum.r != 0) fprintf(stderr,"ix = %d iy =%d dum = %f, %f\n", ix, iy, dum.r, dum.i);
			*/
		}
	}

	memcpy(data, convr, nx*ny*sizeof(float));

	free(convr);
	free(opx);
	free(hopx);

	return;
}

