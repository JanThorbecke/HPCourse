#ifndef GENFFT_H
#define GENFFT_H

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

#if defined(HAVE_LIBCOMPLIB_SGIMATH)
#include <fft.h>
#define COMPLEX
#endif
#ifdef CWP_H
#define COMPLEX
#endif

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#define COMPLEX
#endif/* complex */

#ifndef PI
#define PI (3.141592653589793)
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#define MAX_NUMTHREADS 128
#define MAXNUMTHREAD 128

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

int optncc(int n);
int optncr(int n);

void cc1fft(complex *data, int n, int sign);
void ccmfft(complex *data, int n1, int n2, int ld1, int sign);
void cc2dfft(complex *data, int n1, int n2, int ld1, int sign);
 
void cc1_fft(complex *cdata, int n1, int sign);
void ccm_fft(complex *data, int n1, int n2, int ld1, int sign);
void cc2d_fft(complex *data, int n1, int ld1, int n2, int sign);

void rc1fft(float *rdata, complex *cdata, int n, int sign);
void rcmfft(float *rdata, complex *cdata, int n1, int n2, int ldr, int ldc, int sign);
void rc2dfft(float *rdata, complex *cdata, int nr, int nc, int ldr, int ldc, int sign);

void cr1fft(complex *cdata, float *rdata, int n, int sign);
void crmfft(complex *cdata, float *rdata, int n1, int n2, int ldc, int ldr, int sign);
void cr2dfft(complex *cdata, float *rdata, int nr, int nc, int ldc, int ldr, int sign);

void rc1_fft(float *rdata, complex *cdata, int n, int sign);
void rcm_fft(float *rdata, complex *cdata, int n1, int n2, int ldr, int ldc, int isign);
void cr1_fft(complex *cdata, float *rdata, int n, int sign);
void crm_fft(complex *cdata, float *rdata, int n1, int n2, int ldc, int ldr, int sign);

void xt2wx(float *rdata, complex *cdata, int nt, int nx, int ldr, int ldc);
void xt2wkx(float *rdata, complex *cdata, int nt, int nx, int ldr, int ldc, int xorig);
void wx2xt(complex *cdata, float *rdata, int nt, int nx, int ldc, int ldr);
void wkx2xt(complex *cdata, float *rdata, int nt, int nx, int ldc, int ldr, int xorig);

void dc1_fft(double *rdata, dcomplex *cdata, int n, int sign);
void cd1_fft(dcomplex *cdata, double *rdata, int n, int sign);

int factorized(int n);
double wallclock_time(void);

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
}
#endif

/*----------------------------------------------------------------
** MACRO FNAME for any fortran callable routine name.
**
**  This macro prepends, appends, or does not modify a name
**  passed as a macro parameter to it based on the FNAME_PRE_UNDERSCORE,
**  FNAME_POST_UNDERSCORE macros set for a specific system.
**
**---------------------------------------------------------------*/
#if defined(FNAME_PRE_UNDERSCORE) && defined(FNAME_POST_UNDERSCORE)
#ifdef __STDC__
#   define FNAME(x)     _##x##_
#else
#   define FNAME(x)     _/**/x/**/_
#endif
#endif
#if defined(FNAME_PRE_UNDERSCORE) && !defined(FNAME_POST_UNDERSCORE)
#ifdef __STDC__
#   define FNAME(x)     _##x
#else
#   define FNAME(x)     _/**/x
#endif
#endif
#if !defined(FNAME_PRE_UNDERSCORE) && defined(FNAME_POST_UNDERSCORE)
#ifdef __STDC__
#   define FNAME(x)     x##_
#else
#   define FNAME(x)     x/**/_
#endif
#endif
#if !defined(FNAME_PRE_UNDERSCORE) && !defined(FNAME_POST_UNDERSCORE)
#   define FNAME(x)     x
#endif

#endif /* GENFFT_H */
