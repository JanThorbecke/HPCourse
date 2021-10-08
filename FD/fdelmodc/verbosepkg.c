#include <stdio.h>
#include <stdarg.h>
#include "par.h"
#include <errno.h>
#include <string.h>
#ifdef _CRAYMPP
#include <intrinsics.h>
#endif

void verr(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nverr: fflush failed on stdout");
	}
	fprintf(stderr, "    Error in %s: ", xargv[0]);
#ifdef _CRAYMPP
        fprintf(stderr, "PE %d: ", _my_pe());
#elif defined(SGI)
        fprintf(stderr, "PE %d: ", mp_my_threadnum());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");

	exit(EXIT_FAILURE);
}

void vwarn(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nvwarn: fflush failed on stdout");
	}
	fprintf(stderr, "    Warning in %s: ", xargv[0]);
#ifdef _CRAYMPP
        fprintf(stderr, "PE %d: ", _my_pe());
#elif defined(SGI)
        fprintf(stderr, "PE %d: ", mp_my_threadnum());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}

extern int errno;

void vsyserr(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nsyserr: fflush failed on stdout");
	}
	fprintf(stderr, "    %s: ", xargv[0]);
#ifdef _CRAYMPP
        fprintf(stderr, "PE %d: ", _my_pe());
#elif defined(SGI)
        fprintf(stderr, "PE %d: ", mp_my_threadnum());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "error number %d: %s\n",errno,strerror(errno));
	exit(EXIT_FAILURE);
}

void vmess(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nvmess: fflush failed on stdout");
	}
	fprintf(stderr, "    %s: ", xargv[0]);
#ifdef _CRAYMPP
        fprintf(stderr, "PE %d: ", _my_pe());
#elif defined(SGI)
        fprintf(stderr, "PE %d: ", mp_my_threadnum());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}
