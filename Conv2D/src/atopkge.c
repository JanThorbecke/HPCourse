/* Copyright (c) Colorado School of Mines, 1997.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
ATOPKGE - convert ascii to arithmetic and with error checking

 
eatoh		ascii to short
eatou		ascii to unsigned short
eatoi		ascii to int
eatop		ascii to unsigned
eatol		ascii to long
eatov		ascii to unsigned long
eatof		ascii to float
eatod		ascii to double

****************************************************************************
Function Prototypes:
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

****************************************************************************
Input:
s		string 

Returned:	type indicated
 
****************************************************************************
Notes:
Each of these routines acts like atoi, but has error checking:

This is a major revision of the tedious code required before
vendors implemented the ANSI C strtol, strtoul and strtod.

In addition to the size checks for each integer type, a
specific test on errno is required.  For example, INT_MAX
may (and probably does) equal LONG_MAX.  In this case,
if fed a number exceeding INT_MAX (and LONG_MAX), strtol
will make a quiet return with the wrong answer and it is up
to the user to check if errno == ERANGE.

Size limits are machine dependent and are read from the
ANSI C include files limits.h and float.h.

Bug Report: With NeXT c and Gnucc, when x > DBL_MAX (x <-DBL_MAX),
the return value from strtod was +Infinity (-Infinity), not HUGE_VAL
and more important, errno was not set to ERANGE.  To cope with this,
I put explicit size checks in eatod (which would not be needed if
errno were set as it should be in ANSI C.    jkc 01/29/94

On IBM RS6000, the return value from strtod was +-Inf on
overflow, but errno was set correctly.

****************************************************************************
References:
For old code:
Plum: Reliable Data Structures in C, p. 2-17.
Kernighan and Ritchie: The C Programming Language, p. 58.

CWP: Jack K. Cohen, Brian Sumner
 
For new code:
ANSI C routines with a little help from Jack

****************************************************************************
Author: Jack Cohen, Center for Wave Phenomena, 1994.
***************************************************************************/
/**************** end self doc ********************************/

#include "par.h"
#include <float.h>
#include <limits.h>
#include <stdarg.h>
#include <errno.h>

/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} */
short eatoh(char *s)
{
	long n = strtol(s, NULL, 10);
	
	if ( (n > SHRT_MAX) || (n < SHRT_MIN) || (errno == ERANGE) )
		err("%s: eatoh: overflow", __FILE__);

	return (short) n;
}


/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} */
unsigned short eatou(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if ( (n > USHRT_MAX) || (errno == ERANGE) )
		err("%s: eatou: overflow", __FILE__);

	return (unsigned short) n;
}


/* eatoi - convert string s to integer {INT_MIN:INT_MAX} */
int eatoi(char *s)
{
	long n = strtol(s, NULL, 10);

	if ( (n > INT_MAX) || (n < INT_MIN) || (errno == ERANGE) )
		err("%s: eatoi: overflow", __FILE__);

	return (int) n;
}


/* eatop - convert string s to unsigned integer {0:UINT_MAX} */
unsigned int eatop(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if ( (n > UINT_MAX) || (errno == ERANGE) )
		err("%s: eatop: overflow", __FILE__);

	return (unsigned int) n;
}


/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX} */
long eatol(char *s)
{
	long n = strtol(s, NULL, 10);

	if (errno == ERANGE)
		err("%s: eatol: overflow", __FILE__);

	return n;
}


/* eatov - convert string s to unsigned long {0:ULONG_MAX} */
unsigned long eatov(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if (errno == ERANGE)
		err("%s: eatov: overflow", __FILE__);

	return n;
}


/* eatof - convert string s to float {-FLT_MAX:FLT_MAX} */
float eatof(char *s)
{
	float x = strtod(s, NULL);

	if ( (x > FLT_MAX) || (x < -FLT_MAX) || (errno == ERANGE) )
		err("%s: eatof: overflow", __FILE__);

	return (float) x;
}


/* eatod - convert string s to double {-DBL_MAX:DBL_MAX} */
double eatod(char *s)
{
	double x = strtod(s, NULL);

	/* errno == ERANGE suffices if compiler sets errno on overflow */
	if ( (errno == ERANGE) || (x > DBL_MAX) || (x < -DBL_MAX) )
		err("%s: eatod: overflow", __FILE__);

	return x;
}


/**************************************************************************
ERRPKGE - routines for reporting errors

err	print warning on application program error and die
warn	print warning on application program error
syserr	print warning on application program error using errno and die

***************************************************************************
Function Prototypes:
void err (char *fmt, ...);
void warn (char *fmt, ...);
void syserr (char *fmt, ...);

***************************************************************************
Return: void

***************************************************************************
Notes:
fmt		a printf format string ("\n" not needed)
...		the variables referenced in the format string

Examples:
	err("Cannot divide %f by %f", x, y);
	warn("fmax = %f exceeds half nyquist= %f", fmax, 0.25/dt);
 
	if (NULL == (fp = fopen(xargv[1], "r")))
 		err("can't open %s", xargv[1]);
 	...
 	if (-1 == close(fd))
 		err("close failed");

***************************************************************************
References:
Kernighan and Pike, "The UNIX Programming Environment", page 207.
Also Rochkind, "Advanced UNIX Programming", page 13.

***************************************************************************
Authors:SEP: Jeff Thorson, Stew Levin	CWP: Shuki Ronen, Jack Cohen
**************************************************************************/


void err(char *fmt, ...)
{
	va_list args;

 
	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nerr: fflush failed on stdout");
	}
	fprintf(stderr, "\n%s: ", xargv[0]);
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}


void warn(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nwarn: fflush failed on stdout");
	}
	fprintf(stderr, "\n%s: ", xargv[0]);
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}


void syserr(char *fmt, ...)
{
    va_list args;

    if (EOF == fflush(stdout)) {
        fprintf(stderr, "\nsyserr: fflush failed on stdout");
    }
    fprintf(stderr, "\n%s: ", xargv[0]);
    va_start(args,fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, " (%s)\n", strerror(errno));
    exit(EXIT_FAILURE);
}

#ifdef TEST
main(int argc, char **argv)
{
	char s[BUFSIZ];
	short nh;
	unsigned short nu;
	int ni;
	unsigned int np;
	long nl;
	unsigned long nv;

	initargs(argc, argv);


	/* Test code for eatoh */
	if (SHRT_MAX == LONG_MAX) {
	    warn("Warning: eatoh not used on this machine.\n");
	} else {
	    warn("\n");
	}
	strcpy(s, "0");
	nh = eatoh(s);
	warn("eatoh(%s) = %hd\n", s, nh);

	strcpy(s, "32767");
	nh = eatoh(s);
	warn("eatoh(%s) = %hd\n", s, nh);

	strcpy(s, "-32768");
	nh = eatoh(s);
	warn("eatoh(%s) = %hd\n", s, nh);


	/* Test code for eatou */
	if (USHRT_MAX == ULONG_MAX) {
	    warn("Warning: eatou not used on this machine.\n");
	} else {
	    warn("\n");
	}
	strcpy(s, "0");
	nu = eatou(s);
	warn("eatou(%s) = %hu\n", s, nu);

	strcpy(s, "65535");
	nu = eatou(s);
	warn("eatou(%s) = %hu\n", s, nu);


	/* Test code for eatoi */
	if (INT_MAX == LONG_MAX) {
	    warn("Warning: eatoi not used on this machine.\n");
	} else {
	    warn("\n");
	}
	strcpy(s, "0");
	ni = eatoi(s);
	warn("eatoi(%s) = %d\n", s, ni);

	strcpy(s, "2147483647");
	ni = eatoi(s);
	warn("eatoi(%s) = %d\n", s, ni);


	strcpy(s, "-2147483648");
	ni = eatoi(s);
	warn("eatoi(%s) = %d\n", s, ni);


	/* Test code for eatop */
	if (INT_MAX == LONG_MAX) {
	    warn("Warning: eatop not used on this machine.\n");
	} else {
	    warn("\n");
	}
	strcpy(s, "0");
	np = eatop(s);
	warn("eatop(%s) = %lu\n", s, np);

	strcpy(s, "4294967295");
	np = eatop(s);
	warn("eatop(%s) = %lu\n", s, np);


	/* Test code for eatol */
	warn("\n");
	strcpy(s, "0");
	nl = eatol(s);
	warn("eatol(%s) = %ld\n", s, nl);

	strcpy(s, "2147483647");
	nl = eatol(s);
	warn("eatol(%s) = %ld\n", s, nl);

	strcpy(s, "-2147483648");
	nl = eatol(s);
	warn("eatol(%s) = %ld\n", s, nl);


	/* Test code for eatov */
	strcpy(s, "0");
	nv = eatov(s);
	warn("eatov(%s) = %lu\n", s, nv);

	strcpy(s, "4294967295");
	nv = eatov(s);
	warn("eatov(%s) = %lu\n", s, nv);

	warn("Now we feed in 4294967296, expecting fatal error exit\n");
	strcpy(s, "4294967296");
	nv = eatov(s);
	warn("eatov(%s) = %lu\n", s, nv);

	return EXIT_SUCCESS;
}
#endif
