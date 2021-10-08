/*
 * IA-32 (Pentium or better) CPU hardware cycle counter register.
 * Works with Intel C/C++ (icc) or GNU C (gcc).
 *
 * Compiles as a C or Fortran callable object file by default.
 *
 * Compile with -DBARRIER to force the processor to complete all
 * outstanding instructions before returning the real time counter.
 *
 * Compile with -DTEST to create a test executable.
 */

#include <stdio.h>

#ifdef CRAY_XT3
extern int __cpu_mhz;
#endif

/*
 * _rtc - Processor real time clock.
 * Interface to match the Cray C/C++ intrinsic.
 */
  unsigned long long _rtc (void)
{
  unsigned long high, low;

#ifdef BARRIER
  /* Instruction barrier - complete outstanding instructions */
  /* CPUID instruction is opcode 0x0f 0xa2 */
  unsigned long version, feature;
  asm volatile(".byte 0x0f,0xa2"
    : "=a" (version), "=d" (feature)
    : "a" (1)
    : "ebx", "ecx" );
  /*printf("Version %08lx, feature %08lx\n",version,feature);*/
#endif

  /* read the 64 bit process cycle counter into low/high */ 
  /* RDTSC instruction is opcode 0x0f 0x31 */
  asm volatile(".byte 0x0f,0x31" : "=a" (low), "=d" (high));

  return ((unsigned long long) high << 32) + low;
}


/*
 * irtc_rate - real-time clock number of ticks per second.
 * Interface to match the Cray C/C++ intrinsic.
 * Uses the dubious /proc/cpuinfo file.
 * Linux/IA-32 version.
 */
unsigned long long irtc_rate()
#ifndef CRAY_XT3
{
  FILE *fp = fopen("/proc/cpuinfo","r");
  char buf[1024];
  double r;
  unsigned long long ir = -1;

  if (fp == NULL) return 0;

  while (fgets(buf,1024,fp) != NULL) {
    if (!strncmp(buf,"cpu MHz",7)) {
      sscanf(buf,"cpu MHz\t\t: %lf",&r);
      ir = r * 1e6;
      break;
    }
  }

  fclose(fp);
  return ir;
}
#else
{
  return (((long long)__cpu_mhz)*1000000);
}
#endif

/*
 * irtc - Fortran-callable real time clock.
 * Interface to match the Cray Fortran intrinsic.
 *
 * integer*8 irtc
 */
unsigned long long irtc_ (void)
{
  return _rtc();
}


/*
 * rtc - Fortran-callable real time clock.
 * Interface to match the Cray Fortran intrinsic.
 *
 * real*8 rtc
 */
double rtc_ (void)
{
  return (double) _rtc();
}


/*
 * irtc_rate - Fortran-callable real-time clock cycles per second.
 * Interface to match the Cray Fortran intrinsic.
 *
 * integer*8 irtc_rate
 */
unsigned long long irtc_rate_ (void)
{
  return irtc_rate();
}


#ifdef TEST

#include <stdlib.h>

int main()
{
  const int n = 1000;
  unsigned long long c1, c2, c3;
  int i;
  double x = 3.1415926;
  double rtc_rate = (double) irtc_rate();
  double rtc_mhz = rtc_rate / 1000000.0;
  double resolution = 1e9 / rtc_rate;
  double t, mflops;
  double *a, *b, *c;

  a = (double *) malloc( n * sizeof(double) );
  b = (double *) malloc( n * sizeof(double) );
  c = (double *) malloc( n * sizeof(double) );
  if ((a == NULL) || (b == NULL) || (c == NULL)) {
    fprintf(stderr,"Out of memory\n");
    exit(1);
  }

  for (i=0;i < n;i++) {
    a[i] = drand48();
    b[i] = drand48();
    c[i] = drand48();
  }

  c1 = _rtc();

  for (i=0;i < n;i++) {
    c[i] = c[i] + a[i] * b[i];
  }

  c2 = _rtc();

  c3 = c2 - c1;

  for (i=0;i < n;i++) {
    x = x + c[i];
  }

  printf("Result %lf\n",x);

  printf("CPU clock rate %lf MHz\n",rtc_mhz);
  printf("Resolution     %lf nanoseconds\n",resolution);

  printf("Start   cycles %lld\n",c1);
  printf("End     cycles %lld\n",c2);
  printf("Elapsed cycles %lld\n",c3);

  t = c3/rtc_mhz;
  printf("Elapsed time   %lf microseconds\n",t);

  if (t > 0) {
    mflops = 2 * n / t;
    printf("Elapsed MFLOPS %lf\n",mflops);
  }

  return 0;
}

#endif
