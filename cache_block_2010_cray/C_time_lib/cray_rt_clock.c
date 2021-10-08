#
#ifdef CRAY_X2
  #include <intrinsics.h>
#endif
#include <stdlib.h>

#define FUNCTION_ID Cray_mpi_wtime
#define FUNCTION_ID_F Cray_mpi_wtime_

#ifdef CRAY_USER_WTIME
#define FUNCTION_ID CRAY_USER_WTIME
#endif
#ifdef CRAY_USER_WTIME_F
#define FUNCTION_ID_F CRAY_USER_WTIME_F
#endif

unsigned long long irtc_rate();

#ifndef CRAY_X2
unsigned long long _rtc (void);
#endif

#ifdef CRAY_NO_INIT_CALL
static int first_call=0;
static int64_t clk_rate;
#endif

static double   tckclk, clktck, scale_factor;
static int64_t  iclktck,init_tick;

double FUNCTION_ID(void);
double FUNCTION_ID_F(void);

int64_t Cray_timer_init(void)
{
   iclktck = irtc_rate();
   clktck = (double) iclktck;
   scale_factor = 1.0/1000000.0;
   tckclk = 1.0/(clktck*scale_factor);
   init_tick = _rtc();
   return iclktck;
}

double FUNCTION_ID(void)
{
   double tmpt,tmpf;
#ifdef CRAY_NO_INIT_CALL
   if (first_call==0)  {
      clk_rate=Cray_timer_init();
      first_call=1;
#ifdef CRAY_LIST_CLOCK_RATE
      printf("*** clock rate %f ***\n",(double) clk_rate);
#endif
   }
#endif
   tmpf = ((double) (_rtc()-init_tick) )*tckclk*scale_factor;
   return tmpf;
}
double FUNCTION_ID_F(void)
{
   double tmpt;
#ifdef CRAY_NO_INIT_CALL
   if (first_call==0)  {
      clk_rate=Cray_timer_init();
      first_call=1;
   }
#endif
   tmpt = ((double) (_rtc()-init_tick) )*tckclk*scale_factor;
   return tmpt;
}
