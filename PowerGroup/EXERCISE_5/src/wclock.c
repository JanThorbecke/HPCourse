#include <sys/time.h>
#include <sys/types.h>
#include <stdlib.h>
 
/*
-----------------------------------------------------------------------

    This function returns the wall clock time with micro seconds
    accuracy.
    The data type of the returned value is "double".

    The function can be called from a FORTRAN module. The value
    returned by wall_clock_ should be of type Real (Kind = 8).


-----------------------------------------------------------------------
*/

double wclock()
{
        static  double  c0 = 1.0e-06; /* Conversion constant */

        struct  timeval tp;     /* Structures used by gettimeofday */
        struct  timezone tzp;

        double  wall_time;


                if ( gettimeofday(&tp,&tzp) == -1 )
                        perror("gettimeofday");

                wall_time = tp.tv_sec + c0 * tp.tv_usec;

                return(wall_time);
}

double wclock_()

{
        return wclock();
}
