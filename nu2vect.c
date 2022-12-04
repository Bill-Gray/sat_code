/* Code to read in the Spektr-RG state vectors from .nu files at

ftp://ftp.kiam1.rssi.ru/pub/gps/spectr-rg/nu/

   and output into a form ingestible by Find_Orb.  (Note that these
haven't been updated for a while;  we've been tracking Spektr-RG
solely through observed astrometry.)  Input lines give

0.000    (always zero,  means J2000)
1        (loop number;  ignore)
2.021030100000E+7    ( = 20210301 = 2021 Mar 01)
 1.649335000000E+05  ( = 16:49:33.5 Moscow time (three hours ahead of UTC))
-1.190093387974E+03  (x-coord, geocentric, equatorial J2000,  in thousands of km)
 9.733125144680E+02  (y-coord,  same system)
 4.956578219661E+02  (z-coord,  same system)
-7.639675611173E-02  (x-velocity, in km/s)
-6.678070845750E-02  (y-velocity)
-2.372125806637E-01  (z-velocity)
 0.000000000000E+00  (ballistic coefficient)
 1.246012245177E-05  (solar radiation coefficient, unitless)

   Note that the file name gives the UTC,  and the time given within
the file is three hours ahead of that.  The positions match those
determined by optical astrometry to within a few kilometers.

   The above solar radiation coefficient means that the sun's gravity
is counteracted by a force 1.24601e-5 times as great.  It corresponds
approximately to an area/mass ratio of 0.015 m^2/kg,  which is
quite close to what we've been getting from the optical astrometry
orbit solution.

   Compile with

gcc -Wextra -Wall -O3 -pedantic nu2vect.c -I../include -L../lib -o nu2vect -llunar

   Python code to convert .nu files to STK ephemeris files is at

https://github.com/Satsir/STK/blob/main/nuToEph.py     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "stringex.h"

int main( const int argc, const char **argv)
{
   FILE *ifile = (argc == 2 ? fopen( argv[1], "rb") : NULL);
   char buff[90];
   char command[200];
   int i;

   if( !ifile)
      fprintf( stderr, "See comments at start of 'nu2vect.c' for usage\n");
   assert( ifile);
   for( i = 0; i < 10; i++)
      if( fgets( buff, sizeof( buff), ifile))
         {
         const double ival = atof( buff);

         switch( i)
            {
            case 2:
               {
               const int day = (int)ival;

               snprintf_err( buff, sizeof( buff), "%04d-%02d-%02d",
                        day / 10000, (day / 100) % 100, day % 100);
               printf( "Date : %s\n", buff);
               snprintf_err( command, sizeof( command),
                       "find_orb \"-oSpektr-RG = 2019-040A = NORAD 44432\" -v%sT", buff);
               }
               break;
            case 3:
               {
               const int millisec = (int)( ival * 1000.);

               snprintf_err( buff, sizeof( buff), "%02d:%02d:%02d.%03d",
                        millisec / 10000000,
                        (millisec / 100000) % 100,
                        (millisec / 1000) % 100, millisec % 100);
               printf( "Time : %s\n", buff);
               strlcat_error( command, buff);
               strlcat_error( command, "-3h");   /* correction for Moscow time */
               }
               break;
            case 4: case 5: case 6:
               printf( "  Posn %f km\n", ival * 1000.);
               snprintf_append( command, sizeof( command), ",%.2f", ival * 1000.);
               break;
            case 7: case 8: case 9:
               printf( "  Vel  %f km/s\n", ival);
               snprintf_append( command, sizeof( command), ",%.7f", ival);
               break;
            break;
            }
         }
   fclose( ifile);
   printf( "%s,eq,geo,km,s\n", command);
   return( 0);
}
