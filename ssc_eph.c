/* Code to convert ephems from SSCWeb into the format 'eph2tle'
uses to generate TLEs.  The usefulness here is that you can
get TLEs for the Magnetospheric Multiscale (MMS) satellites that
are good enough to ID the individual satellites,  even though
the four of them are in a tight cluster.  This does fail if the
satellites maneuver during the course of a day;  the TLEs are
then unable to fit the (non-gravitational) motion.  Compile with

gcc -I../include -Wextra -Wall -O3 -pedantic ssc_eph.c -o ssc_eph ../lib/liblunar.a -lm

  Go to the SSCWeb Locator page :

https://sscweb.gsfc.nasa.gov/cgi-bin/Locator.cgi

   Select the 'Standard' interface.  Turn on MMS-1, 2, 3,  and 4.
Set the output frequency to 144 minutes (= 0.1 day).  Under
'Optional Settings',  select kilometers and yy/mm/dd output.
Select the desired time span.

   Under 'Output Options',  select GEI/J2000 XYZ.

   Generate the output and save it as /tmp/mms.txt.  The result
can be fed through 'eph2tle' (see the 'find_orb' repository) to
generate TLEs.       */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"

#define is_power_of_two( X)   (!((X) & ((X) - 1)))

static void dump_to_file( const int mms_no, const double t0, const double dt,
                     const int n_found, const double *xyz)
{
   FILE *ofile;
   char buff[200];
   int i, j;

   snprintf( buff, sizeof( buff), "mms%d.txt", mms_no);
   ofile = fopen( buff, "wb");
   assert( ofile);
   fprintf( ofile, "%f %f %d 0,149597870.700000,86400.000000,2000"
                     " (500) Geocentric: MMS-%d = 2015-011%c = NORAD %d\n",
                            t0, dt, n_found, mms_no, mms_no + '@', mms_no + 40481);
   for( i = 0; i < n_found; i++)
      {
      double vel[3];
      const double *tptr = xyz + i * 3;

      for( j = 0; j < 3; j++)
         vel[j] = (i == n_found - 1 ? tptr : tptr + 3)[j];
      for( j = 0; j < 3; j++)
         {
         vel[j] -= (i ? tptr - 3 : tptr)[j];
         vel[j] /= 86400. * dt;
         if( i && i != n_found - 1)
            vel[j] /= 2.;
         }
      fprintf( ofile, "%f %11.2f %11.2f %11.2f %11.6f %11.6f %11.6f\n",
                  t0 + (double)i * dt,
                  tptr[0], tptr[1], tptr[2], vel[0], vel[1], vel[2]);
      }
   fclose( ofile);
}

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( "/tmp/mms.txt", "rb");
   char buff[100];
   long header_offset = 0;
   int mms_no;

   assert( ifile);
   while( !header_offset && fgets( buff, sizeof( buff), ifile))
      if( !memcmp( buff, "yy/mm/dd", 8))
         header_offset = ftell( ifile);
   if( !header_offset)
      {
      fprintf( stderr, "Couldn't find yy/mm/dd header\n");
      return( -1);
      }
   for( mms_no = 1; mms_no <= 4; mms_no++)
      {
      int n_found = 0;
      double t0 = 0., step = 0., *xyz = NULL;
      const double td_minus_utc = 69.184 / 86400.;
      char *tptr;

      fseek( ifile, header_offset, SEEK_SET);
      buff[0] = '2';       /* assume 2nd millennium */
      buff[1] = '0';       /* assume 21st century */
      while( fgets( buff + 2, sizeof( buff) - 2, ifile))
         if( (tptr = strstr( buff, "mms")) != NULL && tptr[3] == '0' + mms_no)
            {
            int n_fields_read;

            tptr[-1] = '\0';
            if( !n_found)
               t0 = get_time_from_string( 0., buff, FULL_CTIME_YMD, NULL);
            if( n_found == 1)
               step = get_time_from_string( 0., buff, FULL_CTIME_YMD, NULL) - t0;
            n_found++;
            if( is_power_of_two( n_found))
               xyz = (double *)realloc( xyz, 2 * n_found * 3 * sizeof( double));
            n_fields_read = sscanf( tptr + 4, "%lf %lf %lf",
                        xyz + n_found * 3 - 3,
                        xyz + n_found * 3 - 2,
                        xyz + n_found * 3 - 1);
            assert( 3 == n_fields_read);
            }
      dump_to_file( mms_no, t0 + td_minus_utc, step, n_found, xyz);
      free( xyz);
      }
   fclose( ifile);
   return( 0);
}
