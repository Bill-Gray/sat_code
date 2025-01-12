/* See LICENSE.   NOTE that this has been obsoleted by the 'add_off.c'
program in the 'lunar' repository (q.v.).  The only use this code would
have would be for a spacecraft getting astrometric data for which we
don't have Horizons data.  I don't expect that to happen.  With that
disclaimer :

   The Minor Planet Center accepts astrometry from spacecraft using
a modification of their usual 80-column "punched-card" format.  A
second line is used to tell you where the spacecraft was,  relative
to the geocenter.

https://minorplanetcenter.net/iau/info/SatelliteObs.html

   This code can reset those positions for Earth-orbiting spacecraft using
TLEs ('two-line elements';  https://www.projectpluto.com/tle_info.htm.)
This helps when accurate positions are not provided or completely trusted
or appear to be bad.

   This code will read the 80-column astrometry and,  when it finds a
spacecraft position report (the "second line"),  look through the TLE
file for a matching TLE.  If it finds one,  it'll compute the spacecraft
location for that time and modify the position accordingly.

   This was written specifically to handle a problem with a satellite in
LEO.  TLEs don't exist for heliocentric objects,  but I may revise this
code eventually to add positions for them (this could be done by,  for
example,  requesting them from JPL Horizons).  TLEs exist for TESS,
computed by me,  but some work would be required to make those actually
function properly (the numbers are larger and the format is somewhat
different as a result).  So this really works,  at present,  only for
(C51) WISE, (C52) Swift,  and (C53) NEOSSat.             */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include "norad.h"
#include "watdefs.h"
#include "afuncs.h"
#include "mpc_func.h"
#include "stringex.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

int find_tle( tle_t *tle, const char *filename, const int norad_no)
{
   FILE *ifile = fopen( filename, "rb");
   char line0[100], line1[100], line2[100];
   int rval = -1;

   if( !ifile)
      {
      fprintf( stderr, "Couldn't open TLE file '%s'\n", filename);
      exit( -1);
      }
   *line0 = *line1 = '\0';
   while( rval && fgets( line2, sizeof( line2), ifile))
      {
      if( *line1 == '1' && *line2 == '2' && atoi( line1 + 2) == norad_no)
         if( parse_elements( line1, line2, tle) >= 0)
            rval = 0;
      strlcpy_error( line0, line1);
      strlcpy_error( line1, line2);
      }
   fclose( ifile);
   if( rval)
      {
      fprintf( stderr, "Couldn't find TLEs for %5d in TLE file '%s'\n",
                           norad_no, filename);
      exit( -1);
      }
   return( rval);
}

/*
C49 = 29510 = 2006-047A = STEREO-A
C50 = 29511 = 2006-047B = STEREO-B
C51 = 36119 = 2009-071A = WISE   *
C52 = 28485 = 2004-047A = Swift   *
C53 = 39089 = 2013-009D = NEOSSat  *
C54 = 28928 = 2006-001A = New Horizons
C55 = 34380 = 2009-011A = Kepler
C56 = 41043 = 2015-070A = LISA-Pathfinder
C57 = 43435 = 2018-038A = TESS      *
   Note that only the asterisked objects are actually in earth orbit
and therefore have TLEs.  But I may try to figure out some way to add
in positions for the heliocentric spacecraft later.  */

static int mpc_code_to_norad_number( const char *mpc_code)
{
   const char *codes = "C49 C50 C51 C52 C53 C54 C55 C56 C57 ";
   const int norad_numbers[] = { 29510, 29511, 36119, 28485, 39089,
                                        28928, 34380, 41043, 43435 };
   int i;

   for( i = 0; codes[i * 4]; i++)
      if( !memcmp( codes + i * 4, mpc_code, 3))
         return( norad_numbers[i]);
   assert( 1);          /* not supposed to ever happen,  unless a new */
   return( 0);          /* satellite has been added */
}

int main( const int argc, const char **argv)
{
   const char *tle_filename = (argc > 2 ? argv[2] : "all_tle.txt");
   const char *astrometry_filename = argv[1];
   FILE *ifile = (argc > 1 ? fopen( astrometry_filename, "rb") : NULL);
   int is_deep = 0, curr_norad = 0;
   char buff[200];
   double params[N_SAT_PARAMS];
   const double J2000 = 2451545.;
   tle_t tle;

   if( argc > 1 && !ifile)
      fprintf( stderr, "'%s' not found : %s\n", astrometry_filename, strerror( errno));
   if( !ifile)
      {
      fprintf( stderr, "'line2' takes as a command-line argument the name of the input\n"
                       "astrometry file.  It sends astrometry with corrected/added\n"
                       "spacecraft locations to stdout.\n");
      exit( -1);
      }
   while( fgets( buff, sizeof( buff), ifile))
      {
      double jd;

      if( strlen( buff) > 80 && buff[14] == 's'
            && (jd = extract_date_from_mpc_report( buff, NULL)) > 2400000.)
         {
         double state[6], state_j2000[6], precess_matrix[9];
         double t_since;
         const int norad_number = mpc_code_to_norad_number( buff + 77);
         int i;

         if( curr_norad != norad_number)
            {
            curr_norad = norad_number;
            find_tle( &tle, tle_filename, norad_number);
            is_deep = select_ephemeris( &tle);
            if( is_deep)
               SDP4_init( params, &tle);
            else
               SGP4_init( params, &tle);
            }
         t_since = (jd - tle.epoch) * minutes_per_day;
         if( is_deep)
            SDP4( t_since, &tle, params, state, state + 3);
         else
            SGP4( t_since, &tle, params, state, state + 3);
         setup_precession( precess_matrix, 2000. + (jd - J2000) / 365.25, 2000.);
         precess_vector( precess_matrix, state, state_j2000);
         precess_vector( precess_matrix, state + 3, state_j2000 + 3);
         for( i = 0; i < 3; i++)
            {
            char *tptr = buff + 34 + i * 12;

            snprintf_err( tptr, 12, "%11.4f", fabs( state_j2000[i]));
            *tptr = (state_j2000[i] > 0. ? '+' : '-');
            tptr[11] = ' ';
            }
         }
      printf( "%s", buff);
      }
   fclose( ifile);
   return( 0);
}
