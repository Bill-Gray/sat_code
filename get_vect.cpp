/* Copyright (C) 2018, Project Pluto.  See LICENSE.

   This code can read a TLE and compute a geocentric state vector,
formatted such that Find_Orb can then read it in.  I did this
partly to test the hypothesis that if you compute a state vector
from Space-Track TLEs at their epoch,  you get the "actual" motion.
That is to say,  you could numerically integrate it to get a better
result.  This turns out not to be the case.  Space-Track TLEs may
be a best-fit to a set of observations or (as with my own TLEs)
a best fit to a numerically integrated ephemeris,  but there
doesn't seem to be a way to improve them by doing a numerical
integration.

   My second purpose was to be able to feed the state vector created
by this program into Find_Orb as an initial orbit guess.  For that
purpose,  it seems to work.  You see large residuals as a result
of the difference between numerical integration and SGP4/SDP4.
But it gets you close enough that you can then do differential
corrections (least-squares fitting).         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "norad.h"
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   FILE *ifile;
   const char *filename = "all_tle.txt";
   char line0[100], line1[100], line2[100];
   int i;
   const char *norad = NULL, *intl = NULL;
   double jd = 0.;
   int is_j2000 = 1, is_equatorial = 1;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *arg = (i < argc - 1 && !argv[i][2]
                                          ? argv[i + 1] : argv[i] + 2);

         switch( argv[i][1])
            {
            case 'n':
               norad = arg;
               printf( "Looking for NORAD %s\n", norad);
               break;
            case 'i':
               intl = arg;
               printf( "Looking for international ID %s\n", intl);
               break;
            case 't':
               jd = get_time_from_string( 0., arg, FULL_CTIME_YMD, NULL);
               break;
            case 'd':
               is_j2000 = 0;
               break;
            default:
               printf( "'%s': unrecognized option\n", argv[i]);
               return( -1);
               break;
            }
         }
   if( argc > 1 && argv[1][0] != '-')
      filename = argv[1];
   ifile = fopen( filename, "rb");
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open '%s': ", filename);
      perror( "");
      return( -1);
      }
   *line0 = *line1 = '\0';
   while( fgets( line2, sizeof( line2), ifile))
      {
      if( *line1 == '1' && (!norad || !memcmp( line1 + 2, norad, 5))
                        && (!intl || !memcmp( line1 + 9, intl, strlen( intl)))
                   && *line2 == '2')
         {
         tle_t tle;
         const int err_code = parse_elements( line1, line2, &tle);

         if( err_code >= 0)
            {
            const int is_deep = select_ephemeris( &tle);
            double state[6], state_j2000[6], precess_matrix[9];
            double params[N_SAT_PARAMS], t_since;
            const double epoch_tdt =
                        tle.epoch + td_minus_utc( tle.epoch) / seconds_per_day;
            const double J2000 = 2451545.;
            double *state_to_show;

            if( !jd)
               jd = epoch_tdt;
            t_since = (jd - epoch_tdt) * minutes_per_day;
            if( is_deep)
               {
               SDP4_init( params, &tle);
               SDP4( t_since, &tle, params, state, state + 3);
               }
            else
               {
               SGP4_init( params, &tle);
               SGP4( t_since, &tle, params, state, state + 3);
               }
            if( strlen( line0) < 60)
               printf( "%s", line0);
            setup_precession( precess_matrix, 2000. + (jd - J2000) / 365.25, 2000.);
            precess_vector( precess_matrix, state, state_j2000);
            precess_vector( precess_matrix, state + 3, state_j2000 + 3);
            state_to_show = (is_j2000 ? state_j2000 : state);
            printf( " %.6f %.6s\n", jd, line1 + 9);
            printf( " %.5f %.5f %.5f 0408   # Ctr 3 km sec %s %s\n",
                           state_to_show[0], state_to_show[1], state_to_show[2],
                           is_equatorial ? "eq" : "ecl",
                           is_j2000 ? "" : "of_date");
            printf( " %.5f %.5f %.5f 0 0 0\n",
                           state_to_show[3] / seconds_per_minute,
                           state_to_show[4] / seconds_per_minute,
                           state_to_show[5] / seconds_per_minute);
            }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( 0);
}
