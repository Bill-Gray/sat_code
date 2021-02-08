/* Copyright (C) 2018, Project Pluto.  See LICENSE.

   This program will generate simulated geocentric observations
for a given object from a TLE.  In theory,  one can then fit these
pseudo-observations to a higher-quality physical model to get a
considerably more accurate ephemeris for the object.  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "norad.h"
#include "observe.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char line1[100], line2[100];
   const char *intl_id = NULL;
   double step_size = .1;
   int i, n_steps = 100;
   bool show_vectors = false;

   if( !ifile)
      {
      printf( "Couldn't open input file\n");
      exit( -1);
      }
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'i':
               intl_id = argv[i] + 2;
               break;
            case 'n':
               n_steps = atoi( argv[i] + 2);
               break;
            case 's':
               step_size = atof( argv[i] + 2);
               break;
            case 'v':
               show_vectors = true;
               break;
            default:
               printf( "Unrecognized option '%s'\n", argv[i]);
               break;
            }
   *line1 = '\0';
   sxpx_set_implementation_param( SXPX_DUNDEE_COMPLIANCE, 1);
   while( fgets( line2, sizeof( line2), ifile))
      {
      tle_t tle; /* Pointer to two-line elements set for satellite */
      int err_val;

      if( (!intl_id || !memcmp( intl_id, line1 + 9, 6))
                && (err_val = parse_elements( line1, line2, &tle)) >= 0)
         {                  /* hey! we got a TLE! */
         int is_deep = select_ephemeris( &tle);
         double sat_params[N_SAT_PARAMS], observer_loc[3];
         double prev_pos[3];

         if( err_val)
            printf( "WARNING: TLE parsing error %d\n", err_val);
         for( i = 0; i < 3; i++)
            observer_loc[i] = 0.;
         if( is_deep)
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         for( i = 0; i < n_steps; i++)
            {
            double pos[3]; /* Satellite position vector */
            double t_since = (double)( i - n_steps / 2) * step_size;
            double jd = tle.epoch + t_since;

            t_since *= 1440.;
            if( is_deep)
               err_val = SDP4( t_since, &tle, sat_params, pos, NULL);
            else
               err_val = SGP4( t_since, &tle, sat_params, pos, NULL);
            if( err_val)
               printf( "Ephemeris error %d\n", err_val);
            if( show_vectors)
               {
               if( i)
                  printf( "%14.6f %14.6f %14.6f - ", pos[0] - prev_pos[0],
                                                 pos[1] - prev_pos[1],
                                                 pos[2] - prev_pos[2]);
               printf( "%14.6f %14.6f %14.6f\n", pos[0], pos[1], pos[2]);
               memcpy( prev_pos, pos, 3 * sizeof( double));
               }
            else
               {
               double ra, dec, dist_to_satellite;

               get_satellite_ra_dec_delta( observer_loc, pos,
                                 &ra, &dec, &dist_to_satellite);
               epoch_of_date_to_j2000( jd, &ra, &dec);
               printf( "%-14sC%13.5f    %08.4f    %+08.4f",
                     intl_id, jd, ra * 180. / PI, dec * 180. / PI);
               printf( "                    TLEs 500\n");
               }
            }
         }
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( 0);
} /* End of main() */

