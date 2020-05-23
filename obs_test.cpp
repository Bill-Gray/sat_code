/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
    obs_test.cpp     23 September 2002

   (Revised slightly December 2012 to fix compiler warning errors.)

   An example 'main' function illustrating how to get topocentric
ephemerides for artificial satellites,  using the basic satellite
code plus the add-on topocentric functions.  The code reads the
file 'obs_test.txt',  getting commands setting the observer lat/lon
and altitude and time of observation.  When it gets a command
setting a particular JD,  it computes the topocentric RA/dec/dist
and prints them out.

   At present,  'obs_test.txt' sets up a lat/lon/height in Bowdoinham,
Maine,  corporate headquarters of Project Pluto,  and computes the
position of one low-orbit satellite (ISS) and one high-orbit satellite
(Cosmos 1966 rocket booster).  You should get :

Near-Earth type Ephemeris (SGP4) selected:
Object 25544U 98067A, as seen from lat 44.01000 lon -69.90000, JD 2452541.50000
RA 350.1615 (J2000) dec -24.0241 dist 1867.97481 km
Deep-Space type Ephemeris (SDP4) selected:
Object 19448U 88076D, as seen from lat 44.01000 lon -69.90000, JD 2452541.50000
RA 3.5743 (J2000) dec 30.4293 dist 32114.83370 km
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "norad.h"
#include "observe.h"

#define PI 3.141592653589793238462643383279

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( (argc == 1) ? "obs_test.txt" : argv[1], "rb");
   tle_t tle; /* Pointer to two-line elements set for satellite */
   char line1[100], line2[100];
   double lat = 0., lon = 0., ht_in_meters = 0., jd = 0.;
   int ephem = 1;       /* default to SGP4 */

   if( !ifile)
      {
      printf( "Couldn't open input OBS_TEST.TXT file\n");
      exit( -1);
      }
   if( fgets( line1, sizeof( line1), ifile))
      while( fgets( line2, sizeof( line2), ifile))
         {
         int err_val;

         if( !memcmp( line2, "Ephem ", 6))
            ephem = (line2[6] - '0');
         else if( !memcmp( line2, "JD ", 3))
            jd = atof( line2 + 3);
         else if( !memcmp( line2, "ht ", 3))
            ht_in_meters = atof( line2 + 3);
         else if( !memcmp( line2, "lat ", 4))
            lat = atof( line2 + 4) * PI / 180.;     /* cvt degrees to radians */
         else if( !memcmp( line2, "lon ", 4))
            lon = atof( line2 + 4) * PI / 180.;
         else if( (err_val = parse_elements( line1, line2, &tle)) >= 0)
            {                  /* hey! we got a TLE! */
            int is_deep = select_ephemeris( &tle);
            const char *ephem_names[5] = { "SGP ", "SGP4", "SGP8", "SDP4", "SDP8" };
            double sat_params[N_SAT_PARAMS], observer_loc[3];
            double rho_sin_phi, rho_cos_phi;
            double ra, dec, dist_to_satellite, t_since;
            double pos[3]; /* Satellite position vector */

            if( err_val)
               printf( "WARNING: TLE parsing error %d\n", err_val);
            earth_lat_alt_to_parallax( lat, ht_in_meters, &rho_cos_phi,
                                                          &rho_sin_phi);
            observer_cartesian_coords( jd, lon, rho_cos_phi, rho_sin_phi,
                        observer_loc);
            if( is_deep && (ephem == 1 || ephem == 2))
               ephem += 2;    /* switch to an SDx */
            if( !is_deep && (ephem == 3 || ephem == 4))
               ephem -= 2;    /* switch to an SGx */
            if( is_deep)
               printf("Deep-Space type Ephemeris (%s) selected:\n",
                                          ephem_names[ephem]);
            else
               printf("Near-Earth type Ephemeris (%s) selected:\n",
                                          ephem_names[ephem]);

            /* Calling of NORAD routines */
            /* Each NORAD routine (SGP, SGP4, SGP8, SDP4, SDP8)   */
            /* will be called in turn with the appropriate TLE set */
            t_since = (jd - tle.epoch) * 1440.;
            switch( ephem)
               {
               case 0:
                  SGP_init( sat_params, &tle);
                  err_val = SGP( t_since, &tle, sat_params, pos, NULL);
                  break;
               case 1:
                  SGP4_init( sat_params, &tle);
                  err_val = SGP4( t_since, &tle, sat_params, pos, NULL);
                  break;
               case 2:
                  SGP8_init( sat_params, &tle);
                  err_val = SGP8( t_since, &tle, sat_params, pos, NULL);
                  break;
               case 3:
                  SDP4_init( sat_params, &tle);
                  err_val = SDP4( t_since, &tle, sat_params, pos, NULL);
                  break;
               case 4:
                  SDP8_init( sat_params, &tle);
                  err_val = SDP8( t_since, &tle, sat_params, pos, NULL);
                  break;
               default:
                  printf( "? How did we get here? ephem = %d\n", ephem);
                  err_val = 0;
                  break;
               }
            if( err_val)
               printf( "Ephemeris error %d\n", err_val);
            line1[15] = '\0';
            printf( "Object %s, as seen from lat %.5f lon %.5f, JD %.5f\n",
                     line1 + 2, lat * 180. / PI, lon * 180. / PI, jd);
            get_satellite_ra_dec_delta( observer_loc, pos,
                                    &ra, &dec, &dist_to_satellite);
            epoch_of_date_to_j2000( jd, &ra, &dec);
            printf( "RA %.4f (J2000) dec %.4f dist %.5f km\n",
                        ra * 180. / PI, dec * 180. / PI, dist_to_satellite);
            }
         strcpy( line1, line2);
         }
   fclose( ifile);
   return( 0);
} /* End of main() */

