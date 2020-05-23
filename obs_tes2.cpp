/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
    obs_tes2.cpp     12 December 2002

   (Revised slightly December 2012 to fix compiler warning errors.)

   An example 'main' function illustrating how to find which satellite(s)
are within a given radius of a given RA/dec,  as seen from a given
point.  The code reads in a TLE file (name provided as the first
command-line argument).  Details of the observer position,  search
radius,  date/time,  and RA/dec are provided on the command line.
For example:

obs_tes2 alldat.tle -l44.01,-69.9,10 -p90,30 -j2452623.5 -r10

   would hunt through the TLE element file 'alldat.tle' for satellites
visible from latitude +44.01,  longitude -69.9,  altitude 10 metres;
at RA=90 degrees (6h),  dec=+30;  on JD 2452623.5 (UTC);  within a
ten-degree search radius.  (All of these are the default values.)
The output looks like this:

NORAD  Int'l     RA (J2000) dec    Delta Radius  PA Speed
08593U 74089DG  88.7235  22.9622   2293.0  7.15 225 10.75
15830U 85049D   82.5051  34.9711  32143.6  8.99  34  0.24
17642U 81053LQ  88.1471  32.6585   1428.5  3.24 213 17.60
21833U 91088A   80.6400  27.9649  36216.6  9.58  87  0.17

   ...with 'delta'=distance to satellite in km,  'radius'=angular
distance in degrees from the search point,  'PA' = position angle
of motion, 'Speed' = apparent angular rate of motion in
arcminutes/second (or degrees/minute). */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "norad.h"
#include "observe.h"

#define PI 3.141592653589793238462643383279
#define TIME_EPSILON (1./86400.)

int main( const int argc, const char **argv)
{
   const char *tle_file_name = ((argc == 1) ? "alldat.tle" : argv[1]);
   FILE *ifile = fopen( tle_file_name, "rb");
   char line1[100], line2[100];
   double lat = 44.01, lon = -69.9, ht_in_meters = 10.;
   double jd = 2452623.5;   /* 15 Dec 2002 0h UT */
   double search_radius = 10.;     /* default to ten-degree search */
   double target_ra = 90., target_dec = 30.;  /* default search is at RA=6h, dec=+30 */
   double rho_sin_phi, rho_cos_phi, observer_loc[3], observer_loc2[3];
   int i, header_line_shown = 0;

   if( !ifile)
      {
      printf( "Couldn't open input file %s\n", tle_file_name);
      exit( -1);
      }

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'l':
               sscanf( argv[i] + 2, "%lf,%lf,%lf", &lat, &lon, &ht_in_meters);
               break;
            case 'p':
               sscanf( argv[i] + 2, "%lf,%lf", &target_ra, &target_dec);
               break;
            case 'j':
               jd = atof( argv[i] + 2);
               break;
            case 'r':
               search_radius = atof( argv[i] + 2);
               break;
            default:
               printf( "Unrecognized command-line option '%s'\n", argv[i]);
               exit( -2);
               break;
            }

            /* Figure out where the observer _really_ is,  in Cartesian */
            /* coordinates of date: */
   earth_lat_alt_to_parallax( lat * PI / 180., ht_in_meters, &rho_cos_phi,
                                                             &rho_sin_phi);
   observer_cartesian_coords( jd,
                lon * PI / 180., rho_cos_phi, rho_sin_phi, observer_loc);
   observer_cartesian_coords( jd + TIME_EPSILON,
                lon * PI / 180., rho_cos_phi, rho_sin_phi, observer_loc2);
   target_ra *= PI / 180.;
   target_dec *= PI / 180.;

   if( fgets( line1, sizeof( line1), ifile))
      while( fgets( line2, sizeof( line2), ifile))
         {
         tle_t tle;     /* Structure for two-line elements set for satellite */

         if( !parse_elements( line1, line2, &tle))    /* hey! we got a TLE! */
            {
            int is_deep = select_ephemeris( &tle);
            double sat_params[N_SAT_PARAMS], radius, d_ra, d_dec;
            double ra, dec, dist_to_satellite, t_since;
            double pos[3]; /* Satellite position vector */
            double unused_delta2;

            t_since = (jd - tle.epoch) * 1440.;
            if( is_deep)
               {
               SDP4_init( sat_params, &tle);
               SDP4( t_since, &tle, sat_params, pos, NULL);
               }
            else
               {
               SGP4_init( sat_params, &tle);
               SGP4( t_since, &tle, sat_params, pos, NULL);
               }
            get_satellite_ra_dec_delta( observer_loc, pos,
                                    &ra, &dec, &dist_to_satellite);
            epoch_of_date_to_j2000( jd, &ra, &dec);
            d_ra = (ra - target_ra + PI * 4.);
            while( d_ra > PI)
               d_ra -= PI + PI;
            d_dec = dec - target_dec;
            radius = sqrt( d_ra * d_ra + d_dec * d_dec) * 180. / PI;
            if( radius < search_radius)      /* good enough for us! */
               {
               double speed, posn_ang_of_motion;

               line1[16] = '\0';

               if( !header_line_shown)
                  {
                  printf( "NORAD  Int'l     RA (J2000) dec    Delta Radius  PA Speed\n");
                  header_line_shown = 1;
                  }
                                 /* Compute position one second later,  so we */
                                 /* can show speed/PA of motion: */
               t_since += TIME_EPSILON * 1440.;
               if( is_deep)
                  SDP4( t_since, &tle, sat_params, pos, NULL);
               else
                  SGP4( t_since, &tle, sat_params, pos, NULL);
               get_satellite_ra_dec_delta( observer_loc2, pos,
                                       &d_ra, &d_dec, &unused_delta2);
               epoch_of_date_to_j2000( jd, &d_ra, &d_dec);
               d_ra -= ra;
               d_dec -= dec;
               while( d_ra > PI)
                  d_ra -= PI + PI;
               while( d_ra < -PI)
                  d_ra += PI + PI;
               d_ra *= cos( dec);
               posn_ang_of_motion = atan2( d_ra, d_dec);
               if( posn_ang_of_motion < 0.)
                  posn_ang_of_motion += PI + PI;
               speed = sqrt( d_ra * d_ra + d_dec * d_dec) * 180. / PI;
                        /* Put RA into 0 to 2pi range: */
               ra = fmod( ra + PI * 10., PI + PI);
               printf( "%s %8.4f %8.4f %8.1f %5.2f %3d %5.2f\n",
                        line1 + 2, ra * 180. / PI, dec * 180. / PI,
                        dist_to_satellite, radius,
                        (int)(posn_ang_of_motion * 180 / PI),
                        speed * 60.);
                              /* "Speed" is displayed in arcminutes/second,
                                 or in degrees/minute */
               }
            }
         strcpy( line1, line2);
         }
   fclose( ifile);
   return( 0);
} /* End of main() */

