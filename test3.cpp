/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
 *  test3.cpp      8 October 2002
 *
 *  A skeleton main() function to test the speed of SxPx
 *  functions,  and the positional differences between SGPx and SDPx.
 *  (I was curious as to whether the extra code in SDPx made any really
 *  significant difference;  it was starting to look as if SDPx was a
 *  waste of code.  Turns out that it depends greatly on the satellite
 *  in question;  I'll post results in a bit,  when I find the time.)
 *
 *  Also,  this demonstrates how to use the new dynamically-loaded
 *  functions in 'dynamic.cpp'.  This way,  if your program can't
 *  find 'sat_code.dll',  it'll gracefully give you a message about it.
 *  In the case of my own _Guide_ planetarium software,  if it can't
 *  find the DLL,  it can fall back on its own SGP4-only code.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "norad.h"

/* Main program */
int main( int argc, char **argv)
{
   FILE *ifile = fopen( (argc == 1) ? "test.tle" : argv[1], "rb");
   tle_t tle; /* Pointer to two-line elements set for satellite */
   char line1[100], line2[100];
   double t_since = atof( argv[2]);
   double rms_diff = 0.;
   int verbose = 0, n_found = 0, n_runs = 1, show_sdp = 1;
   int ephem = 1;       /* default to SGP4 */
   int i, j;            /* Index for loops etc */
   clock_t t0;

   if( !ifile)
      {
      printf( "Couldn't open input TLE file\n");
      exit( -1);
      }

   for( i = 3; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v':
               verbose = atoi( argv[i] + 2);
               break;
            case 'n':
               n_runs = atoi( argv[i] + 2);
               break;
            case 's':
               show_sdp = atoi( argv[i] + 2);
               break;
            default:
               printf( "Option '%s' ignored\n", argv[i]);
               break;
            }

   fgets( line1, sizeof( line1), ifile);
   t0 = clock( );
   while( fgets( line2, sizeof( line2), ifile))
      {
      if( line1[0] == '1' && line2[0] == '2' &&
               parse_elements( line1, line2, &tle) >= 0)  /* hey! we got a TLE! */
         if( select_ephemeris( &tle) && show_sdp < 2)
            {
            double sgp_sat_params[N_SAT_PARAMS];
            double sdp_sat_params[N_SAT_PARAMS];
            double perigee;

#ifdef STATIC_LINK
            SGP4_init( sgp_sat_params, &tle);
            if( show_sdp)
               SDP4_init( sdp_sat_params, &tle);
#else
            if( SXPX_init( sgp_sat_params, &tle, 1))
               {
               printf( "Couldn't load 'sat_code.dll'\n");
               exit( -1);
               }
            if( show_sdp)
               SXPX_init( sdp_sat_params, &tle, 3);
#endif
            perigee = sgp_sat_params[9] * (1. - tle.eo) - 1.;
            line2[63] = '\0';
            if( verbose)
               printf( "%5ld %s", (long)( perigee * 6378.14), line2);
            for( j = 0; j <= n_runs; j++)
               {
               double posn2[3], pos[3], vel[3];
               double delta, d2 = 0.;

#ifdef STATIC_LINK
               SGP4( t_since * (double)j, &tle, sgp_sat_params, posn2, vel);
               if( show_sdp)
                  SDP4( t_since * (double)j, &tle, sdp_sat_params, pos, vel);
#else
               SXPX( t_since * (double)j, &tle, sgp_sat_params, posn2, vel, 1);
               if( show_sdp)
                  SXPX( t_since * (double)j, &tle, sdp_sat_params, pos, vel, 3);
#endif
               if( !show_sdp)
                  memcpy( pos, posn2, 3 * sizeof( double));

               for( i = 0; i < 3; i++)
                  {
                  delta = pos[i] - posn2[i];
                  d2 += delta * delta;
                  }
               if( j == 1)
                  rms_diff += d2;
               if( verbose)
                  printf( "%8.1lf%s", sqrt( d2), (j == n_runs ? "\n" : ""));
               }
            n_found++;
            }
      strcpy( line1, line2);
      }
   fclose( ifile);
   t0 = clock( ) - t0;
   printf( "%d deep-space objects found\n", n_found);
   printf( "%lf seconds elapsed\n", (double)t0 / (double)CLOCKS_PER_SEC);
   if( n_found)
      printf( "RMS difference = %.1lf\n", sqrt( rms_diff / (double)n_found));
   return(0);
} /* End of main() */

/*------------------------------------------------------------------*/
