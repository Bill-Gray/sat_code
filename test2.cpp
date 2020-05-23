/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
 *  test2.cpp      5 July 2002
 *
 *  A skeleton main() function to demonstrate the use of
 *  the various NORAD ephemerides.  It reads in the file
 *  'test.tle' and computes ephemerides at assorted times for
 *  all elements in the file.  The times used,  and the choice
 *  of SxP4 vs. SxP8,  can be switched with various keywords in
 *  the 'test.tle' file.  The result should match 'test2.txt'.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "norad.h"

/* Main program */
int main( int argc, char **argv)
{
   double vel[3], pos[3]; /* Satellite position and velocity vectors */
   FILE *ifile = fopen( (argc == 1) ? "test.tle" : argv[1], "rb");
   tle_t tle; /* Pointer to two-line elements set for satellite */
   char line1[100], line2[100];
   int n_times = 5, dundee_output = 0;
   double times[400];
   int ephem = TLE_EPHEMERIS_TYPE_SGP4;       /* default to SGP4 */
   int i;               /* Index for loops etc */

   for( i = 2; i < argc; i++)
      switch( argv[i][0])
         {
         case 's': case 'S':
            sxpx_set_dpsec_integration_step( atof( argv[i] + 1));
            break;
         case 'd': case 'D':
            sxpx_set_implementation_param( SXPX_DUNDEE_COMPLIANCE, 1);
            break;
         case 'o': case 'O':
            {
            const int new_order = atoi( argv[i] + 1);

            sxpx_set_implementation_param(
                               SXPX_DPSEC_INTEGRATION_ORDER, new_order);
            printf( "Set to order %d\n", new_order);
            }
            break;
         default:
            break;
         }

   if( !ifile)
      {
      printf( "Couldn't open input TLE file\n");
      exit( -1);
      }
   for( i = 0; i < n_times; i++)
      times[i] = (double)(i * 30);
   if( fgets( line1, sizeof( line1), ifile))
      while( fgets( line2, sizeof( line2), ifile))
         {
         if( !strncmp( line2, "Ephem ", 6))
            ephem = (line2[6] - '0');
         else if( !strncmp( line2, "Dundee ", 7))
            {
            double t0, step;

            sscanf( line2 + 7, "%lf %lf %d", &t0, &step, &n_times);
            for( i = 0; i < n_times; i++)
               times[i] = t0 + (double)i * step;
            sxpx_set_implementation_param( SXPX_DUNDEE_COMPLIANCE, 1);
            dundee_output = 1;
            }
         else if( !strncmp( line2, "Times: ", 7))
            {
            int loc = 7, bytes_read;

            n_times = 0;
            while( sscanf( line2 + loc, "%lf%n", times + n_times, &bytes_read) == 1)
               {
               loc += bytes_read;
               n_times++;
               }
            }
         else if( parse_elements( line1, line2, &tle) >= 0)
            {                              /* hey! we got a TLE! */
            int is_deep = select_ephemeris( &tle);
            const char *ephem_names[6] = { "???", "SGP ", "SGP4", "SDP4", "SGP8", "SDP8" };
            double sat_params[N_SAT_PARAMS];

            if( is_deep)
               if( ephem == TLE_EPHEMERIS_TYPE_SGP4 || ephem == TLE_EPHEMERIS_TYPE_SGP8)
                  ephem++;    /* switch to an SDPx model */
            if( !is_deep)
               if( ephem == TLE_EPHEMERIS_TYPE_SDP4 || ephem == TLE_EPHEMERIS_TYPE_SDP8)
                  ephem--;    /* switch to an SGPx model */
            line1[69] = line2[69] = '\0';
            if( dundee_output)
               printf( "#%s\n#%s\n", line1, line2);
            else
               printf( "%s\n%s\n", line1, line2);
            if( is_deep)
               printf("Deep-Space type Ephemeris (%s) selected:",
                                          ephem_names[ephem]);
            else
               printf("Near-Earth type Ephemeris (%s) selected:",
                                          ephem_names[ephem]);

            /* Print some titles for the results */
            printf("\nEphem:%s   Tsince         "
                "X/Xdot           Y/Ydot           Z/Zdot\n", ephem_names[ephem]);

            /* Calling of NORAD routines */
            /* Each NORAD routine (SGP, SGP4, SGP8, SDP4, SDP8)   */
            /* will be called in turn with the appropriate TLE set */
            switch( ephem)
               {
               case TLE_EPHEMERIS_TYPE_SGP:
                  SGP_init( sat_params, &tle);
                  break;
               case TLE_EPHEMERIS_TYPE_SGP4:
                  SGP4_init( sat_params, &tle);
                  break;
               case TLE_EPHEMERIS_TYPE_SGP8:
                  SGP8_init( sat_params, &tle);
                  break;
               case TLE_EPHEMERIS_TYPE_SDP4:
                  SDP4_init( sat_params, &tle);
                  break;
               case TLE_EPHEMERIS_TYPE_SDP8:
                  SDP8_init( sat_params, &tle);
                  break;
               }

            for( i = 0; i < n_times; i++)
               {
               switch( ephem)
                  {
                  case TLE_EPHEMERIS_TYPE_SGP:
                     SGP(times[i], &tle, sat_params, pos, vel);
                     break;
                  case TLE_EPHEMERIS_TYPE_SGP4:
                     SGP4(times[i], &tle, sat_params, pos, vel);
                     break;
                  case TLE_EPHEMERIS_TYPE_SGP8:
                     SGP8(times[i], &tle, sat_params, pos, vel);
                     break;
                  case TLE_EPHEMERIS_TYPE_SDP4:
                     SDP4(times[i], &tle, sat_params, pos, vel);
                     break;
                  case TLE_EPHEMERIS_TYPE_SDP8:
                     SDP8(times[i], &tle, sat_params, pos, vel);
                     break;
                  }

               /* Calculate and print results */
               vel[0] /= 60.;    /* cvt km/minute to km/second */
               vel[1] /= 60.;
               vel[2] /= 60.;

               if( dundee_output)
                  {
                  printf("%12.4f   %16.8f %16.8f %16.8f",
                                  times[i],pos[0],pos[1],pos[2]);
                  printf(" %16.8f %16.8f %16.8f\n",
                                           vel[0],vel[1],vel[2]);
                  }
               else
                  {
                  printf("%5d  %12.4f   %16.8f %16.8f %16.8f\n", tle.norad_number,
                                  times[i],pos[0],pos[1],pos[2]);
                  printf("                      %16.8f %16.8f %16.8f\n",
                                         vel[0],vel[1],vel[2]);
                  }
               } /* End of for( i = 0; i < n_times; i++)  */
            printf( "\n");
            }
         strcpy( line1, line2);
         }
  fclose( ifile);
  return(0);
} /* End of main() */

/*------------------------------------------------------------------*/
