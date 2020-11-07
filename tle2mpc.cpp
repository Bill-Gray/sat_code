/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "watdefs.h"
#include "norad.h"
#include "afuncs.h"

/* Code to read TLEs of the sort I post on github.com/Bill-Gray/tles,
where there's one TLE per day,  and output ephemerides of the sort
the Minor Planet Center wants for their Distant Artsat Observation
Page (DASO),  at https://www.minorplanetcenter.net/iau/artsats/artsats.html.
For DASO,  we need geocentric ephems with positions only,  in equatorial
J2000,  in AU,  at 0.1-day intervals.  So each time we read a TLE,  ten
positions are output.         */

static void error_exit( void)
{
   fprintf( stderr, "'tle2mpc' needs the name of an input TLE file and an output\n");
   fprintf( stderr, "MPC ephemeris file.  See 'tle2mpc.cpp' for details.\n");
   exit( -1);
}

static const char *err_fgets( char *buff, size_t buffsize, FILE *ifile)
{
   const char *rval = fgets( buff, (int)buffsize, ifile);

   assert( rval);
   return( rval);
}

int main( const int argc, const char **argv)
{
   char buff[200];
   FILE *ifile, *ofile;
   bool object_name_shown = false;

   if( argc < 3)
      error_exit( );
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      fprintf( stderr,  "Couldn't open input file '%s'\n", argv[1]);
      error_exit( );
      }
   ofile = fopen( argv[2], "wb");
   if( !ofile)
      {
      fprintf( stderr,  "Couldn't open output file '%s'\n", argv[2]);
      error_exit( );
      }
   while( fgets( buff, sizeof( buff), ifile))
      if( !memcmp( buff, "# Ephem range:", 14))
         {
         const double mjd1 = atof( buff + 15);
         const double mjd2 = atof( buff + 28);

         fprintf( ofile, "%f 0.1 %d 0,1,1 ", 2400000.5 + mjd1,
                           (int)( (mjd2 - mjd1) * 10. + .5));
         }
      else if( !memcmp( buff, "# MJD ", 6))
         {
         double mjd = atof( buff + 6);
         char line2[80];
         tle_t tle;  /* Structure for two-line elements set for satellite */
         bool is_a_tle;

         err_fgets( buff, sizeof( buff), ifile);
         if( !object_name_shown)
            fprintf( ofile, "%s", buff);
         object_name_shown = true;
         err_fgets( buff, sizeof( buff), ifile);
         err_fgets( line2, sizeof( line2), ifile);
         is_a_tle = (parse_elements( buff, line2, &tle) >= 0);
         assert( is_a_tle);
         if( is_a_tle)
            {
            int i, j;
            double sat_params[N_SAT_PARAMS], precess_matrix[9];
            const bool is_deep = select_ephemeris( &tle);
            const double year = 2000. + (mjd - 51545.) / 365.25;

            if( is_deep)
               SDP4_init( sat_params, &tle);
            else
               SGP4_init( sat_params, &tle);
            setup_precession( precess_matrix, year, 2000.);
            for( i = 0; i < 10; i++, mjd += 0.1)
               {
               double posn[3], j2000_posn[3];
               const double t_since = (mjd - (tle.epoch - 2400000.5)) * minutes_per_day;

               if( is_deep)
                  SDP4( t_since, &tle, sat_params, posn, NULL);
               else
                  SGP4( t_since, &tle, sat_params, posn, NULL);
               for( j = 0; j < 3; j++)
                  posn[j] /= AU_IN_KM;
               precess_vector( precess_matrix, posn, j2000_posn);
               fprintf( ofile, "%.5f%16.10f%16.10f%16.10f\n",
                    mjd + 2400000.5, j2000_posn[0], j2000_posn[1], j2000_posn[2]);
               }
            }
         }
               /* Rewind to start of input & re-read,  looking for comments */
   fseek( ifile, 0L, SEEK_SET);
   while( fgets( buff, sizeof( buff), ifile))
      if( !memcmp( buff, "# Ephem range: ", 14))
         while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "# 1 NoradU", 10))
            fprintf( ofile, "%s", buff);

   fclose( ifile);
   fclose( ofile);
   return( 0);
}
