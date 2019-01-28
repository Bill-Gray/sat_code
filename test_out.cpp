/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "norad.h"
#include "norad_in.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

/* Example code to add BSTAR data using Ted Molczan's method.  It just
   reads in TLEs,  computes BSTAR if possible,  then writes out the
   resulting modified TLE.

   Add the '-v' (verbose) switch,  and it also writes out the orbital
   period and perigee/apogee distances.  Eventually,  I'll probably
   set it up to dump other data that are not immediately obvious
   just by looking at the TLEs... */

int main( const int argc, const char **argv)
{
   FILE *ifile;
   const char *filename;
   char line0[100], line1[100], line2[100];
   int i, verbose = 0;
   const char *norad = NULL, *intl = NULL;
   bool legend_shown = false;

   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v':
               verbose = 1;
               break;
            case 'n':
               norad = argv[i] + 2;
               if( !*norad && i < argc - 1)
                  norad = argv[++i];
               printf( "Looking for NORAD %s\n", norad);
               break;
            case 'i':
               intl = argv[i] + 2;
               if( !*intl && i < argc - 1)
                  intl = argv[++i];
               printf( "Looking for international ID %s\n", intl);
               break;
            default:
               printf( "'%s': unrecognized option\n", argv[i]);
               return( -1);
               break;
            }
   filename = (argc == 1 ? "all_tle.txt" : argv[1]);
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
                        && (!intl || !memcmp( line1 + 9, intl, 5))
                   && *line2 == '2')
         {
         tle_t tle;

         if( parse_elements( line1, line2, &tle) >= 0)
            {
            char obuff[200];
            double params[N_SGP4_PARAMS], c2;

            if( verbose && !legend_shown)
               {
               legend_shown = true;
               printf(
  "1 NoradU COSPAR   Epoch.epoch     dn/dt/2  d2n/dt2/6 BSTAR    T El# C\n"
  "2 NoradU Inclina RAAscNode Eccent  ArgPeri MeanAno  MeanMotion Rev# C\n");

               }
            SGP4_init( params, &tle);
            c2 = params[0];
            if( c2 && tle.xno)
               tle.bstar = tle.xndt2o / (tle.xno * c2 * 1.5);
            write_elements_in_tle_format( obuff, &tle);
            if( strlen( line0) < 60)
               printf( "%s", line0);
            printf( "%s", obuff);
            if( verbose)
               {
               const double a1 = pow(xke / tle.xno, two_thirds);  /* in Earth radii */

               printf( "Inclination: %8.4f     ", tle.xincl * 180. / PI);
               printf( "   Perigee: %.4f km\n",
                    (a1 * (1. - tle.eo) - 1.) * earth_radius_in_km);

               printf( "Asc node:    %8.4f     ", tle.xnodeo * 180. / PI);
               printf( "   Apogee: %.4f km\n",
                    (a1 * (1. + tle.eo) - 1.) * earth_radius_in_km);

               printf( "Arg perigee: %8.4f     ", tle.omegao * 180. / PI);
               printf( "   Orbital period: %.4f min\n",
                    2. * pi / tle.xno);

               printf( "Mean anom:   %8.4f     ", tle.xmo * 180. / PI);
               printf( "   Epoch: JD %.5f\n", tle.epoch);
               }
            }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( 0);
}
