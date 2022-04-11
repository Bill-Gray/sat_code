/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "norad.h"
#include "norad_in.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

void greg_day_to_dmy( const long jd, int *day,
                  int *month, long *year)
{
   const long mar_1_year_0 = 1721120L;     /* JD 1721120 = 1.5 Mar 0 Greg */
   const long one_year = 365L;
   const long four_years = 4 * one_year + 1;
   const long century = 25 * four_years - 1L;  /* days in 100 'normal' yrs */
   const long quad_cent = century * 4 + 1;     /* days in 400 years */
   long days = jd - mar_1_year_0;
   long day_in_cycle = days % quad_cent;

   if( day_in_cycle < 0)
      day_in_cycle += quad_cent;
   *year = ((days - day_in_cycle) / quad_cent) * 400L;
   *year += (day_in_cycle / century) * 100L;
   if( day_in_cycle == quad_cent - 1)    /* extra leap day every 400 years */
      {
      *month = 2;
      *day = 29;
      return;
      }
   day_in_cycle %= century;
   *year += (day_in_cycle / four_years) * 4L;
   day_in_cycle %= four_years;
   *year +=  day_in_cycle / one_year;
   if( day_in_cycle == four_years - 1)    /* extra leap day every 4 years */
      {
      *month = 2;
      *day = 29;
      return;
      }
   day_in_cycle %= one_year;
   *month = 5 * (day_in_cycle / 153L);
   day_in_cycle %= 153L;
   *month += 2 * (day_in_cycle / 61L);
   day_in_cycle %= 61L;
   if( day_in_cycle >= 31)
      {
      (*month)++;
      day_in_cycle -= 31;
      }
   *month += 3;
   *day = day_in_cycle + 1;
   if( *month > 12)
      {
      *month -= 12;
      (*year)++;
      }
}

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
   int legend_shown = 0;

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
                        && (!intl || !memcmp( line1 + 9, intl, strlen( intl)))
                   && *line2 == '2')
         {
         tle_t tle;
         const int err_code = parse_elements( line1, line2, &tle);

         if( err_code >= 0)
            {
            char obuff[200];
            double params[N_SGP4_PARAMS], c2;

            if( verbose && !legend_shown)
               {
               legend_shown = 1;
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
               long year, ijd;
               int month, day;
               double frac;

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
               ijd = (long)( tle.epoch + 0.5);
               frac = tle.epoch + 0.5 - ijd;
               greg_day_to_dmy( ijd, &day, &month, &year);
               printf( "   Epoch: JD %.5f = %ld-%02d-%02d.%05d\n", tle.epoch,
                              year, month, day, (int)( frac * 100000.));
               printf( "NORAD number %7d         ", tle.norad_number);
               printf( "Semimajor axis : %4f km\n", a1 * earth_radius_in_km);
               }
            if( err_code)
               printf( "Checksum error %d\n", err_code);
            }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( 0);
}
