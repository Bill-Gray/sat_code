/* Code to read 'tle_list.txt' and insert # Range: lines for
files that lack them,  and (for files containing only one
object) NORAD/COSPAR identifiers.   The former can be used to
speed up some programs (by skipping 'included' files that
won't cover a desired time span).  The latter can be used in
an ephemeris program to quickly find a desired object.

Important notes :

   -- tle_list.txt will require occasional updating,  of course,
as new objects are added and archival Space-Track and other TLE
sets are added.  My hope is to run this once and then add
# Range: and # ID: lines as needed.

   -- THEMIS-A, D, and E are oddball cases.  They maneuver;  I
get state vectors from UC Berkeley and compute TLEs from those
(see 'up_them' in the 'tles' repository).  Those are updated
roughly weekly.  I don't want to have to update 'tle_list' each
time for that,  so I've pushed the day of reckoning for those
off by a year.  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "norad.h"

static void get_range_info( const char *filename,
               const int range_already_set)
{
   char buff[200], *tptr, id[100];
   int n_ids_found = 0;
   const int is_themis = !memcmp( filename, "07004", 5)
                           && filename[5] < 'g';
   FILE *ifile;

   strcpy( buff, filename);
   tptr = strchr( buff, '\n');
   assert( tptr);
   *tptr = '\0';
   ifile = fopen( buff, "rb");
   assert( ifile);
   *id = '\0';
   while( 2 > n_ids_found && fgets( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "# Ephem range:", 14) && !range_already_set)
         {
         double mjd[2];
         size_t i;
         const int n_read = sscanf( buff + 14, "%lf %lf", mjd, mjd + 1);

         assert( n_read == 2);
         if( is_themis)          /* THEMIS sats get updated/extended */
            mjd[1] += 365.;  /* regularly;  'pad' the end date accordingly */
         printf( "# Range:");
         for( i = 0; i < 2; i++)
            {
            full_ctime( buff, 2400000.5 + mjd[i],
                  FULL_CTIME_YMD | FULL_CTIME_DATE_ONLY
                  | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MONTHS_AS_DIGITS);
            buff[4] = buff[7] = '-';
            printf( " %s", buff);
            }
         printf( "\n");
         }
      else if( *buff == '1' && !tle_checksum( buff))
         {
         buff[7] = ' ';
         buff[17] = '\0';
         if( strcmp( id, buff + 2))
            n_ids_found++;
         strcpy( id, buff + 2);
         }
      }
   if( n_ids_found == 1)
      printf( "# ID: %s\n", id);
   fclose( ifile);
}

int main( const int unused_argc, const char **unused_argv)
{
   char prev_line[100], line[100];
   FILE *ifile = fopen( "tle_list.txt", "rb");

   INTENTIONALLY_UNUSED_PARAMETER( unused_argc);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argv);
   assert( ifile);
   *prev_line = '\0';
   while( fgets( line, sizeof( line), ifile))
      {
      if( !memcmp( line, "# Include ", 10)
                     && memcmp( line + 10, "old_tles", 8))
         get_range_info( line + 10, !memcmp( prev_line, "# Range:", 8));
      printf( "%s", line);
      strcpy( prev_line, line);
      }
   fclose( ifile);
   return( 0);
}
