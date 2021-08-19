/* Copyright (C) 2021, Project Pluto.  See LICENSE.  */
/* Code to test Alpha-5 and 'extended' Alpha-5 (Super-5) designation
schemes.  These enable us to have NORAD numbers up to 34^5-1
= 45435423 while still fitting requirements of the current TLE
format.  See comments about the Alpha-5 and Super-5 schemes in
'get_el.cpp'.

Run without command-line arguments,  this tries out all possible
Super-5 designations with a "round-trip" test : encode,  then
decode to ensure we get the original number.  A couple past
34^5 are run to make sure that they fail as expected.  In this
scheme,  the lexical order of Super-5 designations should
mismatch the numerical order at several points,  and this test
code verifies that. */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "norad.h"

int main( const int argc, const char **argv)
{
   char buff[160], *line1 = buff, line2[80];
   int n;
   tle_t tle;

   strcpy( line1,
      "1 00005U 58002B   21124.86416743 -.00000116  00000-0 -00000-0 0  9998");
   strcpy( line2,
      "2 00005  34.2496 318.0626 1847159 358.7427   0.8906 10.84838792240168");
   parse_elements( line1, line2, &tle);
   if( argc > 1)
      {
      tle.norad_number = atoi( argv[1]);
      write_elements_in_tle_format( buff, &tle);
      printf( "Super-5 version : '%.5s'\n", buff + 2);
      return( 0);
      }

   printf( "Testing all Super-5 designations.  This should take a minute or two.\n");
   for( n = 0; n < 34 * 34 * 34 * 34 * 34 + 2; n++)
      {
      char prev_des[6];

      memcpy( prev_des, buff + 2, 5);
      tle.norad_number = n;
      write_elements_in_tle_format( buff, &tle);
      if( memcmp( prev_des, buff + 2, 5) >= 0)
         printf( "Disorder at %d\n%s\n", n, buff);
      parse_elements( line1, line2, &tle);
      if( tle.norad_number != n)
         printf( "Got back NORAD %d (should be %d)\n%s\n",
                     tle.norad_number, n, buff);
      if( !(n % 1000000))        /* progress indicator */
         {
         printf( ".");
         fflush( stdout);
         }
      }
   return( 0);
}
