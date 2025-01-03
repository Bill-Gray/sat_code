/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/* Code to extract elements for high-flying artsats. Give */
/* it the name of the input file of TLEs and a cutoff of  */
/* the mean motion,  and only TLEs with a lower motion    */
/* will be output.                                        */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "norad.h"

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen((argc > 1 ? argv[1] : "all_tle.txt"), "rb");
   FILE *ofile;
   char line0[200], line1[200], line2[200];
   const double cutoff = (argc > 2 ? atof( argv[2]) : .6);
   const time_t t0 = time( NULL);

   if( !ifile)
      perror( "Input file not opened");
   ofile = (argc > 3 ? fopen( argv[3], "a") : stdout);
   if( !ofile)
      perror( "Output file not opened");
   if( !ifile || !ofile)
      return( -1);
   *line0 = *line1 = '\0';
   fprintf( ofile, "# Added %.24s UTC\n", asctime( gmtime( &t0)));
   while( fgets( line2, sizeof( line2), ifile))
      {
      if( *line2 == '2' && *line1 == '1'
               && !tle_checksum( line1) && !tle_checksum( line2)
               && atof( line2 + 52) < cutoff)
         fprintf( ofile, "%s%s%s", line0, line1, line2);
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( 0);
}
