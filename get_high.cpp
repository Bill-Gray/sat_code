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
   fprintf( ofile, "# Added %s\n", ctime( &t0));
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
