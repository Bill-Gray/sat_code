/* Copyright (C) 2018, Project Pluto.  See LICENSE.

fix_tles : change TLE data and replace with correct checksums

   If one alters a TLE,  one will normally cause a change in
the checksum data at the end of the line.  This program reads
in successive lines from an input file;  if a line and the
preceding line make a TLE,  the checksum is computed for
both lines and is suitably reset.

   I've had instances where I realized that either the COSPAR
or NORAD designation was incorrectly set,  and added options
that let you specify which designation should be used for all
TLEs in the output.  (This is specific to my use case : I usually
compute many TLEs for a particular object,  which each TLE
being fitted to give good state vectors over one day.)

   At some point,  I may need to similarly batch-correct bulletin
numbers or something of that ilk,  but at present,  only the
designations (and checksums) can be batch-corrected with this
program.                                  */

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "norad.h"

/* After being modified,  the checksum byte in TLEs must be recomputed : */

static void set_checksum( char *line)
{
   const int csum =  tle_checksum( line);

   assert( csum >= 0);
   line[68] += csum;
   if( line[68] > '9')
      line[68] -= 10;
}

static void usage( void)
{
   fprintf( stderr,
       "'fix_tles' reads in TLEs and outputs TLEs with corrected checksums.\n"
       "Options are :\n"
       "   -i YYNNNA     Replace COSPAR designation\n"
       "   -n NNNNN      Replace five-digit NORAD designation\n"
       "   -f filename   Specify input file (default = stdin)\n"
       "   -o filename   Specify output file (default = stdout)\n" );
   exit( -1);
}

int main( const int argc, const char **argv)
{
   int i, line_no = 0;
   const char *norad_desig = NULL;
   char intl_desig[10];
   FILE *ifile = stdin, *ofile = stdout;
   char line1[200], line2[200];

   *intl_desig = '\0';
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *arg = argv[i] + 2;

         if( !*arg && i < argc - 1)
            arg = argv[i + 1];
         switch( argv[i][1])
            {
            case 'i':
               assert( strlen( arg) < 9);
               snprintf( intl_desig, sizeof( intl_desig), "%-9s", arg);
               break;
            case 'n':
               norad_desig = arg;
               break;
            case 'o':
               ofile = fopen( arg, "wb");
               if( !ofile)
                  {
                  fprintf( stderr, "'%s' not opened\n", arg);
                  usage( );
                  }
               break;
            case 'f':
               ifile = fopen( arg, "rb");
               if( !ifile)
                  {
                  fprintf( stderr, "'%s' not opened\n", arg);
                  usage( );
                  }
               break;
            default:
               fprintf( stderr, "Unrecognized option '%s'\n", argv[i]);
               usage( );
               break;
            }
         }
   *line1 = '\0';
   while( fgets( line2, sizeof( line2), ifile))
      {
      tle_t unused_tle;

      if( parse_elements( line1, line2, &unused_tle) >= 0)
         {
         if( norad_desig)
            {
            memcpy( line1 + 2, norad_desig, 5);
            memcpy( line2 + 2, norad_desig, 5);
            }
         if( *intl_desig)
            memcpy( line1 + 9, intl_desig, 8);
         set_checksum( line1);
         set_checksum( line2);
         }
      if( line_no)
         fputs( line1, ofile);
      line_no++;
      strcpy( line1, line2);
      }
   if( line_no)
      fputs( line1, ofile);
   return( 0);
}
