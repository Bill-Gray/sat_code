/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int is_position_line( const char *buff)
{
   int rval = 0;

   if( strlen( buff) > 73 && buff[73] < ' ' &&
            buff[29] == '.' && buff[46] == '.' && buff[63] == '.')
      rval = 1;
   return( rval);
}

int main( const int argc, const char **argv)
{
   FILE *ifile1 = fopen( argv[1], "rb");
   FILE *ifile2 = fopen( argv[2], "rb");
   char buff1[80], buff2[80];
   int line = 0, n_valid = 0, worst_line = -1;
   double max_diff2 = 0.;

   if( !ifile1)
      printf( "%s not opened\n", argv[1]);
   if( !ifile2)
      printf( "%s not opened\n", argv[2]);
   if( argc < 3 || !ifile1 || !ifile2)
      {
      printf( "out_comp needs two files of output from test_sat to compare.\n");
      exit( -1);
      }

   while( fgets( buff1, sizeof( buff1), ifile1) &&
          fgets( buff2, sizeof( buff2), ifile2))
      {
      line++;
      if( is_position_line( buff1) && is_position_line( buff2))
         {
         double diff2 = 0., delta;
         int i;

         n_valid++;
         for( i = 22; i < 60; i += 17)
            {
            delta = atof( buff1 + i) - atof( buff2 + i);
            diff2 += delta * delta;
            }
         if( diff2 > max_diff2)
            {
            max_diff2 = diff2;
            worst_line = line;
            }
         }
      }
   fclose( ifile1);
   fclose( ifile2);
   printf( "%d lines read in; %d had positions\n", line, n_valid);
   printf( "Max difference: %.8f km at line %d\n",
                                  sqrt( max_diff2), worst_line);

   return( 0);
}
