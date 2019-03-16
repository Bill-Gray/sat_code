/* Code to check for the existence of certain artsats in Space-Track's
master TLE list.  Occasionally,  they've dropped objects and I didn't
realize it.  The objects ended up on NEOCP and I didn't ID them as
quickly as might be desired,  because I assumed they must be "new".
This should warn me if certain artsats get dropped from 'all_tle.txt'.
*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

int main( const int argc, const char **argv)
{
   static const char *sats[] = {
          "00041B ", "00045A ", "00045B ", "02048A ", "07004D ",
          "07004E ", "07004G ", "15011A ", "15011B ",
          "15011C ", "15011D ", "67040F ", "69046F ",
          "77093E ", "83020A ", "83020D ", "99040B ",
          "99066A ", "99066B ",
          NULL };
   FILE *ifile = fopen( argc == 1 ? "all_tle.txt" : argv[1], "rb");
   char buff[100];
   size_t i;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff == '1' && buff[1] == ' ')
         for( i = 0; sats[i]; i++)
            if( sats[i][0] == buff[9] && !memcmp( sats[i], buff + 9, 7))
               sats[i] = "";
   fclose( ifile);
   printf( "This should list 2002-048A INTEGRAL and 1983-020A and D,  and nothing else.\n");
   printf( "(At times,  one or more of those may be listed,  too,  but you shouldn't\n");
   printf( "rely on that happening.)\n");
   for( i = 0; sats[i]; i++)
      if( sats[i][0])
         printf( "%s\n", sats[i]);
   return( 0);
}
