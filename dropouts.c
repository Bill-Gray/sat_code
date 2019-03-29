/* Code to check for the existence of certain artsats in Space-Track's
master TLE list.  Occasionally,  they've dropped objects and I didn't
realize it.  The objects ended up on NEOCP and I didn't ID them as
quickly as might be desired,  because I assumed they must be "new".
This should warn me if certain artsats get dropped from 'all_tle.txt'.

   The absence of certain artsats is essentially routine.  But for some
objects (marked with an !),  Space-Track is our only source of TLEs.
(Or at least,  I've been relying on them.  I _could_ generate TLES for
CXO,  for example,  based on _Horizons_ ephems.  But I want this
program to squawk loudly if Space-Track stops supplying CXO TLEs.)
*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#define VT_NORMAL        "\033[0m"
#define VT_REVERSE       "\033[7m"

int main( const int argc, const char **argv)
{
   static const char *sats[] = {
          "00041B ! Cluster II-FM7",
          "00045A ! Cluster II-FM5",
          "00045B ! Cluster II-FM8",
          "02048A ! INTEGRAL",
          "07004A   THEMIS-A",
          "07004D   THEMIS-D",
          "07004E   THEMIS-E",
          "07004G   THEMIS booster",
          "15011A ! MMS 1",
          "15011B ! MMS 2",
          "15011C ! MMS 3",
          "15011D ! MMS 4",
          "15019C   YZ-1 booster",
          "67040F ! Titan 3C transtage booster",
          "69046F   Titan 3C transtage booster",
          "77093E   SL-6 R/B(2)",
          "83020A   ASTRON",
          "83020D   ASTRON booster",
          "99040B ! Chandra X-Ray Observatory",
          "99066A ! XMM/Newton",
          "99066B ! XMM/Newton booster",
          NULL };
   FILE *ifile = fopen( argc == 1 ? "all_tle.txt" : argv[1], "rb");
   char buff[100];
   size_t i;
   int trouble_found = 0;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff == '1' && buff[1] == ' ')
         for( i = 0; sats[i]; i++)
            if( sats[i][0] == buff[9] && !memcmp( sats[i], buff + 9, 7))
               sats[i] = "";
   fclose( ifile);
   printf( "This will list high-flying artsats for which TLEs are not provided :\n");
   for( i = 0; sats[i]; i++)
      if( sats[i][0])
         {
         printf( "%s\n", sats[i]);
         if( sats[i][7] == '!')
            {
            trouble_found = 1;
            printf( VT_REVERSE);
            printf( "DANGER!!! We do NOT have an independent source of TLEs\n");
            printf( "for this object.  Please report to "
                        "pluto\x40\x70roject\x70luto\x2e\x63om.\n");
                   /* Above is (slightly) obfuscated address to foil spambots */
            printf( "This needs to be fixed.\n");
            printf( VT_NORMAL);
            }
         }
   if( !trouble_found)
      printf( "Any missing objects are covered by other sources.  Nothing\n"
              "to worry about here.\n");
   return( 0);
}
