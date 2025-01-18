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

   As of 2024 Aug 31,  the program can also be used for updating the
Space-Track TLEs in a slightly more cautious manner.  If you have
downloaded new TLEs as the (default) ALL_TLE.TXT,  and your "usual"
TLEs are at all_tle.txt,  then

./dropouts ALL_TLE.TXT 25000 all_tle.txt

   will check to see if ALL_TLE.TXT actually has 25000 or more TLEs in it.
If it does,  the download is presumed to have succeeded;  all_tle.txt is
unlinked and replaced with ALL_TLE.TXT.  If it fails,  we leave both files
undisturbed.

   This should help in the increasingly frequent situations where new
TLE files are downloaded and then have only an error message,  or a
drastically reduced number of TLEs.        */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#ifdef _WIN32
   #include <io.h>
#else
   #include <unistd.h>
#endif

#define VT_NORMAL        "\033[0m"
#define VT_REVERSE       "\033[7m"

int main( const int argc, const char **argv)
{
   static const char *sats[] = {
          "00041A ! Cluster II-FM7",
          "00045A ! Cluster II-FM5",
          "00045B ! Cluster II-FM8",
          "02048A ! INTEGRAL",
          "07004A ! THEMIS-A",
          "07004D ! THEMIS-D",
          "07004E ! THEMIS-E",
          "09017B ! Atlas 5 Centaur R/B",
          "09068B ! Delta 4 R/B",
          "12003B ! Delta 4 R/B",
          "12011B ! Breeze-M R/B",
          "13024B ! WGS-5 R/B",
          "13026B ! Breeze-M R/B",
          "15005B ! Inmarsat 5F2 booster",
          "15011A ! MMS 1",
          "15011B ! MMS 2",
          "15011C ! MMS 3",
          "15011D ! MMS 4",
          "15019C ! Yuanzheng-1 Y1",
          "15042B ! Breeze-M R/B",
          "16041A ! MUOS 5",
          "18038A   TESS",
          "22110B ! Ariane 5 R/B",
          "22146B ! Falcon 9 R/B",
          "22134B ! Falcon 9 R/B",
          "24048E ! DRO R/B",
          "24059B ! Falcon 9 R/B",
          "24127B ! Falcon 9 R/B",
          "24214B ! Falcon 9 R/B",
          "24233A ! Proba-3",
          "24233B ! Proba-3 booster",
          "63039A ! Vela 1A",
          "64040B ! Vela 2B",
          "65058A ! Vela 3A",
          "65058B   Vela 6",
          "67040A   Vela 4A",
          "67040F ! Titan 3C transtage booster",
          "69046F ! Titan 3C transtage booster",
          "70027C ! Vela 6 booster",
          "72073A   IMP-7",
          "76023C ! SOLRAD-11A",
          "76023H ! SOLRAD-11 debris",
          "77093E ! SL-6 R/B(2)",
          "83020A ! ASTRON",
          "83020D ! ASTRON booster",
          "92044A   GEOTAIL",
          "95062A ! ISO",
          "95062C ! ISO debris",
          "97075B ! Equator S",
          "99040B ! Chandra X-Ray Observatory",
          "99040D ! IUS (for CXO)",
          "99066A ! XMM/Newton",
          "99066B ! XMM/Newton booster",
          NULL };
   FILE *ifile = fopen( argc == 1 ? "all_tle.txt" : argv[1], "rb");
   char buff[100];
   size_t i;
   int trouble_found = 0, n_found = 0;

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff == '1' && buff[1] == ' ' && buff[7] == 'U')
         {
         size_t len = strlen( buff);

         while( len && buff[len - 1] < ' ')
            len--;
         if( 69 == len)
            {
            n_found++;
            for( i = 0; sats[i]; i++)
               if( sats[i][0] == buff[9] && !memcmp( sats[i], buff + 9, 7))
                  sats[i] = "";
            }
         }
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
   printf( "%d found\n", n_found);
   if( 4 == argc && n_found > atoi( argv[2]))
      {
      printf( "Replacing '%s' with '%s'\n", argv[3], argv[1]);
#ifdef _WIN32
      _unlink( argv[3]);
#else
      unlink( argv[3]);
#endif
      rename( argv[1], argv[3]);
      }
   return( 0);
}
