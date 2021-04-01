#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sat_util.h"

char *fgets_trimmed( char *buff, const int buffsize, FILE *ifile)
{
   char *rval = fgets( buff, buffsize, ifile);

   if( rval)
      {
      size_t i = 0;

      while( rval[i] != 10 && rval[i] != 13 && rval[i])
         i++;
      while( i && rval[i - 1] == ' ')
         i--;        /* drop trailing spaces */
      rval[i] = '\0';
      }
   return( rval);
}

#if !defined( _WIN32)
void make_config_dir_name( char *oname, const char *iname)
{
#ifdef ON_LINE_VERSION
   strcpy( oname, getenv( "DOCUMENT_ROOT"));
   strcat( oname, "/");
#else
   const char *home_dir = getenv( "HOME");

   if( home_dir)
      {
      strcpy( oname, home_dir);
      strcat( oname, "/.find_orb/");
      }
   else
      *oname = '\0';
#endif
   strcat( oname, iname);
}

FILE *local_then_config_fopen( const char *filename, const char *permits)
{
   FILE *rval = fopen( filename, permits);

   if( !rval)
      {
      char ext_filename[255];

      make_config_dir_name( ext_filename, filename);
      rval = fopen( ext_filename, "rb");
      }
   return( rval);
}
#else
FILE *local_then_config_fopen( const char *filename, const char *permits)
{
   return( fopen( filename, permits));
}
#endif
