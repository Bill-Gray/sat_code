/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"

/* Code to extract TLEs for a particular date (MJD).  It looks
through all TLEs listed in 'tle_list.txt' and outputs those matching
the MJD given on the command line.  */

char path[100];
int verbose;

static void extract_tle_for_date( const char *fname, const double mjd)
{
   char prev_line[100], buff[100];
   FILE *ifile = fopen( fname, "rb");
   double step = 0.;

   assert( strlen( fname) < 20);
   assert( strlen( path) < 60);
   strcpy( buff, path);
   strcat( buff, fname);
   ifile = fopen( buff, "rb");
   if( !ifile)
      {
      printf( "'%s' not opened\n", fname);
      exit( -1);
      }
   if( verbose)
      printf( "Looking at '%s' for %f\n", fname, mjd);
   *prev_line = '\0';
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "# MJD ", 6) && atof( buff + 6) <= mjd
                     && atof( buff + 6) + step > mjd)
         {
         size_t i;

         printf( "%s", prev_line);    /* show 'worst residual' data */
         for( i = 0; i < 3 && fgets( buff, sizeof( buff), ifile); i++)
            printf( "%s", buff);   /* obj name,  lines 1 and 2 of TLE */
         fclose( ifile);
         return;
         }
      else if( !memcmp( buff, "# Ephem range: ", 15))
         {
         double mjd1, mjd2;

         if( sscanf( buff + 15, "%lf %lf %lf", &mjd1, &mjd2, &step) != 3)
            {
            printf( "Ephem step problem in '%s'\n'%s'\n",
                     fname, buff);
            exit( -2);
            }
         if( mjd < mjd1 || mjd > mjd2)
            {
            if( verbose)
               printf( "'%s': outside range\n", fname);
            fclose( ifile);
            return;
            }
         }
      else if( !memcmp( buff, "# Include ", 10))
         {
         int i = 0;
         char *filename = buff + 10;

         while( buff[i] >= ' ')
            i++;
         buff[i] = '\0';
         if( memcmp( filename, "classfd", 7) && memcmp( filename, "inttles", 7)
               && memcmp( filename, "all_tle", 7)
               && memcmp( filename, "old_tle", 7))
            extract_tle_for_date( filename, mjd);
         }
      strcpy( prev_line, buff);
      }
   fclose( ifile);
}

static void err_exit( void)
{
   printf( "tle_date (MJD)\n");
   exit( -1);
}

#ifdef ON_LINE_VERSION
int dummy_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   const double jan_1_1970 = 2440587.5;
   const double curr_t = jan_1_1970 + (double)time( NULL) / seconds_per_day;
   double mjd;

   if( argc < 2)
      err_exit( );
   mjd = atof( argv[1]);
   mjd = get_time_from_string( curr_t, argv[1], FULL_CTIME_YMD, NULL)
                        - 2400000.5;
   if( mjd < 40000 || mjd > 65000)
      err_exit( );
   if( argc > 2)
      strcpy( path, argv[2]);
   if( argc > 3)
      verbose = 1;
   extract_tle_for_date( "tle_list.txt", mjd);
   return( 0);
}

#ifdef ON_LINE_VERSION
#include <cgi_func.h>

int main( void)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;
   char *buff = (char *)malloc( max_buff_size);
   char field[30], date_text[80];
   FILE *lock_file = fopen( "lock.txt", "w");
   extern char **environ;
   int cgi_status;

   avoid_runaway_process( 15);
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "<p> Server is busy.  Try again in a minute or two. </p>");
      printf( "<p> Your TLEs are very important to us! </p>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
   cgi_status = initialize_cgi_reading( );
   strcpy( date_text, "now");
   fprintf( lock_file, "CGI status %d\n", cgi_status);
   if( cgi_status <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", cgi_status);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, max_buff_size))
      {
      if( !strcmp( field, "date") && strlen( buff) < 70)
         {
         strncpy( date_text, buff, sizeof( date_text));
         date_text[sizeof( date_text) - 1] = '\0';
         }
      }
   free( buff);
   argv[0] = "tle_date";
   argv[1] = date_text;
   argv[2] = "/home/projectp/public_html/tles/";
   argv[3] = NULL;
   dummy_main( 3, argv);
   printf( "</pre> </body> </html>");
   fclose( lock_file);
   return( 0);
}
#endif
