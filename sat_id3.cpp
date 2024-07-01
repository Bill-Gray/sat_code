/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __has_include
   #if __has_include(<cgi_func.h>)
       #include "cgi_func.h"
   #else
       #error   \
         'cgi_func.h' not found.  This project depends on the 'lunar'\
         library.  See www.github.com/Bill-Gray/lunar .\
         Clone that repository,  'make'  and 'make install' it.
   #endif
#else
   #include "cgi_func.h"
#endif
#include "watdefs.h"
#include "stringex.h"

/* This is the code behind

https://www.projectpluto.com/sat_img.htm

   It takes the input data,  creates a temporary file containing a
single 'image field data' line,  then uses the Sat_ID code (see sat_id.cpp)
to figure out which artsats are in that field.  */

int sat_id_main( const int argc, const char **argv);

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   char buff[80];
   char field[30];
   char search_radius[10], date[50], latitude[30], longitude[30];
   char altitude[10], ra[20], dec[20];
   const int argc = 3;
   const char *output_file_name = "field.txt";
   FILE *ofile;
   FILE *lock_file = fopen( "lock.txt", "w");
   size_t i;
   int cgi_status;
#ifndef _WIN32
   extern char **environ;

   avoid_runaway_process( 15);
#endif         /* _WIN32 */
   setbuf( lock_file, NULL);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argc);
   INTENTIONALLY_UNUSED_PARAMETER( unused_argv);
   printf( "Content-type: text/html\n\n");
   printf( "<html> <body> <pre>\n");
   if( !lock_file)
      {
      printf( "<p> Server is busy.  Try again in a minute or two. </p>");
      printf( "<p> Your orbit is very important to us! </p>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
#ifndef _WIN32
   for( i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
#endif
   cgi_status = initialize_cgi_reading( );
   fprintf( lock_file, "CGI status %d\n", cgi_status);
   if( cgi_status <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", cgi_status);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   *search_radius = *date = *latitude = *longitude = *altitude = '\0';
   *ra = *dec = '\0';
   while( !get_cgi_data( field, buff, NULL, sizeof( buff)))
      {
      fprintf( lock_file, "Field '%s'\n", field);
      if( !strcmp( field, "radius"))
         {
         strlcpy( search_radius, buff, sizeof( search_radius));
         if( !atof( search_radius))
            printf( "Search radius must be non-zero\n");
         }
      else if( !strcmp( field, "time"))
         strlcpy( date, buff, sizeof( date));
      else if( !strcmp( field, "lat"))
         strlcpy( latitude, buff, sizeof( latitude));
      else if( !strcmp( field, "lon"))
         strlcpy( longitude, buff, sizeof( longitude));
      else if( !strcmp( field, "alt"))
         strlcpy( altitude, buff, sizeof( altitude));
      else if( !strcmp( field, "ra"))
         strlcpy( ra, buff, sizeof( ra));
      else if( !strcmp( field, "dec"))
         strlcpy( dec, buff, sizeof( dec));
      }
   ofile = fopen( output_file_name, "w");
   fprintf( ofile, "COD XXX\n");
   fprintf( ofile, "COM Long. %s, Lat. %s, Alt. %s, unspecified\n",
                  longitude, latitude, altitude);
   fprintf( ofile, "Field,%s,%s,%s,XXX\n", date, ra, dec);
   fclose( ofile);
   argv[0] = "sat_id";
   argv[1] = output_file_name;
   argv[2] = "-t../../tles/tle_list.txt";
   argv[3] = NULL;
   for( i = 0; argv[i]; i++)
      fprintf( lock_file, "arg %d: '%s'\n", (int)i, argv[i]);
   sat_id_main( argc, argv);
   fprintf( lock_file, "sat_id_main called\n");
   printf( "On-line artsats-in-field-finder compiled"
                            __DATE__ " " __TIME__ " UTC-5h\n");
   printf( "See <a href='https://www.github.com/Bill-Gray/sat_code'>"
               "https://www.github.com/Bill-Gray/sat_code</a> for source code\n");
   printf( "</pre> </body> </html>");
   return( 0);
}
