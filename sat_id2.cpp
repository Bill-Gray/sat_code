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

int sat_id_main( const int argc, const char **argv);

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   const size_t max_buff_size = 400000;       /* room for 5000 obs */
   char *buff = (char *)malloc( max_buff_size), *tptr;
   char field[30];
   const char *temp_obs_filename = "sat_obs.txt";
   double search_radius = 4.;     /* look 4 degrees for matches */
   double motion_cutoff = 60.;  /* up to 60" discrepancy OK */
   double low_speed_cutoff = 0.001;  /* anything slower than this is almost */
   int argc = 0;                     /* certainly not an artsat */
   FILE *lock_file = fopen( "lock.txt", "w");
   size_t bytes_written = 0, i;
   int cgi_status, show_summary = 0;
   extern int verbose;
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
   while( !get_cgi_data( field, buff, NULL, max_buff_size))
      {
      fprintf( lock_file, "Field '%s'\n", field);
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

            fprintf( lock_file, "File opened : %p\n", (void *)ofile);
            if( !ofile)
               {
               printf( "<p> Couldn't open %s : %s </p>\n", temp_obs_filename, strerror( errno));
               fprintf( lock_file, "Couldn't open %s: %s\n", temp_obs_filename, strerror( errno));
               return( -1);
               }
            bytes_written += fwrite( buff, 1, strlen( buff), ofile);
            fclose( ofile);
            }
         }
      if( !strcmp( field, "radius"))
         {
         const char *verbosity = strchr( buff, 'v');

         search_radius = atof( buff);
         if( verbosity)
            verbose = atoi( verbosity + 1) + 1;
         }
      if( !strcmp( field, "motion"))
         motion_cutoff = atof( buff);
      if( !strcmp( field, "low_speed"))
         low_speed_cutoff = atof( buff);
      if( !strcmp( field, "summary"))
         show_summary = 1;
      }
   fprintf( lock_file, "Fields read\n");
// printf( "<p>Fields read</p>\n");
   if( verbose)
      printf( "Searching to %f degrees;  %u bytes read from input\n",
                     search_radius, (unsigned)bytes_written);
   argv[argc++] = "sat_id";
   argv[argc++] = temp_obs_filename;
   argv[argc++] = "-t../../tles/tle_list.txt";
   snprintf_err( field, sizeof( field), "-r%.2f", search_radius);
   argv[argc++] = field;
   snprintf_err( buff, max_buff_size, "-y%f", motion_cutoff);
   argv[argc++] = buff;
   tptr = buff + strlen( buff) + 1;
   snprintf( tptr, 15, "-z%f", low_speed_cutoff);
   argv[argc++] = tptr;
   if( show_summary)
      argv[argc++] = "-u";
   argv[argc] = NULL;
   for( i = 0; argv[i]; i++)
      fprintf( lock_file, "arg %d: '%s'\n", (int)i, argv[i]);
   sat_id_main( argc, argv);
   fprintf( lock_file, "sat_id_main called\n");
   free( buff);
   printf( "On-line Sat_ID compiled " __DATE__ " " __TIME__ " UTC-5h\n");
   printf( "See <a href='https://www.github.com/Bill-Gray/sat_code'>"
               "https://www.github.com/Bill-Gray/sat_code</a> for source code\n");
   printf( "</pre> </body> </html>");
   return( 0);
}
