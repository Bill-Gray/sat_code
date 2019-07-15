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

int sat_id_main( const int argc, const char **argv);

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[20];
   const size_t max_buff_size = 40000;       /* room for 500 obs */
   char *buff = (char *)malloc( max_buff_size);
   char boundary[100], field[30];
   const char *temp_obs_filename = "sat_obs.txt";
   double search_radius = 4.;     /* look 2 degrees for matches */
   double motion_cutoff = 0.015;  /* up to .015'/s motion discrepancy OK */
   double low_speed_cutoff = 0.003;  /* anything slower than this is almost */
   const int argc = 6;               /* certainly not an artsat */
   FILE *lock_file = fopen( "lock.txt", "w");
   size_t bytes_written = 0;
   extern int verbose;
#ifndef _WIN32
   extern char **environ;

   avoid_runaway_process( 15);
#endif         /* _WIN32 */
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
   for( size_t i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
#endif
   if( !fgets( boundary, sizeof( boundary), stdin))
      {
      printf( "<p><b> No info read from stdin</b></p>");
      printf( "<p>This isn't supposed to happen.</p>");
      return( 0);
      }
   while( get_multipart_form_data( boundary, field, buff, NULL, max_buff_size) >= 0)
      {
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

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
      }
   if( verbose)
      printf( "Searching to %f degrees;  %u bytes read from input\n",
                     search_radius, (unsigned)bytes_written);
   argv[0] = "sat_id";
   argv[1] = temp_obs_filename;
   argv[2] = "-t../../tles/tle_list.txt";
   sprintf( field, "-r%.2f", search_radius);
   argv[3] = field;
   sprintf( buff, "-y%f", motion_cutoff);
   argv[4] = buff;
   sprintf( boundary, "-z%f", low_speed_cutoff);
   argv[5] = boundary;
   argv[6] = NULL;
   sat_id_main( argc, argv);
   free( buff);
   printf( "</pre> </body> </html>");
   return( 0);
}
