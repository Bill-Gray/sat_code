/* Copyright (C) 2018, Project Pluto.  See LICENSE.

Code to read one or more files of TLEs,  remove duplicates,
and sort them according to various possible criteria (eccentricity,
orbital period,  etc.)   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TLE struct tle
#define MAX_TLES 200000

/* As of late 2020,  NORAD numbers are stored in five digits.  There's
some possibility of that being bumped up to nine digits,  in which
case I'll have to revisit this code.         */

#define MAX_NORAD_NUMBER 100000

TLE
   {
   char name_line[80], line1[80], line2[80];
   };

int n_duplicates = 0, heavens_above_html_tles = 0;

int load_tles_from_file( FILE *ifile, TLE *tles, char *already_found)
{
   int rval = 0;
   char buff[80], prev_line[80];

   *prev_line = '\0';
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( *buff == '1' && heavens_above_html_tles &&
                            !memcmp( buff + 69, "<B></B>", 7))
         strcpy( buff + 69, "\n");
      if( *buff == '1' && strlen( buff) > 69 && buff[69] < ' ')
         {
         char buff2[80];
         const int norad_number = atoi( buff + 2);

         if( fgets( buff2, sizeof( buff2), ifile))
            if( *buff2 == '2' && strlen( buff2) > 69 && buff2[69] < ' ')
               {
               if( !already_found[norad_number])
                  {
                  if( heavens_above_html_tles)     /* can't use HA names */
                     *prev_line = '\0';
                  strcpy( tles[rval].name_line, prev_line);
                  strcpy( tles[rval].line1, buff);
                  strcpy( tles[rval].line2, buff2);
                  already_found[norad_number] = 1;
                  rval++;
                  *buff = '\0';
                  }
               else
                  n_duplicates++;
               }
         }
      strcpy( prev_line, buff);
      }
   return( rval);
}

FILE *test_fopen( const char *filename, const char *permits)
{
   FILE *rval = fopen( filename, permits);

   if( !rval)
      {
      printf( "%s not opened\n", filename);
      exit( -1);
      }
   return( rval);
}

void show_tle( FILE *ofile, const TLE *tle)
{
   fprintf( ofile, "%s", tle->name_line);
   fprintf( ofile, "%s", tle->line1);
   fprintf( ofile, "%s", tle->line2);
}

static double get_perigee( const TLE *tle)
{
   const double ecc = atof( tle->line2 + 26) * 1e-7;
   const double revs_per_day = atof( tle->line2 + 52);
   const double semimajor = pow( revs_per_day, -2. / 3.);

   return( semimajor * (1. - ecc));
}

static int compare_doubles( const char *buff1, const char *buff2)
{
   const double d1 = atof( buff1);
   const double d2 = atof( buff2);
   int rval;

   if( d1 < d2)
      rval = -1;
   else if( d1 > d2)
      rval = 1;
   else
      rval = 0;
   return( rval);
}

int tle_compare( const TLE *tle1, const TLE *tle2, const char sort_method)
{
   int i, rval = 0;

   switch( sort_method)
      {
      case 'n': case 'N':        /* sort by NORAD number */
         rval = atoi( tle1->line1 + 2) - atoi( tle2->line1 + 2);
         break;
      case 'c': case 'C':        /* sort by COSPAR (international) desig */
         if( tle1->line1[9] >= '5' && tle2->line1[9] < '5')
            rval = -1;
         else if( tle2->line1[9] >= '5' && tle1->line1[9] < '5')
            rval = 1;
         else        /* COSPAR IDs are from the same century */
            rval = memcmp( tle1->line1 + 9, tle2->line1 + 9, 8);
         break;
      case 'm': case 'M':        /* sort by mean motion */
         rval = compare_doubles( tle1->line2 + 52, tle2->line2 + 52);
         break;
      case 'e': case 'E':        /* sort by eccentricity */
         for( i = 26; !rval && i < 33; i++)
            rval = tle1->line2[i] - tle2->line2[i];
         break;
      case 'p': case 'P':        /* sort by epoch */
         for( i = 18; !rval && i < 32; i++)
            rval = tle1->line1[i] - tle2->line1[i];
         break;
      case 'i': case 'I':        /* sort by incl */
         rval = compare_doubles( tle1->line2 + 8, tle2->line2 + 8);
         break;
      case 'o': case 'O':        /* sort by ascending node */
         rval = compare_doubles( tle1->line2 + 17, tle2->line2 + 17);
         break;
      case 'q':
         {
         const double q1 = get_perigee( tle1);
         const double q2 = get_perigee( tle2);

         rval = ( q1 > q2 ? 1 : -1);
         }
         break;
      }
   if( sort_method >= 'A' && sort_method <= 'Z')
      rval = -rval;
   return( rval);
}

      /* MS Windows lacks the re-entrant qsort_r;  we have to use  */
      /* plain old non-re-entrant qsort and a global variable.     */
      /* BSD has it,  but with the wrong order of arguments.       */
#ifdef __linux
#define HAVE_REENTRANT_QSORT
#endif

#ifdef HAVE_REENTRANT_QSORT
int tle_compare_for_qsort_r( const void *a, const void *b, void *c)
{
   return( tle_compare( (const TLE *)a, (const TLE *)b, *(char *)c));
}
#else
static char comparison_method;

int tle_compare_for_qsort( const void *a, const void *b)
{
   return( tle_compare( (const TLE *)a, (const TLE *)b, comparison_method));
}
#endif


static void error_exit( void)
{
   printf( "'mergetle' will merge TLEs from one or more files.  Optionally,\n");
   printf( "the output can be sorted.  For example:\n\n");
   printf( "mergetle geosynch.tle molniya.tle visual.tle -se -oall.tle\n\n");
   printf( "would create a file 'all.tle',  containing elements from the three\n");
   printf( "input .tle files,  sorted by NORAD number.  If a satellite appears\n");
   printf( "in more than one .tle,  the .tle from the first file on the command\n");
   printf( "line is used.  Options are:\n\n");
   printf( "-sn, -sN       Sort output by ascending/descending NORAD number\n");
   printf( "-sc, -sC       Sort output by ascending/descending COSPAR desig\n");
   printf( "-sm, -sM       Sort output by ascending/descending mean motion\n");
   printf( "-se, -sE       Sort output by ascending/descending eccentricity\n");
   printf( "-sp, -sP       Sort output by ascending/descending epoch\n");
   printf( "-si, -sI       Sort output by ascending/descending inclination\n");
   printf( "-so, -sO       Sort output by ascending/descending ascending node\n");
   printf( "-sq            Sort output by ascending perigee\n");
   printf( "-o(filename)   Set name of output .tle file (default is out.tle)\n");
   printf( "-n             Remove names from input TLEs\n");
   printf( "-h             Remove HTML tags from input.  This allows you to extract\n");
   printf( "                 TLEs from certain Web pages.\n");
   exit( -1);
}

int main( const int argc, const char **argv)
{
   TLE *tles = (TLE *)calloc( MAX_TLES, sizeof( TLE));
   char *already_found = (char *)calloc( MAX_NORAD_NUMBER, sizeof( char));
   int n_found = 0, i, strip_names = 0;
   char sort_method = 0;
   const char *output_filename = "out.tle";
   FILE *ofile;

   if( argc < 2)
      error_exit( );
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'o':
               output_filename = argv[i] + 2;
               break;
            case 's':
               sort_method = argv[i][2];
               break;
            case 'n':
               strip_names = 1;
               printf( "Names will be removed in output\n");
               break;
            case 'h':
               heavens_above_html_tles = 1;
               printf( "HTML tags will be removed from input\n");
               break;
            default:
               printf( "Ignoring unknown option '%s'\n", argv[i]);
               break;
            }
      else
         {
         FILE *ifile = test_fopen( argv[i], "rb");
         int n;

         n_duplicates = 0;
         n = load_tles_from_file( ifile, tles + n_found, already_found);
         printf( "%d TLEs added from %s,  with %d duplicates found\n",
                            n, argv[i], n_duplicates);
         n_found += n;
         fclose( ifile);
         }
   free( already_found);
   if( strip_names)
      for( i = 0; i < n_found; i++)
         tles[i].name_line[0] = '\0';
#ifndef HAVE_REENTRANT_QSORT
   comparison_method = sort_method;
   qsort( tles, n_found, sizeof( TLE), tle_compare_for_qsort);
#else
   qsort_r( tles, n_found, sizeof( TLE), tle_compare_for_qsort_r, &sort_method);
#endif
   ofile = test_fopen( output_filename, "wb");
   for( i = 0; i < n_found; i++)
      show_tle( ofile, tles + i);
}
