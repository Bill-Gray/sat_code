/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
    sat_id.cpp     8 March 2003,  with updates as listed below

   An example 'main' function illustrating how to find which satellite(s)
are within a given radius of a given RA/dec,  as seen from a given
point.  The code reads in a file of observations in MPC format (name
provided as the first command-line argument).  For example:

sat_id mpc_obs.txt

   would hunt through the file 'mpc_obs.txt' for MPC-formatted
observations.  It would then read the file 'alldat.tle',  looking
for corresponding satellites within .2 degrees of said observations.
It then spits out the original file,  with satellite IDs added
(when found) after each observation line.  For each IDed satellite,
the international and NORAD designations are given,  along with
its angular distance from the search point,  position angle of
motion,  and apparent angular rate of motion in arcminutes/second
(or,  equivalently,  degrees/minute). */

/* 2 July 2003:  fixed the day/month/year to JD part of 'get_mpc_data()'
so it will work for all years >= 0 (previously,  it worked for years
2000 to 2099... plenty for the practical purpose of ID'ing recently-found
satellites,  but this is also an 'example' program.) */

/* 3 July 2005:  revised the check on the return value for parse_elements().
Now elements with bad checksums won't be rejected. */

/* 23 June 2006:  after comment from Eric Christensen,  revised to use
names 'ObsCodes.html' or 'ObsCodes.htm',  with 'stations.txt' being a
third choice.  Also added the '-a' command line switch to cause the program
to show all lines from input (default is now that only MPC astrometric
input gets echoed.)   */

/* 30 June 2006:  further comment from Eric Christensen:  when computing
object motion from two consecutive observations,  if the second one has
a date/time preceding the first,  you get a negative rate of motion that's
off by 180 degrees.  Fixed this. */

/* 17 Nov 2006:  artificial satellite data is now being provided in a
file named 'ALL_TLE.TXT'.  I've modified the default TLE to match.
(Note : since modified to be read from 'tle_list.txt',  as found in
the https://www.github.com/Bill-Gray/tles repository.  This allows
for multiple TLEs to be read.) */

/* 22 Oct 2012:  minor cosmetic changes,  such as making constant variables
of type 'const',  updating URL for the MPC station code file,  adding a
comment or two. */

/* 7 Jan 2013:  revised output to show satellite name if available,  plus
the eccentricity,  orbital period,  and inclination. */

/* 2013 Dec 8:  revised to pay attention to "# MJD" and "#Ephem start"
lines,  for files that contain many TLEs covering different time spans
for the same object.  I sometimes create such files;  when that happens,
for each observation,  only the TLE(s) covering that observation's time
should be used,  and the others are suppressed.       */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#if defined( _WIN32) || defined( __WATCOMC__)
   #include <malloc.h>     /* for alloca() prototype */
#else
   #include <unistd.h>
#endif

#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

#if defined( __has_include)
   #if !__has_include(<watdefs.h>)
      #error   \
        "You need the 'lunar' library (https://www.github.com/Bill-Gray/lunar).\
 See https://github.com/Bill-Gray/sat_code/issues/2 for a fix to this."
            /* This would normally be followed by dozens of errors.  GCC, */
            /* at least,  stops completely when it can't find a system */
            /* include file.  */
      #include <stop_cascading_errors>
   #endif
#endif

#include "norad.h"
#include "observe.h"
#include "watdefs.h"
#include "mpc_func.h"
#include "afuncs.h"
#include "date.h"
#include "sat_util.h"
#include "stringex.h"

#define OBSERVATION struct observation

OBSERVATION
   {
   char text[81];
   double jd, ra, dec;
   double observer_loc[3];
   };

typedef struct
{
   double dist, ra, dec, motion_rate, motion_pa;
   int norad_number;
   char intl_desig[9];
   char text[80];
   bool in_shadow;
} match_t;

typedef struct
{
   OBSERVATION *obs;
   size_t idx1, idx2, n_obs, n_matches;
   double speed;
   match_t *matches;
} object_t;

/* When we encounter a line for a spacecraft-based observation's offset
from the geocenter (a "second line" in the MPC's punched-card format system),
we store that offset in an offset_t structure.  Once the file has been
read in,  we go back and match the offsets to their corresponding observations.
If the lines aren't in proper order (which happens),  we'll fix it,  and
can detect errors where we get one line but not the other. */

typedef struct
{
   double jd, posn[3];
   char mpc_code[4];
} offset_t;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define TIME_EPSILON (1./86400.)
#define EARTH_MAJOR_AXIS    6378137.
#define EARTH_MINOR_AXIS    6356752.
#define is_power_of_two( X)   (!((X) & ((X) - 1)))

/* For testing purposes,  I sometimes want _only_ "my" TLEs (not the
   'classified' or Space-Watch TLEs) to be used.  This is invoked with -s. */
static bool my_tles_only = false;

#if !defined( ON_LINE_VERSION) && !defined( _WIN32)
   #define CSI "\x1b["
   #define OSC "\x1b]"
   #define REVERSE_VIDEO      CSI "0m" CSI "7m"
   #define NORMAL_VIDEO       CSI "0m"
#else
   #define REVERSE_VIDEO
   #define NORMAL_VIDEO
#endif

/* Sputnik 1 was launched on 1957 Oct 4.  There is no point in
checking for satellites before that date. TLEs,  and therefore
this program,  won't work after 2057 Jan 1. */

const double oct_4_1957 = 2436115.5;
const double jan_1_2057 = 2472364.5;

static int get_mpc_data( OBSERVATION *obs, const char *buff)
{
   obs->jd = extract_date_from_mpc_report( buff, NULL);
   if( obs->jd < oct_4_1957 || obs->jd > jan_1_2057)
      return( -1);            /* not an 80-column MPC record */
   if( get_ra_dec_from_mpc_report( buff, NULL, &obs->ra, NULL,
                                     NULL, &obs->dec, NULL))
      if( 's' != buff[14] && 'v' != buff[14])   /* satellite offsets and */
         return( -2);          /* roving obs lat/lons won't have RA/decs */
   assert( strlen( buff) < sizeof( obs->text));
   strncpy( obs->text, buff, sizeof( obs->text));
   obs->text[80] = '\0';
   return( 0);
}

static int extract_csv_value( char *obuff, const char *ibuff, int idx, size_t obuff_len)
{
   while( idx && *ibuff)
      if( *ibuff++ == ',')
         idx--;
   if( !*ibuff)
      return( -1);
   if( *ibuff == '"')
      ibuff++;
   obuff_len--;
   while( obuff_len-- && *ibuff && *ibuff != ',' && *ibuff != '"')
      *obuff++ = *ibuff++;
   *obuff = '\0';
   return( 0);
}

/* If we encounter 'field' data -- we're determining which artsats
are in a field,  not which ones match moving objects -- the output
is very different : */

static bool field_mode = false;

/* An imaging field is specified as an image name,  date/time, RA and
dec in decimal degrees,  and MPC obscode,  all comma-separated,  such as

MyImage,2023 nov 17 03:14:15.9,271.818,-14.142,T05
"Another image",2023-11-17.314159,21 16 58.21,"+31 41 59.26","V00"

In the first,  the RA and dec are in decimal degrees;  in the second,
sexagesimal notation is used and the RA is assumed to be in hours,
minutes,  and seconds.  You need the commas,  but can skip the
quotation marks.     */

static int get_field_data( OBSERVATION *obs, const char *buff)
{
   int rval = -1;
   size_t n_commas, i;

   for( i = n_commas = 0; buff[i] && n_commas < 5; i++)
      if( buff[i] == ',')
         n_commas++;
   if( n_commas == 4 && strlen( buff) < 80)
      {
      char tbuff[80];
      int bytes_read;

      extract_csv_value( tbuff, buff, 1, sizeof( tbuff));
      obs->jd = get_time_from_string( 0., tbuff, 0, NULL);
      if( obs->jd > oct_4_1957 && obs->jd < jan_1_2057)
         {
         extract_csv_value( tbuff, buff, 4, sizeof( tbuff));
         if( strlen( tbuff) != 3)     /* MPC code must be three characters */
            return( -1);
         strcpy( obs->text + 77, tbuff);
         extract_csv_value( tbuff, buff, 2, sizeof( tbuff));
         obs->ra = get_ra_from_string( tbuff, &bytes_read);
         if( bytes_read <= 0)
            return( -1);
         extract_csv_value( tbuff, buff, 3, sizeof( tbuff));
         obs->dec = get_dec_from_string( tbuff, &bytes_read);
         if( bytes_read <= 0)
            return( -1);
         extract_csv_value( obs->text, buff, 0, sizeof( obs->text));
         field_mode = true;
         rval = 0;
         }
      }
   return( rval);
}

/* This loads up the file 'ObsCodes.html' into memory on its first call.
Then,  given an MPC code,  it finds the corresponding line and copies
it into 'station_code_data'.  It looks in several places for the file;
if you've installed Find_Orb,  it should be able to get it from the
~/.find_orb directory.  It also checks for the truncated 'ObsCodes.htm'
version of the file.      */

int verbose = 0;
static bool check_all_tles = false;

static int get_station_code_data( char *station_code_data,
                  const char *mpc_code)
{
   static char *cached_data, *cached_ptr;

   if( !mpc_code)       /* freeing memory */
      {
      if( cached_data)
         free( cached_data);
      cached_data = cached_ptr = NULL;
      return( 0);
      }
   *station_code_data = '\0';
   if( !cached_data)
      {
      const char *filenames[2] = { "ObsCodes.html", "ObsCodes.htm" };
      FILE *ifile = NULL;
      size_t size;
      int i;

      for( i = 0; !ifile && i < 2; i++)
         ifile = local_then_config_fopen( filenames[i], "rb");
      if( !ifile)
         {
         printf( "Failed to find MPC station list 'ObsCodes.html'\n");
         printf( "This can be downloaded at:\n\n");
         printf( "http://www.minorplanetcenter.org/iau/lists/ObsCodes.html\n");
         exit( -3);
         }
      fseek( ifile, 0L, SEEK_END);
      size = (size_t)ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      cached_data = (char *)malloc( size + 1);
      if( fread( cached_data, 1, size, ifile) != size)
         {
         printf( "Failed to read station file\n");
         exit( -4);
         }
      fclose( ifile);
      if( verbose)
         printf( "Station codes: %u bytes read\n", (unsigned)size);
      ifile = local_then_config_fopen( "rovers.txt", "rb");
      if( ifile)
         {                             /* if 'rovers.txt' is available, */
         size_t size2;               /* append its data to ObsCodes.htm */

         fseek( ifile, 0L, SEEK_END);
         size2 = (size_t)ftell( ifile);
         fseek( ifile, 0L, SEEK_SET);
         cached_data = (char *)realloc( cached_data, size + size2 + 1);
         if( fread( cached_data + size, 1, size2, ifile) != size2)
            {
            printf( "Failed to read 'rovers.txt'\n");
            exit( -4);
            }
         fclose( ifile);
         size += size2;
         }
      cached_data[size] = '\0';
      }
   if( !cached_ptr || memcmp( cached_ptr, mpc_code, 3))
      {
      char search_buff[5];

      snprintf_err( search_buff, sizeof( search_buff), "\n%.3s", mpc_code);
      cached_ptr = strstr( cached_data, search_buff);
      if( cached_ptr)
         cached_ptr++;
      }
   if( cached_ptr)
      {
      size_t i;

      for( i = 0; cached_ptr[i] >= ' '; i++)
         station_code_data[i] = cached_ptr[i];
      station_code_data[i] = '\0';
      }
   else
      {
      static char *codes_already_reported = NULL;
      static size_t n_reported = 0;

      if( codes_already_reported && strstr( codes_already_reported, mpc_code))
         return( -1);      /* we've already reported this as 'unknown' */
      n_reported++;
      codes_already_reported = (char *)realloc( codes_already_reported,
                        n_reported * 4 + 1);
      if( n_reported == 1)    /* realloc( ) doesn't initialize to zero */
         *codes_already_reported = '\0';
      strcat( codes_already_reported, mpc_code);
      strcat( codes_already_reported, " ");
      printf( "Station code '%s' not found.\n", mpc_code);
      if( n_reported == 1)    /* only really need the following text once */
         {
#ifdef ON_LINE_VERSION
         printf( "If this is a new MPC code,  it could be that this service needs to be\n");
         printf( "updated to know about it.  Please contact pluto at projectpluto.com so\n");
         printf( "I can fix this.\n");
#else
         printf( "If this is a new MPC code,  you may need to get this file:\n");
         printf( "http://www.minorplanetcenter.org/iau/lists/ObsCodes.html\n");
         printf( "and replace the existing ObsCodes.html.\n");
#endif
         }
      }
   return( cached_ptr ? 0 : -1);
}

/* Code to check if a 'second line' (v for roving observer or s for
spacecraft observation) matches a 'first line' (the one that has the actual
astrometry).  If the date and obscode match,  and the observation doesn't
already have a non-geocentric position set,  then we have a match.  */

static bool offset_matches_obs( const offset_t *offset, const OBSERVATION *obs)
{
   if( offset->jd == obs->jd && !obs->observer_loc[0]
               && !strcmp( offset->mpc_code, obs->text + 77)
               && !obs->observer_loc[1] && !obs->observer_loc[2])
      return( true);
   else
      return( false);
}

/* (XXX) locations are specified with text such as

COM Long. 239 18 45 E, Lat. 33 54 11 N, Alt. 100m, Google Earth        */

static char xxx_location[80];

static int set_observer_location( OBSERVATION *obs)
{
   mpc_code_t code_data;
   int rval;

   if( memcmp( obs->text + 77, "XXX", 3))
      {
      char station_data[100];

      rval = get_station_code_data( station_data, obs->text + 77);
      if( !rval)
         get_mpc_code_info( &code_data, station_data);
      }
   else
      {
      static bool first_time = true;

      rval = get_xxx_location_info( &code_data, xxx_location);
      if( rval && first_time)
         printf( "WARNING: (XXX) observations won't be handled,  because the\n"
               "observatory location was either missing or incorrectly formatted.\n"
               "See https://www.projectpluto.com/xxx.htm for information on how\n"
               "this line should be handled.\n");
      first_time = false;
      }
   if( !rval)
      {
      observer_cartesian_coords( obs->jd, code_data.lon,
            code_data.rho_cos_phi, code_data.rho_sin_phi, obs->observer_loc);
      j2000_to_epoch_of_date( obs->jd, &obs->ra, &obs->dec);
      }
   return( rval);
}

/* Loads up MPC-formatted 80-column observations from a file.  Makes
a pass to find out how many observations there are,  allocates space
for them,  then reads them again to actually load the observations. */

static const char *_target_desig;

static OBSERVATION *get_observations_from_file( FILE *ifile, size_t *n_found,
         const double t_low, const double t_high)
{
   OBSERVATION *rval = NULL, obs;
   void *ades_context = init_ades2mpc( );
   char buff[400];
   size_t count = 0, n_allocated = 0, n_offsets = 0, i;
   offset_t *offsets = NULL;
   int n_errors_found = 0;

   assert( ades_context);
   memset( &obs, 0, sizeof( OBSERVATION));
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      if( (!get_mpc_data( &obs, buff) || !get_field_data( &obs, buff))
                   && obs.jd > t_low && obs.jd < t_high
                   && (!_target_desig || strstr( buff, _target_desig)))
         {
         if( buff[14] == 's' || buff[14] == 'v')
            {                 /* satellite obs or roving observer */
            offset_t toff;

            toff.jd = obs.jd;
            strlcpy_error( toff.mpc_code, buff + 77);
            if( buff[14] == 's')
               {
               int err_code = get_satellite_offset( buff, toff.posn);
               for( i = 0; i < 3; i++)
                  toff.posn[i] *= AU_IN_KM;
               ecliptic_to_equatorial( toff.posn);
               if( err_code < 0 && n_errors_found++ < 10)
                     fprintf( stderr, REVERSE_VIDEO
                                    "Malformed satellite offset; err %d\n%s\n"
                                    NORMAL_VIDEO, err_code, buff);
               }
            else
               {
               double lat, lon, alt_in_meters, rho_sin_phi, rho_cos_phi;

               if( sscanf( buff + 34, "%lf %lf %lf", &lon, &lat, &alt_in_meters) != 3)
                  if( n_errors_found++ < 10)
                     fprintf( stderr, "Couldn't parse roving observer:\n%s\n", buff);
               lon *= PI / 180.;
               lat *= PI / 180.;
               lat_alt_to_parallax( lat, alt_in_meters,
                          &rho_cos_phi, &rho_sin_phi,
                          EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS);
               observer_cartesian_coords( obs.jd, lon, rho_cos_phi,
                                        rho_sin_phi, toff.posn);
               }
            n_offsets++;
            offsets = (offset_t *)realloc( offsets, n_offsets * sizeof( offset_t));
            offsets[n_offsets - 1] = toff;
            }
         else if( !set_observer_location( &obs))
            {
            if( count == n_allocated)
               {
               n_allocated += 10 + n_allocated / 2;
               rval = (OBSERVATION *)realloc( rval,
                               (n_allocated + 1) * sizeof( OBSERVATION));
               }
            rval[count] = obs;
            count++;
            }
         }
      else if( !strncmp( buff, "COM verbose ", 12))
         verbose = atoi( buff + 12) + 1;
      else if( !strcmp( buff, "COM check tles"))
         check_all_tles = true;
      else if( !memcmp( buff, "COM Long.", 9))
         strlcpy_error( xxx_location, buff);
      else if( !strcmp( buff, "COM ignore obs"))
         while( fgets_trimmed( buff, sizeof( buff), ifile))
            if( !strcmp( buff, "COM end ignore obs"))
               break;
   *n_found = count;
   free_ades2mpc_context( ades_context);

            /* for each spacecraft offset,  look for the corresponding
            observation.  If we don't find one,  emit a warning. */
   for( i = 0; i < n_offsets; i++)
      {
      size_t j = 0;

      while( j < count && !offset_matches_obs( offsets + i, rval + j))
         j++;
      if( j == count)
         {
         if( n_errors_found++ < 10)
            fprintf( stderr, REVERSE_VIDEO "Unmatched artsat obs for JD %f\n"
                     NORMAL_VIDEO, offsets[i].jd);
         }
      else
         memcpy( rval[j].observer_loc, offsets[i].posn, 3 * sizeof( double));
      }
          /* Warn about obs with no offset data,  too. */
   for( i = 0; i < count; i++)
      if( !rval[i].observer_loc[0] && !rval[i].observer_loc[1] && !rval[i].observer_loc[2])
         if( n_errors_found++ < 10)
            fprintf( stderr, REVERSE_VIDEO "No position for this observation :\n%s\n"
                           NORMAL_VIDEO, rval[i].text);
   free( offsets);
   if( n_errors_found >= 10)
      fprintf( stderr, "Showing first ten of %d errors\n", n_errors_found);
   return( rval);
}

static int id_compare( const OBSERVATION *a, const OBSERVATION *b)
{
   return( memcmp( a->text, b->text, 12));
}

static int compare_obs( const void *a, const void *b, void *context)
{
   const OBSERVATION *aptr = (const OBSERVATION *)a;
   const OBSERVATION *bptr = (const OBSERVATION *)b;
   int rval = id_compare( aptr, bptr);

   INTENTIONALLY_UNUSED_PARAMETER( context);
   if( !rval)        /* same IDs?  Then sort by JD of observation */
      rval = (aptr->jd > bptr->jd ? 1 : -1);
   return( rval);
}

/* Copied straight from 'shellsor.cpp' in Find_Orb.  See comments there. */

void shellsort_r( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *, void *), void *context)
{
#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux)
   qsort_r( base, n_elements, elem_size, compare, context);
#else
   size_t gap = 250104703;
   char *data = (char *)base;
   char *pivot = (char *)alloca( elem_size);

   while( gap < n_elements)
      gap = gap * 8 / 3 + 1;
   while( (gap = gap * 3 / 8) != 0)
      {
      size_t j;
      const size_t spacing = elem_size * gap;

      for( j = gap; j < n_elements; j++)
         {
         char *tptr = data + j * elem_size;
         char *tptr2 = tptr - spacing;

         if( (compare)( tptr2, tptr, context) > 0)
            {
            memcpy( pivot, tptr, elem_size);
            memcpy( tptr, tptr2, elem_size);
            tptr = tptr2;
            tptr2 -= spacing;
            while( tptr2 >= base && (compare)( tptr2, pivot, context) > 0)
               {
               memcpy( tptr, tptr2, elem_size);
               tptr = tptr2;
               tptr2 -= spacing;
               }
            memcpy( tptr, pivot, elem_size);
            }
         }
      }
#endif
}

/* Given a unit vector,  this creates a perpendicular xi_vect
in the xy plane and an eta_vect perpendicular to them both. */

static void create_orthogonal_vects( const double *v, double *xi_vect, double *eta_vect)
{
   const double tval = sqrt( v[0] * v[0] + v[1] * v[1]);

   xi_vect[2] = 0.;
   if( !tval)     /* 'mid' is directly at a celestial pole */
      {
      xi_vect[0] = 1.;
      xi_vect[1] = 0.;
      }
   else
      {
      xi_vect[0] =  v[1] / tval;
      xi_vect[1] = -v[0] / tval;
      }
   vector_cross_product( eta_vect, v, xi_vect);
}

/* relative_motion() is intended for situations where you've
   computed that an object moved from p1 to p2,  and want to
   compare that motion vector to an observed motion from p3 to
   p4.  Initial use cases are in Sat_ID (we have a computed motion
   from TLEs and two observed RA/decs) and astcheck (similar,  but
   the computed positions are from orbital elements).  In each case,
   you're trying to determine if the observed and computed motions
   match to within some threshhold.

   I used to do this by computing (delta_RA * cos_dec, delta_dec)
   for both observed and computed motions.  That works well as long as
   the two points are very close together.  As they're separated,  you
   get distortions.  This should be more accurate.

   Given four points
         (ra_dec[0], ra_dec[1]) = starting point object 1
         (ra_dec[2], ra_dec[3]) = ending point object 1
         (ra_dec[4], ra_dec[5]) = starting point object 2
         (ra_dec[6], ra_dec[7]) = ending point object 2
   we compute their (xi, eta) sky plane coordinates,  using a plane
   tangent to the midpoint of the starting locations.  That way,  any
   distortion will affect both ends equally.  Then we compute how far
   each object moved in the sky plane coordinates.  Then we compute
   the differences in speed.              */

double relative_motion( const double *ra_dec)
{
   double mid[3], xi_vect[3], eta_vect[3];
   double xi[4], eta[4], v[4][3];
   double delta_xi, delta_eta;
   int i;

   for( i = 0; i < 4; i++)
      polar3_to_cartesian( v[i], ra_dec[i + i], ra_dec[i + i + 1]);

   for( i = 0; i < 3; i++)
      mid[i] = v[0][i] + v[1][i] + v[2][i] + v[3][i];
   normalize_vect3( mid);
   create_orthogonal_vects( mid, xi_vect, eta_vect);
   for( i = 0; i < 4; i++)
      {
      const double dist = dot_product( mid, v[i]);

      xi[i] = dot_product( xi_vect, v[i]) / dist;
      eta[i] = dot_product( eta_vect, v[i]) / dist;
      }
   delta_xi =  (xi[0]  - xi[1])  - (xi[2]  - xi[3]);
   delta_eta = (eta[0] - eta[1]) - (eta[2] - eta[3]);
   return( sqrt( delta_xi * delta_xi + delta_eta * delta_eta));
}

static double angular_sep( const double delta_ra, const double dec1,
            const double dec2, double *posn_ang)
{
   double p1[2], p2[2], dist;

   p1[0] = 0.;
   p1[1] = dec1;
   p2[0] = delta_ra;
   p2[1] = dec2;
   calc_dist_and_posn_ang( p1, p2, &dist, posn_ang);
   if( posn_ang)
      *posn_ang *= 180. / PI;
   return( dist);
}

/* Out of all observations for a given object,  this function will pick
two that "best" describe the object's motion.  For that purpose,  we look
for a pair closest to 'optimal_dist' apart.  We also limit the separation
in time to 'max_time_sep';  that's to avoid a situation where the observations
are really close to the optimal distance apart,  but are actually from
different orbits.

This code thinks in terms of pairs of observations.  If somebody insists
on providing a single observation,  we duplicate it.  (Unless the '-1'
switch is specified,  in which case single observations are just dropped.)
*/

static bool include_singletons = true;

static double find_good_pair( OBSERVATION *obs, const size_t n_obs,
                   size_t *idx1, size_t *idx2)
{
   size_t a, b;
   double speed = 0., dt;
   const double max_time_sep = 0.1;  /* .1 day = 2.4 hr */
   double best_score = 1e+30;

   *idx1 = *idx2 = 0;
   for( b = 0; b < n_obs; b++)
      for( a = b + 1; a < n_obs
                        && (dt = obs[a].jd - obs[b].jd) < max_time_sep; a++)
         if( !memcmp( &obs[a].text[77], &obs[b].text[77], 3))
            {
            const double optimal_dist = PI / 180.;   /* one degree */
            const double dist = angular_sep( obs[b].ra - obs[a].ra,
                              obs[b].dec, obs[a].dec, NULL);
            const double score = fabs( dist - optimal_dist);

            assert( dt >= .0);
            if( best_score > score)
               {
               best_score = score;
               *idx2 = a;
               *idx1 = b;
               speed = dist / dt;
               }
            }
   speed *= 180. / PI;    /* cvt speed from radians/day to deg/day */
   speed /= hours_per_day;    /* ...then to deg/hour = arcmin/minute */
   return( speed);
}

/* Aberration from the Ron-Vondrak method,  from Meeus'
_Astronomical Algorithms_, p 153,  just the leading terms */

static void compute_aberration( const double t_cen, double *ra, double *dec)
{
   const double l3 = 1.7534703 + 628.3075849 * t_cen;
   const double sin_l3 = sin( l3), cos_l3 = cos( l3);
   const double sin_2l3 = 2. * sin_l3 * cos_l3;
   const double cos_2l3 = 2. * cos_l3 * cos_l3 - 1.;
   const double x = -1719914. * sin_l3 - 25. * cos_l3
                       +6434. * sin_2l3 + 28007 * cos_2l3;
   const double y = 25. * sin_l3 + 1578089 * cos_l3
                +25697. * sin_2l3 - 5904. * cos_2l3;
   const double z = 10. * sin_l3 + 684185. * cos_l3
                +11141. * sin_2l3 - 2559. * cos_2l3;
   const double c = 17314463350.;    /* speed of light is 173.1446335 AU/day */
   const double sin_ra = sin( *ra), cos_ra = cos( *ra);

   *ra -= (y * cos_ra - x * sin_ra) / (c * cos( *dec));
   *dec += ((x * cos_ra + y * sin_ra) * sin( *dec) - z * cos( *dec)) / c;
}

static void error_exit( const int exit_code)
{
   printf(
"sat_id takes the name of an input file of MPC-formatted (80-column)\n\
astrometry as a command-line argument.  It searches for matches between\n\
the observation data and satellites in TLEs specified in 'tle_list.txt'.\n\
By default,  matches within .2 degrees are shown.\n\n\
Additional command-line arguments are:\n\
   -a YYYYMMDD  Only use observations after this time\n\
   -b YYYYMMDD  Only use observations before this time\n\
   -c           Check all TLEs for existence\n\
   -m (nrevs)   Only consider objects with fewer # revs/day (default=6)\n\
   -n (NORAD)   Only consider objects with this NORAD identifier\n\
   -r (radius)  Only show matches within this radius in degrees (default=4)\n\
   -t (fname)   Get TLEs from this filename\n\
   -v           Verbose output. '-v2' gets still more verboseness.\n\
   -y           Set tolerance for apparent motion mismatch\n\
   -z (rate)    Only consider observations above 'rate' deg/hr (default=.001)\n\
   \n\
   See 'sat_id.txt' for additional details.\n");
   exit( exit_code);
}

static int compute_artsat_ra_dec( double *ra, double *dec, double *dist,
         const OBSERVATION *optr, tle_t *tle,
         const double *sat_params, bool *in_shadow)
{
   double pos[3]; /* Satellite position vector */
   double sun_xyzr[4], tval;
   double t_since = (optr->jd - tle->epoch) * minutes_per_day;
   const double j2000 = 2451545.;      /* JD 2451545 = 2000 Jan 1.5 */
   int sxpx_rval;

   if( select_ephemeris( tle))
      sxpx_rval = SDP4( t_since, tle, sat_params, pos, NULL);
   else
      sxpx_rval = SGP4( t_since, tle, sat_params, pos, NULL);

   if( sxpx_rval == SXPX_WARN_PERIGEE_WITHIN_EARTH)
      sxpx_rval = 0;
   if( verbose > 2 && sxpx_rval)
      printf( "TLE failed for JD %f: %d\n", optr->jd, sxpx_rval);
   get_satellite_ra_dec_delta( optr->observer_loc, pos, ra, dec, dist);
   compute_aberration( (optr->jd - j2000) / 36525., ra, dec);
   lunar_solar_position( optr->jd, NULL, sun_xyzr);
   ecliptic_to_equatorial( sun_xyzr);
   tval = dot_product( sun_xyzr, pos);
   if( tval < 0. && in_shadow)    /* elongation greater than 90 degrees; */
      {                           /* may be in earth's shadow */
      double tvect[3], r2 = 0.;
      const double earth_r = EARTH_MAJOR_AXIS / 1000.;   /* in km */
      size_t i;

      tval /= sun_xyzr[3] * sun_xyzr[3];
      for( i = 0; i < 3; i++)
         {
         tvect[i] = pos[i] - sun_xyzr[i] * tval;
         r2 += tvect[i] * tvect[i];
         }
      *in_shadow = (r2 < earth_r * earth_r);
      }
   else if( in_shadow)
      *in_shadow = false;
   return( sxpx_rval);
}

static bool is_in_range( const double jd, const double tle_start,
                                             const double tle_range)
{
   return( !tle_range || !tle_start ||
            (jd >= tle_start && jd <= tle_start + tle_range));
}

/* Determines if we have _any_ observations between the given JDs.  If we
don't,  we can skip an individual TLE or an entire file.  */

static bool got_obs_in_range( const object_t *objs, size_t n_objects,
               const double jd_start, const double jd_end)
{
   while( n_objects--)
      {
      OBSERVATION *obs = objs->obs;
      size_t i;

      if( obs[0].jd < jd_end && obs[objs->n_obs - 1].jd > jd_start)
         for( i = 0; i < objs->n_obs; i++)
            if( obs[i].jd > jd_start && obs[i].jd < jd_end)
               return( true);
      objs++;
      }
   return( false);
}

static int _pack_intl_desig( char *desig_out, const char *desig)
{
   size_t i = 0;
   int rval = 0;

   while( desig[i] && isdigit( desig[i]))
      i++;
   memset( desig_out, ' ', 8);
   desig_out[8] = '\0';
   if( i == 5 && isupper( desig[5]))      /* already in YYYNNAaa form */
      {
      memcpy( desig_out, desig, 5);
      desig += 5;
      }
   else if( i == 4 && desig[4] == '-' && isdigit( desig[5])
                  && isdigit( desig[6]) && isdigit( desig[7]) && isupper( desig[8]))
      {
      desig_out[0] = desig[2];      /* desig is in YYYY-NNNAaa form */
      desig_out[1] = desig[3];
      desig_out[2] = desig[5];
      desig_out[3] = desig[6];
      desig_out[4] = desig[7];
      desig += 8;
      }
   else                 /* not a recognized intl desig form */
      rval = -1;
   if( !rval)
      {
      desig_out[5] = *desig++;
      if( isupper( *desig))
         desig_out[6] = *desig++;
      if( isupper( *desig))
         desig_out[7] = *desig++;
      }
   return( rval);
}

static int _compare_intl_desigs( const char *desig1, const char *desig2)
{
   char odesig1[9], odesig2[9];

   _pack_intl_desig( odesig1, desig1);
   _pack_intl_desig( odesig2, desig2);
   return( strcmp( odesig1, odesig2));
}

/* Code to look through 'sat_xref.txt',  if available,  and assign NORAD
and international designations to TLEs with only default designations.
See 'sat_xref.txt' and 'eph2tle.cpp' in the Find_Orb repository. */

static int look_up_extended_identifiers( const char *line0, tle_t *tle)
{
   static char buff[100];
   int match_found = !strcmp( line0, buff + 21);
   static bool got_sat_xref_txt = true;  /* until proven otherwise */

   if( !match_found && got_sat_xref_txt)
      {
      FILE *ifile = local_then_config_fopen( "sat_xref.txt", "rb");

      if( !ifile)    /* don't look for it again */
         got_sat_xref_txt = false;
      else
         {
         while( !match_found && fgets_trimmed( buff, sizeof( buff), ifile))
            match_found = !strcmp( line0, buff + 21);
         fclose( ifile);
         }
      }
   if( match_found)
      {
      tle->norad_number = atoi( buff);
      memcpy( tle->intl_desig, buff + 12, 8);
      }
   return( match_found);
}

/* The international (YYYY-NNNletter(s)) designation and the NORAD
five-digit designation are always shown for a match.  If they're in
the name as well,  we should remove them so that more of the actual
name gets displayed.  */

static void remove_redundant_desig( char *name, const char *desig)
{
   size_t len = strlen( desig), i;

   while( len && desig[len - 1] == ' ')
      len--;
   for( i = 0; name[i]; i++)
      if( !memcmp( name + i, desig, len))
         {
         size_t n = len;

         if( i >= 3 && name[i - 2] == '=' && name[i - 1] == ' '
                                          && name[i - 3] == ' ')
            {
            i -= 3;     /* remove preceding '=' */
            n += 3;
            }
         else if( name[i + n] == ' ' && name[i + n + 1] == '='
                                     && name[i + n + 2] == ' ')
            n += 3;
         memmove( name + i, name + i + n, strlen( name + i + n) + 1);
         i--;
         }
}

/* We check astrometry first against TLEs from github.com/Bill-Gray/tles,
then some other sources such as the amateur community's TLEs,  and
only then against Space-Track TLEs.  If we've already checked an
object against the previous sources,  then we really ought not to
check the Space-Track TLEs.  For one thing,  the mere fact that we've
gone to the effort of computing "our own" TLEs means the Space-Track
TLEs are unreliable (some have poor accuracy;  others are not consistently
available).

Therefore,  we maintain a list of '# ID:' numbers from tle_list.txt,  and when we
get to '# ID off',  we figure we're in Space-Track TLE territory.  If we
encounter a NORAD number that's in our list,  we essentially say :
we already found this object and have handled it. */

static bool already_found_desig( const int curr_norad, size_t n_norad_ids, const int *norad_ids)
{
   bool rval = false;

   while( !rval && n_norad_ids)
      rval = (curr_norad == norad_ids[--n_norad_ids]);
   return( rval);
}

   /* By default,  we warn you if TLEs are going to run out next week. */

static double lookahead_warning_days = 7.;

          /* The computed and observed motions should match,  but
          (obviously) only to some tolerance.  A tolerance of 60
          arcseconds seems to work. (May seem large,  but these
          objects often 'streak' and have large residuals.)         */

double motion_mismatch_limit = 60.;

/* Given a set of MPC observations and a TLE file,  this function looks at
each TLE in the file and checks to see if that satellite came close to any
of the observations.  The function is called for each TLE file.
*/

int norad_id = 0;
const char *intl_desig = NULL;

static int search_norad = 0;
static char search_intl[12];

double tle_start = 0., tle_range = 1e+9;

/* The following may be reset for a particular object where we have TLEs
that we expect have a certain maximum expected error,  using a line
starting with # Max error.  This lets us use a large search radius for
some objects with poor TLEs,  but restrict the radius tightly for objects
that we know we have a good handle on.  */

static double max_expected_error = 180.;
static int n_tles_expected_in_file = 0;

static int add_tle_to_obs( object_t *objects, const size_t n_objects,
             const char *tle_file_name, const double search_radius,
             const double max_revs_per_day)
{
   char line0[100], line1[100], line2[100];
   FILE *tle_file;
   int rval = 0, n_tles_found = 0;
   bool check_updates = true;
   bool look_for_tles = true;
   static bool error_check_date_ranges = true;
   const clock_t time_started = clock( );
   static int *norad_ids = NULL;
   static size_t n_norad_ids = 0;
   const double mjd_1970 = 40587.;     /* MJD for 1970 Jan 1 */
   const double curr_mjd = mjd_1970 + (double)time( NULL) / 86400.;

   if( !tle_file_name)     /* flag to free internal memory at shutdown */
      {
      if( norad_ids)
         free( norad_ids);
      n_norad_ids = 0;
      return( 0);
      }
   tle_file = fopen( tle_file_name, "rb");
   if( !tle_file)
      {
#ifdef ON_LINE_VERSION
      printf( "<h1> WARNING : '%s' not opened<br>\n", tle_file_name);
      printf( "Please e-mail pluto\x40projectpluto\x2e\x63om about this. </h1>\n");
#else
      fprintf( stderr, "WARNING : '%s' not opened\n", tle_file_name);
      fprintf( stderr, "Please e-mail pluto\x40projectpluto\x2e\x63om about this.\n");
#endif
      return( -1);
      }
   if( verbose)
      printf( "Looking through TLE file '%s', %u objs, radius %f, max %f revs/day\n",
                 tle_file_name, (unsigned)n_objects, search_radius, max_revs_per_day);
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), tle_file))
      {
      tle_t tle;  /* Structure for two-line elements set for satellite */
      const double mins_per_day = 24. * 60.;
      bool is_a_tle = false;

      if( verbose > 3)
         printf( "%s\n", line2);
      if( look_for_tles && parse_elements( line1, line2, &tle) >= 0)
         {
         is_a_tle = true;
         n_tles_found++;
         if( tle.norad_number == 99999)
            look_up_extended_identifiers( line0, &tle);
         if( (search_norad && search_norad != tle.norad_number)
               || (*search_intl && strcmp( search_intl, tle.intl_desig)))
            {
            fprintf( stderr, REVERSE_VIDEO "WARNING: desigs mismatch the 'ID' field : %s"
                        NORMAL_VIDEO "\n", tle_file_name);
            fprintf( stderr, "NORAD IDs are %d, %d\n", search_norad, tle.norad_number);
            fprintf( stderr, "Int'l IDs are '%s', '%s'\n", search_intl, tle.intl_desig);
            search_norad = 0;       /* suppress cascading errors */
            *search_intl = '\0';
            }
         }
      if( is_a_tle && (tle.ephemeris_type == 'H'
                 || tle.xno < 2. * PI * max_revs_per_day / mins_per_day)
                 && (!norad_id || norad_id == tle.norad_number)
                 && (!intl_desig || !_compare_intl_desigs( tle.intl_desig, intl_desig))
                 && (search_norad || !already_found_desig( tle.norad_number, n_norad_ids, norad_ids)))
         {                           /* hey! we got a TLE! */
         double sat_params[N_SAT_PARAMS];
         size_t idx;

         if( verbose > 1)
            printf( "TLE found:\n%s\n%s\n", line1, line2);
         if( select_ephemeris( &tle))
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         for( idx = 0; idx < n_objects; idx++)
            {
            object_t *obj_ptr = objects + idx;
            const OBSERVATION *optr1 = obj_ptr->obs + obj_ptr->idx1;
            const OBSERVATION *optr2 = obj_ptr->obs + obj_ptr->idx2;

            assert( obj_ptr->idx1 <= obj_ptr->idx2);
            assert( optr2->jd >= optr1->jd);
            if( is_in_range( optr1->jd, tle_start, tle_range))
               {
               double radius;
               double ra, dec, dist_to_satellite;
               int sxpx_rval;
               size_t i = 0;
               bool in_shadow;

               sxpx_rval = compute_artsat_ra_dec( &ra, &dec, &dist_to_satellite,
                              optr1, &tle, sat_params, &in_shadow);
               radius = angular_sep( ra - optr1->ra, dec, optr1->dec, NULL) * 180. / PI;
               while( i < obj_ptr->n_matches
                       && obj_ptr->matches[i].norad_number != tle.norad_number
                       && obj_ptr->matches[i].norad_number)
                  i++;
               if( !sxpx_rval && radius < search_radius      /* good enough for us! */
                       && radius < max_expected_error
                       && i == obj_ptr->n_matches)
                  {
                  double dt = optr2->jd - optr1->jd;
                  const double min_dt = 1e-6;   /* 0.0864 seconds */
                  double motion_diff, ra2, dec2;
                  double temp_array[8];
                  bool show_computed_motion = true;

                  assert( dt >= 0.);
                  if( !dt)
                     {
                     OBSERVATION temp_obs = *optr2;

                     temp_obs.jd += min_dt;
                     set_observer_location( &temp_obs);
                     if( vector3_length( optr2->observer_loc) > 6400.)
                        show_computed_motion = false;   /* spacecraft-based obs */
                     compute_artsat_ra_dec( &ra2, &dec2, &dist_to_satellite,
                              &temp_obs, &tle, sat_params, NULL);
                     }
                  else
                     compute_artsat_ra_dec( &ra2, &dec2, &dist_to_satellite,
                              optr2, &tle, sat_params, NULL);
                  temp_array[0] = ra;     /* starting point (computed) */
                  temp_array[1] = dec;
                  temp_array[2] = ra2;    /* ending point (computed) */
                  temp_array[3] = dec2;
                  temp_array[4] = optr1->ra;  /* starting point (observed) */
                  temp_array[5] = optr1->dec;
                  temp_array[6] = optr2->ra;  /* ending point (observed) */
                  temp_array[7] = optr2->dec;
                  if( !dt)
                     motion_diff = 0.;
                  else
                     motion_diff = relative_motion( temp_array);
                  motion_diff *= 3600. * 180. / PI;  /* cvt to arcseconds */
                  if( motion_diff < motion_mismatch_limit)
                     {
                     char obuff[200];
                     char full_intl_desig[20];
                     double motion_rate = 0., motion_pa = 0.;
                     const double arcminutes_per_radian = 60. * 180. / PI;

                     motion_rate = angular_sep( optr1->ra - optr2->ra,
                                          optr1->dec, optr2->dec, &motion_pa);
                     motion_rate *= arcminutes_per_radian;
                     if( dt)
                        motion_rate /= dt * minutes_per_day;
                     line1[8] = line1[16] = '\0';
                     memcpy( line1 + 30, line1 + 11, 6);
                     line1[11] = '\0';
                     i = 0;
                     while( i < obj_ptr->n_matches && radius > obj_ptr->matches[i].dist)
                        i++;
                     obj_ptr->matches = (match_t *)realloc( obj_ptr->matches,
                                    (obj_ptr->n_matches + 1) * sizeof( match_t));
                     memmove( obj_ptr->matches + i + 1, obj_ptr->matches + i,
                              (obj_ptr->n_matches - i) * sizeof( match_t));
                     memset( obj_ptr->matches + i, 0, sizeof( match_t));
                     obj_ptr->matches[i].dist = radius;
                     obj_ptr->matches[i].norad_number = tle.norad_number;
                     obj_ptr->n_matches++;
                     strncpy( obj_ptr->matches[i].intl_desig,
                                                      tle.intl_desig, 9);
                     snprintf_err( full_intl_desig, sizeof( full_intl_desig), "%s%.2s-%s",
                              (tle.intl_desig[0] < '5' ? "20" : "19"),
                              tle.intl_desig, tle.intl_desig + 2);
                     snprintf_err( obuff, sizeof( obuff), "     %05dU = %-11s ",
                           tle.norad_number, full_intl_desig);
                     if( tle.ephemeris_type != 'H')
                        snprintf_append( obuff, sizeof( obuff),
                               "e=%.2f; P=%.1f min; i=%.1f",
                               tle.eo, 2. * PI / tle.xno,
                               tle.xincl * 180. / PI);
                     if( tle_checksum( line0))         /* object name given... */
                        {
                        char norad_desig[20];

                        remove_redundant_desig( line0, full_intl_desig);
                        snprintf( norad_desig, sizeof( norad_desig),
                                           "NORAD %05d", tle.norad_number);
                        remove_redundant_desig( line0, norad_desig);
                        snprintf_append( obuff, sizeof( obuff), ": %s", line0);
                        }
                     strlcpy( obj_ptr->matches[i].text, obuff + 26, sizeof( obj_ptr->matches[i].text));
                     obuff[79] = '\0';    /* avoid buffer overrun */
//                   snprintf_append( obuff, sizeof( obuff), " motion %f", motion_diff);
                     strlcat_error( obuff, "\n");
                     if( !dt)
                        strlcat_error( obuff,
                              "             no observed motion (single obs) ");
                     else
                        snprintf_append( obuff, sizeof( obuff),
                              "             motion %7.4f\"/sec at PA %5.1f;",
                              motion_rate, motion_pa);
                     snprintf_append( obuff, sizeof( obuff),
                                   " dist=%8.1f km; offset=%7.4f deg\n",
                                   dist_to_satellite, radius);
                              /* "Speed" is displayed in arcminutes/second,
                                  or in degrees/minute */
                     if( verbose || !field_mode)
                        {
                        printf( "%s\n", optr1->text);
                        printf( "%s", obuff);
#ifdef SHOW_RA_DEC_OFFSETS
                        printf( "dRA = %.3f  dDec = %.3f\n",
                                    (ra - optr1->ra) * 180. / PI,
                                    (dec - optr1->dec) * 180. / PI);
#endif
                        }
                     motion_rate = angular_sep( ra - ra2, dec, dec2, &motion_pa);
                     motion_rate *= arcminutes_per_radian;
                     if( dt)
                        motion_rate /= dt * minutes_per_day;
                     else
                        motion_rate /= min_dt * minutes_per_day;
                     if( show_computed_motion && (verbose || !field_mode))
                        printf( "             motion %7.4f\"/sec at PA %5.1f (computed)\n\n",
                            motion_rate, motion_pa);
                     obj_ptr->matches[i].ra = ra;
                     obj_ptr->matches[i].dec = dec;
                     obj_ptr->matches[i].motion_rate = motion_rate;
                     obj_ptr->matches[i].motion_pa = motion_pa;
                     obj_ptr->matches[i].in_shadow = in_shadow;
                     }
                  }
               }
            }
         }
      else if( !strncmp( line2, "# No updates", 12))
         check_updates = false;
      else if( !strncmp( line2, "# Max error",  11))
         max_expected_error = atof( line2 + 12);
      else if( !strncmp( line2, "# TLEs expected:", 16))
         n_tles_expected_in_file = atoi( line2 + 17);
      else if( !strncmp( line2, "# Range_check ", 14))
         error_check_date_ranges = (line2[14] != '0');
      else if( !strncmp( line2, "# Ephem range:", 14))
         {
         double mjd_start, mjd_end, tle_step;
         int n_read;

         n_read = sscanf( line2 + 14, "%lf %lf %lf\n", &mjd_start, &mjd_end, &tle_step);
         assert( n_read == 3);
         if( error_check_date_ranges && tle_start)   /* do date ranges specified */
            {              /* in tle_list.txt and in the TLE file itself match? */
            const double tolerance = 0.001;     /* an allowance for roundoff */

            if( fabs( mjd_start + 2400000.5 - tle_start) > tolerance)
               {
               fprintf( stderr, REVERSE_VIDEO "WARNING: starting date for TLEs in '%s' "
                        "mismatches that in tle_list.txt\n" NORMAL_VIDEO, tle_file_name);
               fprintf( stderr, "TLE list start MJD %f\n", tle_start - 2400000.5);
               fprintf( stderr, "'Range:' line start MJD %f\n", mjd_start);
               fprintf( stderr, "diff = %f\n",
                        mjd_start + 2400000.5 - tle_start);
               }
            if( fabs( mjd_end + 2400000.5 - tle_start - tle_range) > tolerance)
               {
               fprintf( stderr, REVERSE_VIDEO "WARNING: ending date for TLES in '%s' "
                        "mismatches that in tle_list.txt\n" NORMAL_VIDEO, tle_file_name);
               fprintf( stderr, "TLE list ends MJD %f\n", tle_start + tle_range - 2400000.5);
               fprintf( stderr, "'Range:' line ends MJD %f\n", mjd_end);
               fprintf( stderr, "diff = %f\n",
                        mjd_end + 2400000.5 - tle_start - tle_range);
               }
            }
         tle_range = tle_step;
         if( check_updates && mjd_end < curr_mjd + lookahead_warning_days)
            {
            char time_buff[40];

            full_ctime( time_buff, mjd_end + 2400000.5, FULL_CTIME_YMD);
            fprintf( stderr, REVERSE_VIDEO "WARNING: TLEs in '%s' run out on %s (%.2f days)\n"
                           NORMAL_VIDEO, tle_file_name, time_buff, mjd_end - curr_mjd);
            }
         if( !got_obs_in_range( objects, n_objects, mjd_start + 2400000.5,
                                            mjd_end + 2400000.5) && !check_all_tles)
            {
            if( verbose)
               fprintf( stderr, REVERSE_VIDEO "'%s' contains no TLEs for our time range\n"
                               NORMAL_VIDEO, tle_file_name);
            fclose( tle_file);
            return( 0);
            }
         }
      else if( !memcmp( line2, "# ID: ", 6))
         {
         int n_read = sscanf( line2 + 6, "%d %s", &search_norad, search_intl);
         size_t i;

         assert( 2 == n_read);
         for( i = strlen( search_intl); i < 8; i++)
            search_intl[i] = ' ';
         search_intl[8] = '\0';
         }
      else if( !memcmp( line2, "# ID off", 8))
         {
         search_norad = 0;
         *search_intl = '\0';
         if( my_tles_only)
            break;
         }
      else if( !memcmp( line2, "# Range: ", 9))
         {
         int n_read;
         char start[20], end[20];

         n_read = sscanf( line2 + 9, "%19s %19s", start, end);
         assert( n_read == 2);
         tle_start = get_time_from_string( 0, start, FULL_CTIME_YMD, NULL);
         tle_range = get_time_from_string( 0, end, FULL_CTIME_YMD, NULL) - tle_start;
         if( !check_all_tles)
            look_for_tles = got_obs_in_range( objects, n_objects, tle_start,
                                 tle_start + tle_range);
         }
      else if( !memcmp( line2, "# MJD ", 6))
         {
         tle_start = atof( line2 + 6) + 2400000.5;
//       look_for_tles = got_obs_in_range( objects, n_objects, tle_start,
//                               tle_start + tle_range);
         }
      else if( !memcmp( line2, "# Include ", 10))
         {
         if( look_for_tles)
            {
            char iname[255];
            size_t i = strlen( tle_file_name);
            const double saved_max_expected_error = max_expected_error;

            while( i && tle_file_name[i - 1] != '/' && tle_file_name[i - 1] != '\\')
               i--;
            memcpy( iname, tle_file_name, i);
            strlcpy_err( iname + i, line2 + 10, sizeof( iname) - i);
            if( verbose > 1)
               printf( "Including '%s'\n", iname);
            i = 0;
            while( i < n_norad_ids && search_norad > norad_ids[i])
               i++;
            if( i == n_norad_ids || search_norad != norad_ids[i])
               {                           /* insert newly-found ID */
               if( !n_norad_ids)
                  norad_ids = (int *)malloc( sizeof( int));
               else if( is_power_of_two( n_norad_ids))
                  norad_ids = (int *)realloc( norad_ids, 2 * n_norad_ids * sizeof( int));
               memmove( norad_ids + i + 1, norad_ids + i, (n_norad_ids - i) * sizeof( int));
               norad_ids[i] = search_norad;
               n_norad_ids++;
               }
            rval = add_tle_to_obs( objects, n_objects, iname, search_radius,
                                    max_revs_per_day);
            max_expected_error = saved_max_expected_error;
            }
         tle_start = 0.;
         tle_range = 1e+9;
         look_for_tles = true;
         }
      else if( !memcmp( line2, "# Warn ", 7))
         {
         char *tptr = strchr( line2, '|');
         double jd_warn;

         assert( tptr);
         *tptr = '\0';
         jd_warn = get_time_from_string( 0, line2 + 7, FULL_CTIME_YMD, NULL);
         if( jd_warn < curr_mjd + 2400000.5 + lookahead_warning_days)
            fprintf( stderr, REVERSE_VIDEO "WARNING: %s" NORMAL_VIDEO "\n",
                           tptr + 1);
         }
      strlcpy_error( line0, line1);
      strlcpy_error( line1, line2);
      }
   if( verbose)
      printf( "%d TLEs read from '%s', %.3f seconds\n", n_tles_found,
                tle_file_name,
                (double)( clock( ) - time_started) / (double)CLOCKS_PER_SEC);
   if( n_tles_found < n_tles_expected_in_file)
      {
      const char *err_message = REVERSE_VIDEO
               "**** WARNING : %d TLEs were read from '%s'.  At least %d were expected.\n"
                                NORMAL_VIDEO;

#ifdef ON_LINE_VERSION
      printf( err_message, n_tles_found, tle_file_name, n_tles_expected_in_file);
#else
      fprintf( stderr, err_message, n_tles_found, tle_file_name, n_tles_expected_in_file);
#endif
      printf( "Please e-mail the author (pluto at projectpluto dot com) about this.\n");
      }
   fclose( tle_file);
   return( rval);
}

/* Output punch-card formatted astrometry is time-stamped when possible.
The scheme is the same as used for time-stamping NEOCP observations in
'neocp.cpp' and 'neocp2.cpp' in the 'miscell' repository (q.v) and makes
use of the fact that columns 57 to 65 default to being spaces. */

static void time_tag( char *tag)
{
   if( !memcmp( tag, "     ", 5))
      {
      const time_t t0 = time( NULL);
      struct tm tm;

#if defined( _WIN32) || defined( __WATCOMC__)
      memcpy( &tm, gmtime( &t0), sizeof( tm));
#else
      gmtime_r( &t0, &tm);
#endif
      tag[0] = '~';
      tag[1] = int_to_mutant_hex_char( tm.tm_mon + 1);
      tag[2] = int_to_mutant_hex_char( tm.tm_mday);
      tag[3] = int_to_mutant_hex_char( tm.tm_hour);
      tag[4] = int_to_mutant_hex_char( tm.tm_min);
      }
}

/* Given a 'packed' international designation such as 92044A or
10050BZ,  this outputs the 'unpacked/ desig 1992-044A or 2010-050BZ. */

static char *unpack_intl( const char *packed, char *unpacked)
{
   snprintf( unpacked, 12, "%s%.2s-%s",
            (atoi( packed) > 57000 ? "19" : "20"), packed, packed + 2);
   return( unpacked);
}

/* The "on-line version",  sat_id2,  gathers data from a CGI multipart form,
   puts it into a file,  possibly adds in some options,  puts together the
   command-line arguments,  and then calls sat_id_main.  See 'sat_id2.cpp'.

   By default,  TLE information is drawn from 'tle_list.txt';  this can
   be overridden with the -t (filename) option.  The specified file is first
   searched for in the current directory,  then in ~/.find_orb,  then in
   ~/.find_orb/tles,  then in ~/tles.    */

#ifdef ON_LINE_VERSION
int sat_id_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   char tle_file_name[256];
   const char *tname = "tle_list.txt";
   const char *output_astrometry_filename = NULL;
   bool output_only_matches = false;
   FILE *ifile;
   OBSERVATION *obs;
   object_t *objects;
   size_t n_obs, n_objects;
   double search_radius = 4;     /* default to 4-degree search */
            /* Asteroid searchers _sometimes_ report Molniyas to me,
            which make two revolutions a day.  This limit could safely
            be set to three,  but few artsats are between 6 and 3 revs/day
            (i.e.,  in four to eight-hour orbits).  So this doesn't result
            in much extra checking. */
   double max_revs_per_day = 6.;
            /* Anything slower than 0.001'/sec is almost certainly not an
            artsat.  We don't even bother to check those (unless the -z
            option is used to reset this limit).  */
   double speed_cutoff = 0.001;
   double t_low = oct_4_1957;
   double t_high = jan_1_2057;
   int rval, i, j, prev_i;
   bool show_summary = false, add_new_line = false, all_single = false;

   if( argc == 1)
      {
      fprintf( stderr, "No input file of astrometry specified on command line\n\n");
      error_exit( -2);
      }

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *param = argv[i] + 2;

         if( !*param && i < argc - 1)
            param = argv[i + 1];
         switch( argv[i][1])
            {
            case '1':
               include_singletons = false;
               break;
            case 'a':
               t_low = get_time_from_string( 0, param, FULL_CTIME_YMD, NULL);
               break;
            case 'b':
               t_high = get_time_from_string( 0, param, FULL_CTIME_YMD, NULL);
               break;
            case 'c':
               check_all_tles = true;
               break;
            case 'd':
               _target_desig = param;
               break;
            case 'i':
               intl_desig = param;
               break;
            case 'l':
               lookahead_warning_days = atof( param);
               break;
            case 'r':
               search_radius = atof( param);
               break;
            case 'y':
               motion_mismatch_limit = atof( param);
               break;
            case 'm':
               max_revs_per_day = atof( param);
               break;
            case 'n':
               norad_id = atoi( param);
               break;
            case 'N':
               add_new_line = true;
               break;
            case 'o':
            case 'O':
               output_astrometry_filename = param;
               output_only_matches = (argv[i][1] == 'o');
               break;
            case 's':
               my_tles_only = true;
               break;
            case 'S':      /* treat each observation individually */
               all_single = true;
               break;
            case 't':
               tname = param;
               break;
            case 'u':
               show_summary = true;
               break;
            case 'v':
               verbose = atoi( param) + 1;
               break;
            case 'z':
               speed_cutoff = atof( param);
               break;
            default:
               fprintf( stderr, "Unrecognized command-line option '%s'\n", argv[i]);
               error_exit( -3);
               break;
            }
         }
   if( verbose)
      for( i = 0; i < argc; i++)
         printf( "Arg %d: '%s'\n", i, argv[i]);

   strlcpy_error( tle_file_name, tname);
#if !defined( _WIN32) && !defined( __WATCOMC__)
   if( access( tle_file_name, F_OK))
      {
      if( verbose)
         fprintf( stderr, "'%s' failed\n", tle_file_name);
      make_config_dir_name( tle_file_name, tname);
      if( access( tle_file_name, F_OK))
         {
         char buff[256];

         if( verbose)
            fprintf( stderr, "'%s' failed\n", tle_file_name);
         strlcpy_error( buff, "tles/");
         strlcat_error( buff, tname);
         make_config_dir_name( tle_file_name, buff);
         if( access( tle_file_name, F_OK))
            {
            if( verbose)
               fprintf( stderr, "'%s' failed\n", tle_file_name);
            strlcpy_error( tle_file_name, getenv( "HOME"));
            strlcat_error( tle_file_name, "/tles/");
            strlcat_error( tle_file_name, tname);
            }
         }
      }
#endif

   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open input file %s\n", argv[1]);
      perror( NULL);
      return( -1);
      }
   obs = get_observations_from_file( ifile, &n_obs, t_low, t_high);
   printf( "%d observations found\n", (int)n_obs);
   if( !obs || !n_obs)
      return( -2);
   shellsort_r( obs, n_obs, sizeof( obs[0]), compare_obs, NULL);

   for( n_objects = i = 0; (size_t)i < n_obs; i++)
      if( !i || id_compare( obs + i - 1, obs + i) || field_mode || all_single)
         n_objects++;
   objects = (object_t *)calloc( n_objects, sizeof( object_t));
   assert( objects);
   if( !field_mode)
      printf( "%d objects\n", (int)n_objects);
   for( n_objects = i = prev_i = 0; (size_t)i < n_obs; i++)
      if( !i || id_compare( obs + i - 1, obs + i) || field_mode || all_single)
         {
         objects[n_objects].obs = obs + i;
         if( n_objects)
            objects[n_objects - 1].n_obs = i - prev_i;
         n_objects++;
         prev_i = i;
         }
   objects[n_objects - 1].n_obs = i - prev_i;
   for( i = j = 0; (size_t)i < n_objects; i++)
      {
      objects[i].speed = find_good_pair( objects[i].obs,
                  objects[i].n_obs, &objects[i].idx1, &objects[i].idx2);
      if( (objects[i].speed >= speed_cutoff && objects[i].n_obs > 1)
            || (include_singletons && objects[i].n_obs == 1))
         {               /* fast enough to be considered */
         objects[j] = objects[i];
         j++;
         }
      }
   n_objects = j;
   if( !field_mode)
      printf( "%u objects after removing slow ones\n", (unsigned)n_objects);
   else
      max_revs_per_day = 20.;   /* for field-finding,  list everything */
   rval = add_tle_to_obs( objects, n_objects, tle_file_name, search_radius,
                                    max_revs_per_day);
   if( rval)
      fprintf( stderr, "Couldn't open TLE file %s\n", tname);
   else if( show_summary)
      {
      int n_matched = 0;
      size_t n_unmatched = 0;

      for( i = 0; (size_t)i < n_objects; i++)
         if( objects[i].n_matches)        /* first show those that were IDed */
            {
            char buff[30];

            printf( "\n%.12s ", objects[i].obs->text);
            for( j = 0; j < (int)objects[i].n_matches; j++)
               printf( " %05d %s", objects[i].matches[j].norad_number,
                            unpack_intl( objects[i].matches[j].intl_desig, buff));
            if( objects[i].n_matches == 1)
               {
               char *tptr = strchr( objects[i].matches[0].text, ':');

               assert( tptr);              /* if only one match,  output */
               printf( " %s", tptr + 1);   /* the object name */
               }
            if( objects[i].n_matches)
               n_matched++;
            }
      for( i = 0; (size_t)i < n_objects; i++)
         if( !objects[i].n_matches)          /* now show those _not_ IDed */
            {
            if( !n_unmatched)
               printf( "\nUnmatched objects :");
            printf( "%s%.12s", (n_unmatched % 5 == 0 ? "\n" : ""), objects[i].obs->text);
            n_unmatched++;
            }
      if( n_objects)
         printf( "\n%d matched of %d tracklets (%.1f%%)\n", n_matched, (int)n_objects,
                  (double)n_matched * 100. / (double)n_objects);
      }
   else if( field_mode)
      {
      const char *header =
           "Field        RA  (J2000) dec  '/min   PA   NORAD Int'l desig Name";

      printf( "%s", header);
      for( i = 0; (size_t)i < n_objects; i++)
         for( j = 0; j < (int)objects[i].n_matches; j++)
            {
            double ra = objects[i].matches[j].ra;
            double dec = objects[i].matches[j].dec;
            char buff[30];

            epoch_of_date_to_j2000( objects[i].obs->jd, &ra, &dec);
            printf( "\n%-12.12s %7.3f %7.3f %7.2f %5.1f", objects[i].obs->text,
                           ra * 180. / PI, dec * 180. / PI,
                           objects[i].matches[j].motion_rate,
                           objects[i].matches[j].motion_pa);
            printf( " %c %05d %-11s",
                         objects[i].matches[j].in_shadow ? '*' : ' ',
                         objects[i].matches[j].norad_number,
                         unpack_intl( objects[i].matches[j].intl_desig, buff));
            printf( "%s", strchr( objects[i].matches[j].text, ':') + 1);
            }
      printf( "\n%s", header);
      }
   if( output_astrometry_filename)
      {
      FILE *ofile = fopen( output_astrometry_filename, "wb");

      if( !ofile)
         {
         fprintf( stderr, "Couldn't open '%s' :", output_astrometry_filename);
         perror( NULL);
         }
      else
         {
         char buff[256];

         fseek( ifile, 0L, SEEK_SET);
         while( fgets( buff, sizeof( buff), ifile))
            {
            if( strlen( buff) > 80 && buff[80] < ' ')
               {
               char tbuff[30];
               bool was_matched = false;

               time_tag( buff + 59);
               for( i = 0; (size_t)i < n_objects; i++)
                  if( !memcmp( objects[i].obs->text, buff, 12))
                     {
                     if( objects[i].n_matches > 0
                             && objects[i].matches[0].norad_number > 0)
                        {
                        const char *common_name =
                                 strchr( objects[i].matches[0].text, ':') + 1;

                        was_matched = true;
                        unpack_intl( objects[i].matches[0].intl_desig, tbuff);
                        if( add_new_line)
                           if( !memcmp( tbuff + 5, "999", 3) || !memcmp( tbuff + 5, "000", 3))
                              {
                              fprintf( ofile, "\nCOM =%s   %05dU = %s\n",
                                           common_name, objects[i].matches[0].norad_number,
                                           tbuff);
                              *tbuff = '\0';
                              }
                        if( *tbuff)
                           fprintf( ofile, "%sCOM %05dU = %s   %s\n",
                               (add_new_line ? "\n" : ""),
                               objects[i].matches[0].norad_number,
                               tbuff, common_name);
                        objects[i].matches[0].norad_number = -1;
                        }
                     }
               if( was_matched && output_only_matches)
                  fputs( buff, ofile);
               }
            if( !output_only_matches)
               fputs( buff, ofile);
            }
         fclose( ofile);
         }
      }
   fclose( ifile);
   free( obs);
   for( i = 0; (size_t)i < n_objects; i++)
      if( objects[i].matches)
         free( objects[i].matches);
   free( objects);
   get_station_code_data( NULL, NULL);
   add_tle_to_obs( NULL, 0, NULL, 0., 0.);
   printf( "\n%.1f seconds elapsed\n", (double)clock( ) / (double)CLOCKS_PER_SEC);
   return( rval);
}     /* End of main() */
