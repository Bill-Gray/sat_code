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
file named 'ALL_TLE.TXT'.  I've modified the default TLE to match. */

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

#define OBSERVATION struct observation
#define MAX_MATCHES 7

OBSERVATION
   {
   char text[81];
   double jd, ra, dec;
   double lon, rho_cos_phi, rho_sin_phi;
   };

typedef struct
{
   double dist;
   int norad_number;
   char intl_desig[9];
} match_t;

typedef struct
{
   OBSERVATION *obs;
   size_t idx1, idx2, n_obs;
   double speed;
   match_t matches[MAX_MATCHES];
} object_t;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define TIME_EPSILON (1./86400.)

static char *fgets_trimmed( char *buff, const int buffsize, FILE *ifile)
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

static int get_mpc_data( OBSERVATION *obs, const char *buff)
{
   obs->jd = extract_date_from_mpc_report( buff, NULL);
   if( !obs->jd)           /* not an 80-column MPC record */
      return( -1);
   get_ra_dec_from_mpc_report( buff, NULL, &obs->ra, NULL,
                                     NULL, &obs->dec, NULL);
   assert( strlen( buff) < sizeof( obs->text));
   strncpy( obs->text, buff, sizeof( obs->text));
   obs->text[80] = '\0';
   return( 0);
}

#if !defined( _WIN32)
void make_config_dir_name( char *oname, const char *iname)
{
   strcpy( oname, getenv( "HOME"));
   strcat( oname, "/.find_orb/");
   strcat( oname, iname);
}

static FILE *local_then_config_fopen( const char *filename, const char *permits)
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
static FILE *local_then_config_fopen( const char *filename, const char *permits)
{
   return( fopen( filename, permits));
}
#endif

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
      cached_data[size] = '\0';
      if( verbose)
         printf( "Station codes: %u bytes read\n", (unsigned)size);
      }
   if( !cached_ptr || memcmp( cached_ptr, mpc_code, 3))
      {
      char search_buff[5];

      sprintf( search_buff, "\n%.3s", mpc_code);
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

/* Loads up MPC-formatted 80-column observations from a file.  Makes
a pass to find out how many observations there are,  allocates space
for them,  then reads them again to actually load the observations. */

static OBSERVATION *get_observations_from_file( FILE *ifile, size_t *n_found,
         const double t_low, const double t_high)
{
   OBSERVATION *rval = NULL, obs;
   void *ades_context = init_ades2mpc( );
   char buff[400];
   size_t count = 0, n_allocated = 0;

   assert( ades_context);
   memset( &obs, 0, sizeof( OBSERVATION));
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      if( !get_mpc_data( &obs, buff) && obs.jd > t_low
                                     && obs.jd < t_high)
         {
         char station_data[100];

         j2000_to_epoch_of_date( obs.jd, &obs.ra, &obs.dec);
         if( !get_station_code_data( station_data, obs.text + 77))
            {
            sscanf( station_data + 3, "%lf %lf %lf", &obs.lon,
                                     &obs.rho_cos_phi, &obs.rho_sin_phi);
            obs.lon *= PI / 180.;
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
      else if( !strcmp( buff, "COM ignore obs"))
         while( fgets_trimmed( buff, sizeof( buff), ifile))
            if( !strcmp( buff, "COM end ignore obs"))
               break;
   *n_found = count;
   free_ades2mpc_context( ades_context);
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

static void error_exit( const int exit_code)
{
   printf(
"sat_id takes the name of an input file of MPC-formatted (80-column)\n\
astrometry as a command-line argument.  It searches for matches between\n\
the observation data and satellites in 'ALL_TLE.TXT'.  By default,  matches\n\
within .2 degrees are shown.\n\n\
Additional command-line arguments are:\n\
   -a YYYYMMDD  Only use observations after this time\n\
   -b YYYYMMDD  Only use observations before this time\n\
   -m (nrevs)   Only consider objects with fewer # revs/day (default=6)\n\
   -n (NORAD)   Only consider objects with this NORAD identifier\n\
   -r (radius)  Only show matches within this radius in degrees (default=4)\n\
   -t (fname)   Get TLEs from this filename\n\
   -v           Verbose output. '-v2' gets still more verboseness.\n\
   -z (rate)    Only consider observations above 'rate' deg/hr (default=.001)\n");
   exit( exit_code);
}

static int compute_artsat_ra_dec( double *ra, double *dec, double *dist,
         const OBSERVATION *optr, tle_t *tle,
         const double *sat_params)
{
   double observer_loc[3];
   double pos[3]; /* Satellite position vector */
   double t_since = (optr->jd - tle->epoch) * minutes_per_day;
   int sxpx_rval;

   observer_cartesian_coords( optr->jd, optr->lon,
           optr->rho_cos_phi, optr->rho_sin_phi, observer_loc);
   if( select_ephemeris( tle))
      sxpx_rval = SDP4( t_since, tle, sat_params, pos, NULL);
   else
      sxpx_rval = SGP4( t_since, tle, sat_params, pos, NULL);

   if( sxpx_rval == SXPX_WARN_PERIGEE_WITHIN_EARTH)
      sxpx_rval = 0;
   if( verbose > 2 && sxpx_rval)
      printf( "TLE failed for JD %f: %d\n", optr->jd, sxpx_rval);
   get_satellite_ra_dec_delta( observer_loc, pos, ra, dec, dist);
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

   /* By default,  we warn you if TLEs are going to run out next week. */

static double lookahead_warning_days = 7.;

          /* The computed and observed motions should match,  but
          (obviously) only to some tolerance.  A tolerance of 20
          arcseconds seems to work. (May seem large,  but these
          objects often 'streak' and have large residuals.)         */

double motion_mismatch_limit = 20.;

/* Given a set of MPC observations and a TLE file,  this function looks at
each TLE in the file and checks to see if that satellite came close to any
of the observations.  The function is called for each TLE file.
*/

int norad_id = 0;

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
   FILE *tle_file = fopen( tle_file_name, "rb");
   int rval = 0, n_tles_found = 0;
   bool check_updates = true;
   bool look_for_tles = true;
   const clock_t time_started = clock( );

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
         }
      if( is_a_tle && (tle.ephemeris_type == 'H'
                 || tle.xno < 2. * PI * max_revs_per_day / mins_per_day)
                 && (!norad_id || norad_id == tle.norad_number))
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
               int sxpx_rval, n_matches = 0;

               sxpx_rval = compute_artsat_ra_dec( &ra, &dec, &dist_to_satellite,
                              optr1, &tle, sat_params);
               radius = angular_sep( ra - optr1->ra, dec, optr1->dec, NULL) * 180. / PI;
               while( n_matches < MAX_MATCHES - 1
                       && obj_ptr->matches[n_matches].norad_number != tle.norad_number
                       && obj_ptr->matches[n_matches].norad_number)
                  n_matches++;
               if( !sxpx_rval && radius < search_radius      /* good enough for us! */
                       && radius < max_expected_error
                       && !obj_ptr->matches[n_matches].norad_number)
                  {
                  const double dt = optr2->jd - optr1->jd;
                  double motion_diff, ra2, dec2;
                  double temp_array[8];

                  assert( dt >= 0.);
                  compute_artsat_ra_dec( &ra2, &dec2, &dist_to_satellite,
                              optr2, &tle, sat_params);
                  temp_array[0] = ra;     /* starting point (computed) */
                  temp_array[1] = dec;
                  temp_array[2] = ra2;    /* ending point (computed) */
                  temp_array[3] = dec2;
                  temp_array[4] = optr1->ra;  /* starting point (observed) */
                  temp_array[5] = optr1->dec;
                  temp_array[6] = optr2->ra;  /* ending point (observed) */
                  temp_array[7] = optr2->dec;
                  motion_diff = relative_motion( temp_array);
                  motion_diff *= 3600. * 180. / PI;  /* cvt to arcseconds */
                  if( motion_diff < motion_mismatch_limit)
                     {
                     char obuff[200];
                     char full_intl_desig[20];
                     double motion_rate = 0., motion_pa = 0.;
                     size_t i;

                     motion_rate = angular_sep( optr1->ra - optr2->ra,
                                                 optr1->dec, optr2->dec, &motion_pa);
                     if( dt)
                        motion_rate /= dt;
                     line1[8] = line1[16] = '\0';
                     memcpy( line1 + 30, line1 + 11, 6);
                     line1[11] = '\0';
                     i = 0;
                     while( (int)i < n_matches && radius > obj_ptr->matches[i].dist)
                        i++;
                     memmove( obj_ptr->matches + i + 1, obj_ptr->matches + i,
                              (n_matches - i) * sizeof( match_t));
                     obj_ptr->matches[i].dist = radius;
                     obj_ptr->matches[i].norad_number = tle.norad_number;
                     strncpy( obj_ptr->matches[i].intl_desig,
                                                      tle.intl_desig, 9);
                     sprintf( full_intl_desig, "%s%.2s-%s",
                              (tle.intl_desig[0] < '5' ? "20" : "19"),
                              tle.intl_desig, tle.intl_desig + 2);
                     sprintf( obuff, "      %05dU = %-11s",
                           tle.norad_number, full_intl_desig);
                     if( tle.ephemeris_type != 'H')
                        sprintf( obuff + strlen( obuff),
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
                        sprintf( obuff + strlen( obuff), ": %s", line0);
                        }
                     obuff[79] = '\0';    /* avoid buffer overrun */
//                   sprintf( obuff + strlen( obuff), " motion %f", motion_diff);
                     strcat( obuff, "\n");
                     sprintf( obuff + strlen( obuff),
                        "             motion %7.4f\"/sec at PA %5.1f; dist=%8.1f km; offset=%6.3f deg\n",
                            motion_rate, motion_pa,
                            dist_to_satellite, radius);
                              /* "Speed" is displayed in arcminutes/second,
                                  or in degrees/minute */
                     printf( "%s\n", optr1->text);
                     printf( "%s", obuff);
                     if( dt)
                        {
                        motion_rate = angular_sep( ra - ra2, dec, dec2, &motion_pa);
                        if( dt)
                           motion_rate /= dt;
                        printf( "             motion %7.4f\"/sec at PA %5.1f (computed)\n",
                            motion_rate, motion_pa);
                        }
//                   printf( "             mismatch %5.2f\"\n", motion_diff);
                     printf( "\n");
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
      else if( !strncmp( line2, "# Ephem range:", 14))
         {
         const double mjd_1970 = 40587.;     /* MJD for 1970 Jan 1 */
         double mjd_start, mjd_end;
         double curr_mjd = mjd_1970 + (double)time( NULL) / 86400.;

         sscanf( line2 + 14, "%lf %lf %lf\n", &mjd_start, &mjd_end, &tle_range);
         if( check_updates && mjd_end < curr_mjd + lookahead_warning_days)
            {
            char time_buff[40];

            full_ctime( time_buff, mjd_end + 2400000.5, FULL_CTIME_YMD);
            printf( "WARNING: TLEs in '%s' run out on %s\n",
                           tle_file_name, time_buff);
            }
         if( !got_obs_in_range( objects, n_objects, mjd_start + 2400000.5,
                                            mjd_end + 2400000.5))
            {
            if( verbose)
               printf( "'%s' contains no TLEs for our time range\n", tle_file_name);
            fclose( tle_file);
            return( 0);
            }
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
            strcpy( iname + i, line2 + 10);
            if( verbose > 1)
               printf( "Including '%s'\n", iname);
            rval = add_tle_to_obs( objects, n_objects, iname, search_radius,
                                    max_revs_per_day);
            max_expected_error = saved_max_expected_error;
            }
         tle_start = 0.;
         tle_range = 1e+9;
         look_for_tles = true;
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   if( verbose)
      printf( "%d TLEs read from '%s', %.3f seconds\n", n_tles_found,
                tle_file_name,
                (double)( clock( ) - time_started) / (double)CLOCKS_PER_SEC);
   if( n_tles_found < n_tles_expected_in_file)
      {
      printf( "**** WARNING : %d TLEs were read from '%s'.  This is unexpected.\n",
                  n_tles_found, tle_file_name);
#ifdef ON_LINE_VERSION
      printf( "Please e-mail the author (pluto at projectpluto dot com) about this.\n");
#endif
      }
   fclose( tle_file);
   return( rval);
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
   double t_low = 2435839.5;  /* no satellites before 1957 Jan 1 */
   double t_high = 3000000.5;
   int rval, i, j, prev_i;
   bool show_summary = false;

   if( argc == 1)
      {
      printf( "No input file of astrometry specified on command line\n\n");
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
            case 'o':
               output_astrometry_filename = param;
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
//          case 'd':
//             debug_level = atoi( param);
//             break;
            case 'z':
               speed_cutoff = atof( param);
               break;
            default:
               printf( "Unrecognized command-line option '%s'\n", argv[i]);
               error_exit( -3);
               break;
            }
         }
   if( verbose)
      for( i = 0; i < argc; i++)
         printf( "Arg %d: '%s'\n", i, argv[i]);

   strcpy( tle_file_name, tname);
#if !defined( _WIN32)
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
         strcpy( buff, "tles/");
         strcat( buff, tname);
         make_config_dir_name( tle_file_name, buff);
         if( access( tle_file_name, F_OK))
            {
            if( verbose)
               fprintf( stderr, "'%s' failed\n", tle_file_name);
            strcpy( tle_file_name, getenv( "HOME"));
            strcat( tle_file_name, "/tles/");
            strcat( tle_file_name, tname);
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
      if( !i || id_compare( obs + i - 1, obs + i))
         n_objects++;
   objects = (object_t *)calloc( n_objects, sizeof( object_t));
   assert( objects);
   printf( "%d objects\n", (int)n_objects);
   for( n_objects = i = prev_i = 0; (size_t)i < n_obs; i++)
      if( !i || id_compare( obs + i - 1, obs + i))
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
      if( objects[i].speed >= speed_cutoff &&
              (objects[i].n_obs > 1 || include_singletons))
         {               /* fast enough to be considered */
         objects[j] = objects[i];
         j++;
         }
      }
   n_objects = j;
   printf( "%u objects after removing slow ones\n", (unsigned)n_objects);
   rval = add_tle_to_obs( objects, n_objects, tle_file_name, search_radius,
                                    max_revs_per_day);
   if( rval)
      fprintf( stderr, "Couldn't open TLE file %s\n", tname);
   else if( show_summary)
      for( i = 0; (size_t)i < n_objects; i++)
//       if( objects[i].matches[0])
            {
            char buff[30];

            printf( "\n%.12s ", objects[i].obs->text);
            for( size_t j = 0; objects[i].matches[j].norad_number; j++)
               printf( " %05d %s", objects[i].matches[j].norad_number,
                         unpack_intl( objects[i].matches[j].intl_desig, buff));
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

               for( i = 0; (size_t)i < n_objects; i++)
                  if( !memcmp( objects[i].obs->text, buff, 12) &&
                               objects[i].matches[0].norad_number)
                     {
                     fprintf( ofile, "COM %05dU = %s\n",
                         objects[i].matches[0].norad_number,
                         unpack_intl( objects[i].matches[0].intl_desig, tbuff));
                     objects[i].matches[0].norad_number = 0;
                     }
               }
            fputs( buff, ofile);
            }
         fclose( ofile);
         }
      }
   fclose( ifile);
   free( obs);
   free( objects);
   get_station_code_data( NULL, NULL);
   printf( "\n%.1f seconds elapsed\n", (double)clock( ) / (double)CLOCKS_PER_SEC);
   return( rval);
}     /* End of main() */
