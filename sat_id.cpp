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
#endif
#include "norad.h"
#include "observe.h"
#include "mpc_func.h"

#ifdef FRENCH_REPUBLICAN_CLOCK
   #define  hours_per_day      10
   #define minutes_per_hour   100
   #define seconds_per_minute 100
#else
   #define hours_per_day       24
   #define minutes_per_hour    60
   #define seconds_per_minute  60
#endif

#define seconds_per_hour   (seconds_per_minute * minutes_per_hour)
#define seconds_per_day    (seconds_per_hour * hours_per_day)
#define minutes_per_day    (minutes_per_hour * hours_per_day)

#define OBSERVATION struct observation

OBSERVATION
   {
   char text[81];
   double jd, ra, dec;
   double lon, rho_cos_phi, rho_sin_phi;
   };

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
      rval[i] = '\0';
      }
   return( rval);
}

static double dmy_to_jd( const int year, const int month, const double day)
{
   double rval = 1721059.5
         + (double)( year * 365 + year / 4 - year / 100 + year / 400) + day;
   int i;
   static const char month_len[12] = { 31, 28, 31, 30, 31, 30,
                                       31, 31, 30, 31, 30, 31 };

   assert( month >= 1 && month <= 12);
   assert( day >= 1. && day < 32.);
   assert( year > 1956 && year < 2050);
   for( i = 0; i < month - 1; i++)
      rval += (double)month_len[i];
   if( month < 3 && !(year % 4))    /* leap years,  January and February */
      if( !(year % 400) || (year % 100))
         rval--;
   return( rval);
}

static int get_mpc_data( OBSERVATION *obs, const char *buff)
{
   if( strlen( buff) != 80)
      return( -1);

   obs->jd = extract_date_from_mpc_report( buff, NULL);
   get_ra_dec_from_mpc_report( buff, NULL, &obs->ra, NULL,
                                     NULL, &obs->dec, NULL);
   strcpy( obs->text, buff);
   return( 0);
}

void make_config_dir_name( char *oname, const char *iname)
{
   strcpy( oname, getenv( "HOME"));
   strcat( oname, "/.find_orb/");
   strcat( oname, iname);
}

/* This loads up the file 'ObsCodes.html' into memory on its first call.
Then,  given an MPC code,  it finds the corresponding line and copies
it into 'station_code_data'.  It looks in several places for the file;
if you've installed Find_Orb,  it should be able to get it from the
~/.find_orb directory.  It also checks for the truncated 'ObsCodes.htm'
version of the file.      */

int verbose = 0;

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
         ifile = fopen( filenames[i], "rb");
      for( i = 0; !ifile && i < 2; i++)
         {
         char filename[255];

         make_config_dir_name( filename, filenames[i]);
         ifile = fopen( filename, "rb");
         }

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
      printf( "Station code '%s' not found.\n", mpc_code);
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
   return( cached_ptr ? 0 : -1);
}

/* Loads up MPC-formatted 80-column observations from a file.  Makes
a pass to find out how many observations there are,  allocates space
for them,  then reads them again to actually load the observations. */

static OBSERVATION *get_observations_from_file( FILE *ifile, size_t *n_found,
         const double t_low, const double t_high)
{
   int pass;
   OBSERVATION *rval = NULL, obs;

   memset( &obs, 0, sizeof( OBSERVATION));
   for( pass = 0; pass < 2; pass++)
      {
      char buff[100];
      size_t count = 0;

      fseek( ifile, 0L, SEEK_SET);
      while( fgets_trimmed( buff, sizeof( buff), ifile))
         if( !get_mpc_data( &obs, buff) && obs.jd > t_low
                                        && obs.jd < t_high)
            {
            if( rval)
               {
               char station_data[100];

               j2000_to_epoch_of_date( obs.jd, &obs.ra, &obs.dec);
               if( get_station_code_data( station_data, obs.text + 77))
                  printf( "FAILED to find MPC code %s\n", obs.text + 77);
               sscanf( station_data + 3, "%lf %lf %lf", &obs.lon,
                                        &obs.rho_cos_phi, &obs.rho_sin_phi);
               obs.lon *= PI / 180.;
               rval[count] = obs;
               }
            count++;
            }
         else if( !strcmp( buff, "COM ignore obs"))
            while( fgets_trimmed( buff, sizeof( buff), ifile))
               if( !strcmp( buff, "COM end ignore obs"))
                  break;
      if( !pass && count)
         rval = (OBSERVATION *)calloc( count, sizeof( OBSERVATION));
      *n_found = count;
      }
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

   if( !rval)        /* same IDs?  Then sort by JD of observation */
      rval = (aptr->jd > bptr->jd ? 1 : -1);
   return( rval);
}

/* Copied straight from 'mpc_obs.cpp' in Find_Orb.  See comments there. */

void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context)
{
#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux)
   qsort_r( base, n_elements, esize, compare, context);
#else
   size_t gap = 1;
   char *data = (char *)base;
   char *pivot = (char *)alloca( esize);

   while( gap < n_elements)
      gap = gap * 3 + 1;
   while( gap /= 3)
      {
      size_t j;
      const size_t spacing = esize * gap;

      for( j = gap; j < n_elements; j++)
         {
         char *tptr = data + j * esize;
         char *tptr2 = tptr - spacing;

         if( (compare)( tptr2, tptr, context) > 0)
            {
            memcpy( pivot, tptr, esize);
            memcpy( tptr, tptr2, esize);
            tptr = tptr2;
            tptr2 -= spacing;
            while( tptr2 >= base && (compare)( tptr2, pivot, context) > 0)
               {
               memcpy( tptr, tptr2, esize);
               tptr = tptr2;
               tptr2 -= spacing;
               }
            memcpy( tptr, pivot, esize);
            }
         }
      }
#endif
}

static void compute_offsets( double *dx, double *dy,
                double delta_ra, const double dec1, const double dec2)
{
   while( delta_ra > PI)
      delta_ra -= PI + PI;
   while( delta_ra < -PI)
      delta_ra += PI + PI;
   *dy = dec2 - dec1;
   *dx = delta_ra * cos( (dec2 + dec1) / 2.);
}

static double angular_sep( const double delta_ra, const double dec1,
            const double dec2)
{
   double dx, dy;

   compute_offsets( &dx, &dy, delta_ra, dec1, dec2);
   return( sqrt( dx * dx + dy * dy));
}

/* Out of all observations for a given object,  this function will pick
two that "best" describe the object's motion.  For that purpose,  we look
for a pair closest to 'optimal_dist' apart.  We also limit the separation
in time to 'max_time_sep';  that's to avoid a situation where the observations
are really close to the optimal distance apart,  but are actually from
different orbits.

This code thinks in terms of pairs of observations.  If somebody insists
on providing a single observation,  we duplicate it.

We also drop objects if they're moving slower than 'speed_cutoff',
set to help us ignore the slow guys that are almost certainly rocks. */

static size_t drop_extra_obs( OBSERVATION *obs, const size_t n_obs,
                  const double speed_cutoff)
{
   size_t i = 0, rval = 0;

   while( i < n_obs)
      {
      size_t j = 0;
      OBSERVATION *optr = obs + i;

      while( j < n_obs - i && !id_compare( optr, optr + j))
         j++;
      if( j == 1)    /* singleton observation */
         {
         obs[rval++] = *optr;
         obs[rval++] = *optr;
         }
      else        /* two or more obs:  pick two best */
         {
         size_t a, b, best_a = 1, best_b = 0;
         double speed = 0.;
         double best_score = 1e+30;

         for( a = 1; a < j; a++)
            for( b = 0; b < a; b++)
               {
               const double max_time_sep = 0.1;  /* .1 day = 2.4 hr */
               const double optimal_dist = PI / 180.;   /* one degree */
               const double dist = angular_sep( optr[b].ra - optr[a].ra,
                                 optr[b].dec, optr[a].dec);
               const double score = fabs( dist - optimal_dist);
               const double dt = optr[a].jd - optr[b].jd;

               assert( dt >= .0);
               if( best_score > score
                        && dt < max_time_sep && dt > 0.
                        && !memcmp( &optr[a].text[77], &optr[b].text[77], 3))
                  {
                  best_score = score;
                  best_a = a;
                  best_b = b;
                  speed = dist / dt;
                  }
               }
         speed *= 180. / PI;    /* cvt speed from radians/day to deg/day */
         speed /= minutes_per_day;    /* ...then to deg/hour = arcmin/minute */
         if( speed > speed_cutoff)    /* omit slow objects */
            {
            obs[rval++] = optr[best_b];
            obs[rval++] = optr[best_a];
            }
         }
      i += j;
      }
   return( rval);
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

static bool got_obs_in_range( const OBSERVATION *obs, size_t n_obs,
               const double jd_start, const double jd_end)
{
// if( obs[0].jd < jd_end && obs[n_obs - 1].jd > jd_start)
      while( n_obs--)
         {
         if( obs->jd > jd_start && obs->jd < jd_end)
            return( true);
         obs++;
         }
   return( false);
}

          /* The computed and observed motions should match,  but (obviously)
          only to some tolerance.  A tolerance of 0.001'/s seems to work. */
double motion_mismatch_limit = .001;

/* Given a set of MPC observations and a TLE file,  this function looks at
each TLE in the file and checks to see if that satellite came close to any
of the observations.  The function is called for each TLE file.
*/

int norad_id = 0;

double tle_start = 0., tle_range = 1e+9;

static int add_tle_to_obs( OBSERVATION *obs, const size_t n_obs,
             const char *tle_file_name, const double search_radius,
             const double max_revs_per_day)
{
   char line0[100], line1[100], line2[100];
   FILE *tle_file = fopen( tle_file_name, "rb");
   int rval = 0, n_tles_found = 0;
   bool check_updates = true;

   if( !tle_file)
      {
      printf( "Couldn't open TLE file %s\n", tle_file_name);
      return( -1);
      }
   if( verbose)
      printf( "Looking through TLE file '%s', %u obs, radius %f, max %f revs/day\n",
                 tle_file_name, (unsigned)n_obs, search_radius, max_revs_per_day);
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), tle_file))
      {
      tle_t tle;  /* Structure for two-line elements set for satellite */
      const double mins_per_day = 24. * 60.;

      if( got_obs_in_range( obs, n_obs, tle_start, tle_start + tle_range)
                 && parse_elements( line1, line2, &tle) >= 0
                 && (tle.ephemeris_type == 'H'
                 || tle.xno < 2. * PI * max_revs_per_day / mins_per_day)
                 && (!norad_id || norad_id == tle.norad_number))
         {                           /* hey! we got a TLE! */
         double sat_params[N_SAT_PARAMS];

         if( verbose > 1)
            printf( "TLE found:\n%s\n%s\n", line1, line2);
         n_tles_found++;
         if( select_ephemeris( &tle))
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         for( size_t idx = 0; idx < n_obs; idx += 2)
            if( is_in_range( obs[idx].jd, tle_start, tle_range))
               {
               OBSERVATION *optr = obs + idx;
               double dx, dy;
               double radius;
               double ra, dec, dist_to_satellite;
               int sxpx_rval;

               sxpx_rval = compute_artsat_ra_dec( &ra, &dec, &dist_to_satellite,
                              optr, &tle, sat_params);
               compute_offsets( &dx, &dy, ra - optr->ra, dec, optr->dec);
               radius = sqrt( dx * dx + dy * dy) * 180. / PI;
               if( !sxpx_rval && radius < search_radius)      /* good enough for us! */
                  {
                  double dx1, dy1;
                  const double dt = optr[1].jd - optr[0].jd;
                  double motion_diff;

                  compute_artsat_ra_dec( &ra, &dec, &dist_to_satellite,
                              optr + 1, &tle, sat_params);
                  compute_offsets( &dx1, &dy1, ra - optr[1].ra, dec, optr[1].dec);
                  dx1 -= dx;
                  dy1 -= dy;
                  motion_diff = sqrt( dx1 * dx1 + dy1 * dy1);
                  if( dt)           /* convert separations/dist into speeds */
                     {
                     dx /= dt;
                     dy /= dt;
                     motion_diff /= dt;
                     }
                  assert( dt >= 0.);
                  motion_diff *= 180. / PI;  /* now in degrees/day */
                  motion_diff /= minutes_per_day;   /* now in arcmin/second */
                  if( motion_diff < motion_mismatch_limit)
                     {
                     char obuff[200];
                     char full_intl_desig[20];
                     double xvel, yvel;
                     double motion_rate = 0., motion_pa = 0.;

                     compute_offsets( &xvel, &yvel, optr[1].ra - optr[0].ra,
                                                    optr[1].dec, optr[0].dec);
                     if( dt)
                        {
                        motion_pa = atan2( yvel, xvel) * 180. / PI + 90.;
                        if( motion_pa < 0.)
                           motion_pa += 180.;
                        motion_rate = sqrt( xvel * xvel + yvel * yvel);
                        motion_rate /= dt;
                        motion_rate *= 180. / PI;        /* now in degrees/day */
                        motion_rate /= 24.;              /* now in degrees/hr = arcsec/second */
                        }
                     line1[8] = line1[16] = '\0';
                     memcpy( line1 + 30, line1 + 11, 6);
                     line1[11] = '\0';
                     sprintf( full_intl_desig, "%s%.2s-%s",
                              (tle.intl_desig[0] < '5' ? "20" : "19"),
                              tle.intl_desig, tle.intl_desig + 2);
                     sprintf( obuff, "      %5dU = %-9s",
                           tle.norad_number, full_intl_desig);
                     sprintf( obuff + strlen( obuff),
                               "e=%.2f; P=%.1f min; i=%.1f",
                               tle.eo, 2. * PI / tle.xno,
                               tle.xincl * 180. / PI);
                     if( tle_checksum( line0))         /* object name given... */
                        sprintf( obuff + strlen( obuff), ": %s", line0);
                     obuff[79] = '\0';    /* avoid buffer overrun */
//                   sprintf( obuff + strlen( obuff), " motion %f", motion_diff);
                     strcat( obuff, "\n");
                     sprintf( obuff + strlen( obuff),
                        "             motion %5.2f\"/sec at PA %.1f; dist=%8.1f km; offset=%5.2f deg\n",
                            motion_rate, motion_pa,
                            dist_to_satellite, radius);
                              /* "Speed" is displayed in arcminutes/second,
                                  or in degrees/minute */
                     printf( "%s\n", obs[idx].text);
                     printf( "%s\n", obuff);
                     }
                  }
               }
         }
      else if( !memcmp( line2, "# No updates", 12))
         check_updates = false;
      else if( !memcmp( line2, "# Ephem range:", 14))
         {
         const double mjd_1970 = 40587.;     /* MJD for 1970 Jan 1 */
         double mjd_start, mjd_end;
         double curr_mjd = mjd_1970 + (double)time( NULL) / 86400.;

         sscanf( line2 + 14, "%lf %lf %lf\n", &mjd_start, &mjd_end, &tle_range);
         if( check_updates && mjd_end < curr_mjd + 7.)
            printf( "WARNING:  Update TLEs in '%s'\n", tle_file_name);
         if( !got_obs_in_range( obs, n_obs, mjd_start + 2400000.5,
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
         int year1, year2, month1, month2;
         double day1, day2;
         int n_read;

         n_read = sscanf( line2 + 9, "%d-%d-%lf %d-%d-%lf",
                     &year1, &month1, &day1,
                     &year2, &month2, &day2);
         assert( n_read == 6);
         tle_start = dmy_to_jd( year1, month1, day1);
         tle_range = dmy_to_jd( year2, month2, day2) - tle_start;
         }
      else if( !memcmp( line2, "# MJD ", 6))
         tle_start = atof( line2 + 6) + 2400000.5;
      else if( !memcmp( line2, "# Include ", 10))
         {
         if( got_obs_in_range( obs, n_obs, tle_start, tle_start + tle_range))
            {
            char iname[255];
            size_t i = strlen( tle_file_name);

            while( i && tle_file_name[i - 1] != '/')
               i--;
            memcpy( iname, tle_file_name, i);
            strcpy( iname + i, line2 + 10);
            rval = add_tle_to_obs( obs, n_obs, iname, search_radius,
                                    max_revs_per_day);
            }
         tle_start = 0.;
         tle_range = 1e+9;
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   if( verbose)
      printf( "%d TLEs read from '%s'\n", n_tles_found, tle_file_name);
   fclose( tle_file);
   return( rval);
}

/* The "on-line version",  sat_id2,  gathers data from a CGI multipart form,
   puts it into a file,  possibly adds in some options,  puts together the
 command-line arguments,  and then calls sat_id_main.  See 'sat_id2.cpp'. */

static double get_time( const char *buff)
{
   long year = atol( buff) / 10000L;
   const long month = (atol( buff) / 100L) % 100L;
   const int day = atol( buff) % 100L;

   if( year < 40)       /* before 2040 */
      year += 2000;
   if( year < 100)       /* 1940-1999 */
      year += 1900;
   return( dmy_to_jd( (int)year, (int)month, (double)day));
}

#ifdef ON_LINE_VERSION
int sat_id_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   const char *tle_file_name = "tle_list.txt";
   FILE *ifile = fopen( argv[1], "rb");
   OBSERVATION *obs;
   size_t n_obs;
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
   int rval;

   if( argc == 1)
      {
      printf( "No input file of astrometry specified on command line\n\n");
      error_exit( -2);
      }

   for( int i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         {
         const char *param = argv[i] + 2;

         if( !*param && i < argc - 1)
            param = argv[i + 1];
         switch( argv[i][1])
            {
            case 'a':
               t_low = get_time( param);
               break;
            case 'b':
               t_high = get_time( param);
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
            case 't':
               tle_file_name = param;
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
      for( int i = 0; i < argc; i++)
         printf( "Arg %d: '%s'\n", i, argv[i]);

   if( !ifile)
      {
      printf( "Couldn't open input file %s\n", argv[1]);
      return( -1);
      }
   obs = get_observations_from_file( ifile, &n_obs, t_low, t_high);
   fclose( ifile);
   printf( "%d observations found\n", (int)n_obs);
   shellsort_r( obs, n_obs, sizeof( obs[0]), compare_obs, NULL);
   n_obs = drop_extra_obs( obs, n_obs, speed_cutoff);
   printf( "%d observations left after dropping extras\n", (int)n_obs);
   if( !obs || !n_obs)
      return( -2);

   rval = add_tle_to_obs( obs, n_obs, tle_file_name, search_radius,
                                    max_revs_per_day);
   free( obs);
   get_station_code_data( NULL, NULL);
   printf( "%.1f seconds elapsed\n", (double)clock( ) / (double)CLOCKS_PER_SEC);
   return( rval);
}     /* End of main() */

