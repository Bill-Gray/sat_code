#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "watdefs.h"
#include "afuncs.h"
#include "comets.h"
#include "date.h"
#include "norad.h"
#include "mpc_func.h"
#include "observe.h"

/* Code to generate topocentric ephemerides from TLE data,  mostly focussed
on the TLEs provided in https://www.github.com/Bill-Gray/tles. The program
can be compiled for standalone use or for use with the on-line artsat
ephemeris service at https://www.projectpluto.com/sat_eph.htm (q.v.). */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

typedef struct
{
   double lat, lon, alt, rho_sin_phi, rho_cos_phi;
   double jd_start, jd_end, step_size;
   int n_steps;
   const char *desig;
} ephem_t;

static int verbose = 0;

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

static void show_base_60( char *buff, const unsigned n_millisec)
{
   snprintf( buff, 15, "%03u %02u %02u.%03u",
               n_millisec / 3600000u, (n_millisec / 60000u) % 60u,
               (n_millisec / 1000u) % 60u, n_millisec % 1000u);
}

static void put_ra_in_buff( char *buff, double ra)
{
   ra = fmod( ra, 2. * PI);
   if( ra < 0.)
      ra += PI + PI;
   show_base_60( buff, (unsigned)( 3600. * 1000. * ra * 12. / PI));
   memmove( buff, buff + 1, strlen( buff));     /* remove leading zero */
}

static void put_dec_in_buff( char *buff, const double dec)
{
   show_base_60( buff, (unsigned)( 3600. * 1000. * fabs( dec) * 180. / PI));
   *buff = (dec > 0. ? '+' : '-');
}

static double angle_between( const double *a, const double *b)
{
   const double cos_ang = dot_product( a, b) /
                              sqrt( dot_product( a, a) * dot_product( b, b));
   double rval = acose( cos_ang);

   return( rval * 180. / PI);
}

static inline bool desig_match( const tle_t *tle, const char *desig)
{
   size_t i = 0;
   bool rval = false;

   while( isdigit( desig[i]))
      i++;
   if( i == 5)
      {
      if( !desig[i])       /* desig is all digits -> it's the NORAD # */
         rval = (atoi( desig) == tle->norad_number);
      else
         {
         i = strlen( desig);
         if( i > 5 && i < 9)
            rval = !memcmp( tle->intl_desig, desig, i) && tle->intl_desig[i] <= ' ';
         }
      }
   return( rval);
}

/* Generates unit vector in the direction of ivect and stores it in z_vect;
a unit vector perpendicular to that in the xy plane,  stored in x_vect;
and a unit vector perpendicular to both,  stored in y_vect.   */

static double make_orthogonal_basis( const double *ivect,
              double *x_vect, double *y_vect, double *z_vect)
{
   double rval, len;

   memcpy( z_vect, ivect, 3 * sizeof( double));
   rval = normalize_vect3( z_vect);
   len = hypot( z_vect[0], z_vect[1]);
   x_vect[0] = z_vect[1] / len;
   x_vect[1] = -z_vect[0] / len;
   x_vect[2] = 0.;
   vector_cross_product( y_vect, z_vect, x_vect);
   return( rval);
}

/* 'obs_pos' = observer position relative to the geocenter, km
   'topo_posn' = artsat position relative to the observer, km
   'sat_vel' = artsat vel, km/minute (usual,  though somewhat oddball,  SxPx units)
   '*motion_pa' = returned posn angle of motion,  degrees
   apparent total angular motion is returned               */

static double compute_angular_rates( const double *obs_pos, const double *topo_posn,
              const double *sat_vel, double *motion_pa,
              double *ra_motion, double *dec_motion)
{
   double vel[3];         /* velocity of sat relative to observer */
   double x_vect[3];    /* unit vector in equatorial plane perpendicular to topo_posn */
   double y_vect[3];    /* unit vector perpendicular to topo_posn & x_vect */
   double z_vect[3];    /* unit vector version of topo_posn */
   double xmotion, ymotion, total_motion, dist;
   const double omega_E = 1.00273790934;
                   /* Earth rotations per sidereal day (non-constant) */
   const double omega = omega_E * 2. * PI / minutes_per_day;
                   /* Earth rotational rate in radians/minute */

   vel[0] = sat_vel[0] + omega * obs_pos[1];
   vel[1] = sat_vel[1] - omega * obs_pos[0];
   vel[2] = sat_vel[2];
   dist = make_orthogonal_basis( topo_posn, x_vect, y_vect, z_vect);
   xmotion = dot_product( vel, x_vect) / dist;      /* all in radians/minute */
   ymotion = dot_product( vel, y_vect) / dist;
   total_motion = hypot( xmotion, ymotion);
   *motion_pa = PI + atan2( xmotion, ymotion);
   *motion_pa *= 180. / PI;
   *ra_motion = xmotion * (180. / PI) * 60.;
   *dec_motion = ymotion * (180. / PI) * 60.;
   return( total_motion * (180. / PI) * 60.);      /* cvt to arcmin/min = arcsec/sec */
}

static char _header[200];
static int motion_units = 1;     /* default to '/minute = degrees/hr = "/sec */
static bool show_separate_motions = false;
static bool output_state_vectors = false;
static bool output_mjd = false;

static int show_ephems_from( const char *path_to_tles, const ephem_t *e,
                                  const char *filename, int start_line)
{
   FILE *ifile;
   char line0[100], line1[100], line2[100];
   int show_it = 1;
   double jd_tle = 0., tle_range = 1e+10, abs_mag = 0.;
   const bool is_geocentric = (e->rho_sin_phi == 0. && e->rho_cos_phi == 0.);
   static const char *header_text =
           "Date (UTC)  Time       R.A. (J2000)  decl   Azim   Alt  Elong"
           "  LuElo  Dist(km) \"/sec     PA";
   static const char *geo_header_text =
           "Date (UTC)  Time       R.A. (J2000)  decl   Elong  LuElo  Dist(km) \"/sec     PA";

   if( verbose)
      printf( "Should examine '%s'; start line %d\n", filename, start_line);
   snprintf( line0, sizeof( line0), "%s/%s", path_to_tles, filename);
   ifile = fopen( line0, "rb");
   if( !ifile)
      {
      fprintf( stderr, "'%s' not opened\n", line0);
      exit( 0);
      }
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), ifile))
      {
      tle_t tle;

      if( *line2 == '#')
         {
         char *tptr = strstr( line2, " H ");

         if( !memcmp( line2, "# MJD ", 6))
            {
            jd_tle = atof( line2 + 6) + 2400000.5;
            tle_range = 1.;
            show_it = (jd_tle < e->jd_end && jd_tle + tle_range > e->jd_start);
            }
         else if( tptr)
            {
            abs_mag = atof( tptr + 2);
            if( verbose)
               printf( "H = %.3f\n", abs_mag);
            }
         }
      else if( show_it && parse_elements( line1, line2, &tle) >= 0
                     && desig_match( &tle, e->desig))
         {
         double sat_params[N_SAT_PARAMS], jd = e->jd_start;
         size_t i, j;
         const int is_deep_type = select_ephemeris( &tle);

         if( is_deep_type)
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         if( verbose > 1)
            {
            printf( "Got TLEs for %f :\n", jd);
            printf( "%s\n%s\n%s\n", line0, line1, line2);
            }
         for( i = 0; i < (size_t)e->n_steps; i++,
                                         jd = e->jd_start + (double)i * e->step_size)
            if( (int)i >= start_line && jd >= jd_tle && jd < jd_tle + tle_range)
               {
               char buff[90], dec_buff[20], ra_buff[20], alt_buff[17];
               double pos[3], vel[3], obs_pos[3], ra, dec, dist;
               const double t_since = (jd - tle.epoch) * minutes_per_day;
               double solar_xyzr[4], lunar_xyzr[4], topo_posn[3], elong;
               double motion_rate, motion_pa;
               double ra_motion, dec_motion;
               const char *format_string;

               if( !i)
                  {
                  char *tptr;

                  printf( "\nEphemerides for %05d = %s%.2s-%s\n",
                              tle.norad_number,
                              (atoi( tle.intl_desig) > 57000) ? "19" : "20",
                              tle.intl_desig, tle.intl_desig + 2);
                  snprintf( _header, sizeof( _header),
                          "%s\n%s", line0, (is_geocentric ? geo_header_text : header_text));
                  if( show_separate_motions)
                     strcat( _header, "    RA \"/sec  dec");
                  if( motion_units == 60)
                     while( NULL != (tptr = strstr( _header, "/sec ")))
                        memcpy( tptr, "/min", 4);
                  strcat( _header, abs_mag ? "      Mag\n" : "\n");
                  if( output_state_vectors)
                     strcpy( _header, "Date (UTC)  Time"
                               "          x            y            z"
                               "          vx           vy           vz\n");
                  printf( "%s", _header);
                  }
               if( output_mjd)
                  snprintf( buff, sizeof( buff), "%.5f", jd - 2400000.5);
               else
                  full_ctime( buff, jd, FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
                                 | FULL_CTIME_LEADING_ZEROES);
               if( is_deep_type)
                  SDP4( t_since, &tle, sat_params, pos, vel);
               else
                  SGP4( t_since, &tle, sat_params, pos, vel);
               observer_cartesian_coords( jd, e->lon, e->rho_cos_phi,
                                        e->rho_sin_phi, obs_pos);
               get_satellite_ra_dec_delta( obs_pos, pos, &ra, &dec, &dist);
               epoch_of_date_to_j2000( jd, &ra, &dec);
               if( output_state_vectors)
                  {
                  const double year = 2000. + (jd - 2451545.) / 365.25;
                  double matrix[9];

                  setup_precession( matrix, year, 2000.);
                  precess_vector( matrix, pos, pos);
                  precess_vector( matrix, vel, vel);
                  printf( "%s %14.5f%14.5f%14.5f%11.5f%11.5f%11.5f\n",
                              buff, pos[0], pos[1], pos[2],
                              vel[0] / 60., vel[1] / 60., vel[2] / 60.);
                  }
               put_ra_in_buff( ra_buff, ra);
               put_dec_in_buff( dec_buff, dec);
               ra_buff[10] = dec_buff[9] = '\0';
               for( j = 0; j < 3; j++)
                  topo_posn[j] = pos[j] - obs_pos[j];
               motion_rate = compute_angular_rates( obs_pos, topo_posn, vel, &motion_pa,
                           &ra_motion, &dec_motion);
               lunar_solar_position( jd, lunar_xyzr, solar_xyzr);
               ecliptic_to_equatorial( solar_xyzr);
               ecliptic_to_equatorial( lunar_xyzr);
               if( !is_geocentric)
                  {
                  double x_vect[3], y_vect[3], z_vect[3], alt, az;

                  make_orthogonal_basis( obs_pos, x_vect, y_vect, z_vect);
                  az = PI + atan2( dot_product( x_vect, topo_posn),
                                   dot_product( y_vect, topo_posn));
                  az *= 180. / PI;
                  alt = 90. - angle_between( topo_posn, obs_pos);
                  snprintf( alt_buff, sizeof( alt_buff), " %5.1f %+05.1f",
                              az,  alt);
                  }
               else
                  *alt_buff = '\0';
               elong = angle_between( topo_posn, solar_xyzr);
               if( !output_state_vectors)
                  {
                  printf( "%s  %s  %s%s %6.1f %6.1f %8.0f", buff, ra_buff, dec_buff,
                     alt_buff, elong, angle_between( topo_posn, lunar_xyzr), dist);
                  motion_rate *= (double)motion_units;
                  if( motion_rate < 9.999)
                     format_string = "  %6.4f %6.1f";
                  else if( motion_rate < 99.99)
                     format_string = "  %6.3f %6.1f";
                  else if( motion_rate < 999.9)
                     format_string = "  %6.2f %6.1f";
                  else if( motion_rate < 9999.)
                     format_string = "  %6.1f %6.1f";
                  else
                     format_string = "  %6.0f %6.1f";
                  printf( format_string, motion_rate, motion_pa);
                  if( show_separate_motions)
                     {
                     const char precision = format_string[5];

                     snprintf( buff, sizeof( buff),
                                 "  %%+7.%cf %%+7.%cf", precision, precision);
                     printf( buff, ra_motion * (double)motion_units,
                                  dec_motion * (double)motion_units);
                     }

                  if( !abs_mag)
                     printf( "\n");
                  else
                     {
                     const double phase_ang = (180. - elong) * (PI / 180.);
                     double mag = abs_mag + 5. * log10( dist / AU_IN_KM)
                              + phase_angle_correction_to_magnitude(
                                       phase_ang, 0.15);
                     printf( "%8.1f\n", mag);
                     }
                  }
               start_line = (int)i + 1;
               }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( start_line);
}

static const char *tle_list_filename = "tle_list.txt";

int generate_artsat_ephems( const char *path_to_tles, const ephem_t *e)
{
   FILE *ifile;
   char buff[100];
   int is_in_range = 0, id_matches = 1, start_line = 0;

   snprintf( buff, sizeof( buff), "%s/%s", path_to_tles, tle_list_filename);
   if( verbose > 1)
      printf( "Opening '%s', looking for '%s'\n", buff, e->desig);
   ifile = fopen( buff, "rb");
   if( !ifile)
      {
      fprintf( stderr, "'%s' not opened\n", buff);
      exit( 0);
      }
   while( start_line != e->n_steps &&
                          fgets_trimmed( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "# Range:", 8))
         {
         char t_start[40], t_end[40];
         int n_scanned = sscanf( buff + 8, "%39s %39s", t_start, t_end);

         assert( 2 == n_scanned);
         if( get_time_from_string( 0., t_start, FULL_CTIME_YMD, NULL) < e->jd_end
                  && get_time_from_string( 0., t_end, FULL_CTIME_YMD, NULL) > e->jd_start)
            is_in_range = 1;
         }
      if( !memcmp( buff, "# ID:", 5))
         {
         if( buff[5] != ' ' || buff[11] != ' ' || buff[12] != ' ')
            fprintf( stderr, "BAD LINE %s\n", buff);
         for( int i = 6; i < 10; i++)
            if( !isdigit( buff[i]) || !isdigit( buff[i + 7]))
               {
               printf( "BAD LINE (2) %s\n", buff);
               i = 99;
               }
         if( strcmp( e->desig, buff + 13) && atoi( buff + 5) != atoi( e->desig))
            id_matches = 0;
         }
      if( !memcmp( buff, "# Include ", 10))
         {
         if( is_in_range && id_matches)
            start_line = show_ephems_from( path_to_tles, e, buff + 10, start_line);
         is_in_range = 0;
         id_matches = 1;
         }
      }
   fclose( ifile);
   if( start_line)
      printf( "%s", _header);
   return( start_line);
}

static int set_location( ephem_t *e, const char *mpc_code, const char *obscode_file_name)
{
   mpc_code_t c;
   int rval = get_lat_lon_info( &c, mpc_code);

   if( rval)
      {
      FILE *ifile = fopen( obscode_file_name, "rb");
      char buff[200];

      if( !ifile)
         {
         fprintf( stderr, "'%s' not found\n", obscode_file_name);
         exit( 0);
         }
      while( rval && fgets_trimmed( buff, sizeof( buff), ifile))
         if( !memcmp( mpc_code, buff, 3))
            {
            const int planet = get_mpc_code_info( &c, buff);

            if( planet != 3)
               {
               fprintf( stderr, "MPC code '%s' is for planet %d\n",
                           mpc_code, planet);
               exit( 0);
               }
            rval = 0;
            printf( "%s\n", c.name);
            }
      fclose( ifile);
      }
   if( !rval)
      {
      e->lat = c.lat;
      e->lon = c.lon;
      e->alt = c.alt;
      e->rho_cos_phi = c.rho_cos_phi;
      e->rho_sin_phi = c.rho_sin_phi;
      if( c.lon > PI)
         c.lon -= PI + PI;
      if( c.rho_sin_phi || c.rho_cos_phi)
         printf( "Latitude %c %f, Longitude %c %f\nAltitude %.1f meters (above WGS84 ellipsoid)\n",
                  (c.lat > 0. ? 'N' : 'S'), fabs( c.lat) * 180. / PI,
                  (c.lon > 0. ? 'E' : 'W'), fabs( c.lon) * 180. / PI,
                  c.alt);
      }
   return( rval);
}

#ifdef ON_LINE_VERSION
   #define OBSCODES_DOT_HTML_FILENAME  "/home/projectp/public_html/cgi-bin/fo/ObsCodes.htm"
   #define ROVERS_DOT_TXT_FILENAME     "/home/projectp/public_html/cgi-bin/fo/rovers.txt"
   #define PATH_TO_TLES                "/home/projectp/public_html/tles"
#else
   #define OBSCODES_DOT_HTML_FILENAME  "/home/phred/.find_orb/ObsCodes.htm"
   #define ROVERS_DOT_TXT_FILENAME     "/home/phred/.find_orb/rovers.txt"
   #define PATH_TO_TLES                "/home/phred/tles"
#endif

static const char *get_arg( const char **argv)
{
   const char *rval;

   if( argv[0][0] == '-' && argv[0][1])
      {
      if( !argv[0][2] && argv[1])
         rval = argv[1];
      else
         rval = argv[0] + 2;
      }
   else
      rval = NULL;
   if( !rval)
      {
      fprintf( stderr, "Can't get an argument : '%s'\n", argv[0]);
      exit( 0);
      }
   return( rval);
}

static void fix_desig( char *desig)
{
   size_t i;
   int bitmask = 0;

   for( i = 0; i < 10 && desig[i]; i++)
      if( isdigit( desig[i]))
         bitmask |= (1 << (int)i);
   if( i >= 9 && bitmask == 0xef && desig[4] == '-')
      {
      desig[0] = desig[2];                   /* it's in YYYY-NNNA form; */
      desig[1] = desig[3];                   /* cvt to YYNNNA form */
      for( i = 5; desig[i - 1]; i++)
         desig[i - 3] = desig[i];
      }
}

static void error_help( void)
{
   printf( "'sat_eph' arguments:\n"
           "   -c(MPC code) : specify location (default = geocentric)\n"
           "   -t(date/time) : starting time of ephemeris (default = now)\n"
           "   -n(#) : number of ephemeris steps (default = 20)\n"
           "   -s(#) : ephemeris step size in days (default = 1h)\n"
           "   -o(#) : five digit NORAD number or YYNNNA international designation\n"
           "   -r    : do _not_ round times to nearest step size\n"
           "   -u    : show motions in \"/min = degrees/hr (default is \"/sec)\n"
           "   -m    : show times as MJD\n"
           "   -V    : output state vectors instead of observables\n"
           "   -v(#) : level of verbosity\n");
}

int dummy_main( const int argc, const char **argv)
{
   int i;
   ephem_t e;
   bool round_to_nearest_step = true;
   const char *mpc_code = "500";
   const char *override_tle_filename = NULL;

   if( argc < 2)
      {
      error_help( );
      return( 0);
      }
   memset( &e, 0, sizeof( ephem_t));
   e.jd_start = current_jd( );
   e.n_steps = 20;
   e.step_size = 1. / 24.;
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1])
         {
         const char *arg = get_arg( argv + i);

         switch( argv[i][1])
            {
            case 'c':
               mpc_code = arg;
               break;
            case 'f':
               tle_list_filename = arg;
               break;
            case 'F':
               override_tle_filename = arg;
               break;
            case 't':
               e.jd_start = get_time_from_string( e.jd_start, arg, FULL_CTIME_YMD, NULL);
               break;
            case 'm':
               output_mjd = true;
               break;
            case 'n':
               e.n_steps = atoi( arg);
               break;
            case 'r':
               round_to_nearest_step = false;
               break;
            case 's':
               if( arg && *arg)
                  {
                  const char end_char = arg[strlen( arg) - 1];

                  e.step_size = atof( arg);
                  switch( end_char)
                     {
                     case 'h':
                        e.step_size /= hours_per_day;
                        break;
                     case 'm':
                        e.step_size /= minutes_per_day;
                        break;
                     case 's':
                        e.step_size /= seconds_per_day;
                        break;
                     }
                  }
               break;
            case 'S':
               show_separate_motions = true;
               break;
            case 'u':
               motion_units = 60;
               break;
            case 'o':
               /* Will handle below */
               break;
            case 'v':
               verbose = 1 + atoi( arg);
               break;
            case 'V':
               output_state_vectors = 1;
               break;
            default:
               fprintf( stderr, "Unrecognized option '%s'\n", argv[i]);
               error_help( );
               return( 0);
            }
         }
   if( set_location( &e, mpc_code, OBSCODES_DOT_HTML_FILENAME))
      if( set_location( &e, mpc_code, ROVERS_DOT_TXT_FILENAME))
         fprintf( stderr, "WARNING: Could not parse location '%s'\n", mpc_code);
   if( round_to_nearest_step && e.step_size)
      e.jd_start = floor( (e.jd_start - 0.5) / e.step_size) * e.step_size + 0.5;
   e.jd_end   = e.jd_start + (double)e.n_steps * e.step_size;
   if( verbose)
      printf( "arguments parsed;  JDs %f to %f\n", e.jd_start, e.jd_end);
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1] == 'o')
         {
         char desig[30];

         strncpy( desig, get_arg( argv + i), 29);
         fix_desig( desig);
         e.desig = desig;
         if( override_tle_filename)
            show_ephems_from( PATH_TO_TLES, &e, override_tle_filename, 0);
         else
            generate_artsat_ephems( PATH_TO_TLES, &e);
         }
   return( 0);
}

#ifndef ON_LINE_VERSION
int main( const int argc, const char **argv)
{
   return( dummy_main( argc, argv));
}
#else

#include <errno.h>
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

int main( const int unused_argc, const char **unused_argv)
{
   const char *argv[2000];
   const size_t max_buff_size = 40000;       /* room for 500 obs */
   char *buff = (char *)malloc( max_buff_size);
   char field[30], time_text[80];
   char num_steps[30], step_size[30], obs_code[40];
   FILE *lock_file = fopen( "lock.txt", "w");
   int cgi_status, i, argc = 9;
   bool round_step = false;
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
      printf( "<p> Your artsat ephemerides are very important to us! </p>");
      return( 0);
      }
   fprintf( lock_file, "We're in\n");
   *time_text = *num_steps = *step_size = *obs_code = '\0';
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
      if( !strcmp( field, "time") && strlen( buff) < sizeof( time_text))
         strcpy( time_text, buff);
      if( !strcmp( field, "obj_name"))
         {
         char *obj_name = (char *)malloc( strlen( buff) + 1);

         strcpy( obj_name, buff);
         argv[argc++] = "-o";
         argv[argc++] = obj_name;
         }
      else if( !memcmp( field, "obj_", 4))      /* selected an object check-box */
         {
         char *obj_name = (char *)malloc( strlen( field) - 3);

         strcpy( obj_name, field + 4);
         argv[argc++] = "-o";
         argv[argc++] = obj_name;
         }
      else if( !strcmp( field, "motion60"))
         motion_units = 60;
      if( !strcmp( field, "num_steps") && strlen( buff) < sizeof( num_steps))
         strcpy( num_steps, buff);
      if( !strcmp( field, "step_size") && strlen( buff) < sizeof( step_size))
         {
         char *tptr = strchr( buff, 'v');

         if( tptr)
            {
            verbose = atoi( tptr + 1);
            *tptr = '\0';
            }
         strcpy( step_size, buff);
         }
      if( !strcmp( field, "obs_code") && strlen( buff) < sizeof( obs_code))
         strcpy( obs_code, buff);
      if( !strcmp( field, "round_step"))
         round_step = true;
      if( !strcmp( field, "vectors"))
         output_state_vectors = true;
      if( !strcmp( field, "mjd"))
         output_mjd = true;
      if( !strcmp( field, "show_separate_motions"))
         argv[argc++] = "-S";
      }
   fprintf( lock_file, "Fields read\n");
   if( !round_step)
      argv[argc++] = "-r";
   argv[0] = "sat_eph";
   argv[1] = "-t";
   argv[2] = time_text;
   argv[3] = "-c";
   argv[4] = obs_code;
   argv[5] = "-n";
   argv[6] = num_steps;
   argv[7] = "-s";
   argv[8] = step_size;
   argv[argc] = NULL;
   for( i = 0; argv[i]; i++)
      fprintf( lock_file, "arg %d: '%s'\n", (int)i, argv[i]);
   dummy_main( argc, argv);
   fprintf( lock_file, "dummy_main called\n");
   free( buff);
   printf( "On-line Sat_eph compiled " __DATE__ " " __TIME__ " UTC-5h\n");
   printf( "See <a href='https://www.github.com/Bill-Gray/sat_code'>"
               "https://www.github.com/Bill-Gray/sat_code</a> for source code\n");
   printf( "</pre> </body> </html>");
   return( 0);
}
#endif
