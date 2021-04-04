#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"
#include "norad.h"
#include "mpc_func.h"
#include "observe.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


typedef struct
{
   double lat, lon, alt, rho_sin_phi, rho_cos_phi;
   double jd_start, jd_end, step_size;
   int n_steps, norad_number;
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

static const char *header_text =
           "Date (UTC)  Time       R.A. (J2000)  decl   Alt   Elong   Dist(km)\n";

static int show_ephems_from( const char *path_to_tles, const ephem_t *e,
                                  const char *filename)
{
   FILE *ifile;
   char line0[100], line1[100], line2[100];
   int show_it = 1, n_lines_generated = 0;
   double jd_tle = 0., tle_range = 1e+10;

   if( verbose)
      printf( "Should examine '%s'\n", filename);
   snprintf( line0, sizeof( line0), "%s/%s", path_to_tles, filename);
   ifile = fopen( line0, "rb");
   assert( ifile);
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), ifile))
      {
      tle_t tle;

      if( !memcmp( line2, "# MJD ", 6))
         {
         jd_tle = atof( line2 + 6) + 2400000.5;
         tle_range = 1.;
         show_it = (jd_tle < e->jd_end && jd_tle + tle_range > e->jd_start);
         }
      else if( show_it && parse_elements( line1, line2, &tle) >= 0
                     && tle.norad_number == e->norad_number)
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
         for( i = 0; i < (size_t)e->n_steps; i++, jd += e->step_size)
            if( jd >= jd_tle && jd < jd_tle + tle_range)
               {
               char buff[90], dec_buff[20], ra_buff[20];
               double pos[3], vel[3], obs_pos[3], ra, dec, dist;
               const double t_since = (jd - tle.epoch) * minutes_per_day;
               double solar_xyzr[4], topo_posn[3];

               if( !i)
                  {
                  printf( "Ephemerides for %05d = %s%.2s-%s\n",
                              tle.norad_number,
                              (atoi( tle.intl_desig) > 57000) ? "19" : "20",
                              tle.intl_desig, tle.intl_desig + 2);
                  printf( "%s", header_text);
                  }
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
               put_ra_in_buff( ra_buff, ra);
               put_dec_in_buff( dec_buff, dec);
               ra_buff[10] = dec_buff[9] = '\0';
               for( j = 0; j < 3; j++)
                  topo_posn[j] = pos[j] - obs_pos[j];
               lunar_solar_position( jd, NULL, solar_xyzr);
               ecliptic_to_equatorial( solar_xyzr);
               printf( "%s  %s  %s %+05.1f %6.1f %8.0f\n", buff, ra_buff, dec_buff,
                     90. - angle_between( topo_posn, obs_pos),     /* alt */
                     angle_between( topo_posn, solar_xyzr), dist);
               n_lines_generated++;
               }
         }
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   fclose( ifile);
   return( n_lines_generated);
}

int generate_artsat_ephems( const char *path_to_tles, const ephem_t *e)
{
   FILE *ifile;
   char buff[200];
   int is_in_range = 0, id_matches = 1, ephem_lines_generated = 0;

   snprintf( buff, sizeof( buff), "%s/tle_list.txt", path_to_tles);
   ifile = fopen( buff, "rb");
   assert( ifile);
   while( ephem_lines_generated != e->n_steps &&
                          fgets_trimmed( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "# Range:", 8))
         {
         buff[19] = '\0';
         if( get_time_from_string( 0., buff + 9, FULL_CTIME_YMD, NULL) < e->jd_end
                  && get_time_from_string( 0., buff + 20, FULL_CTIME_YMD, NULL) > e->jd_start)
            is_in_range = 1;
         }
      if( !memcmp( buff, "# ID:", 5) && atoi( buff + 5) != e->norad_number)
         id_matches = 0;
      if( !memcmp( buff, "# Include ", 10))
         {
         if( is_in_range && id_matches)
            ephem_lines_generated += show_ephems_from( path_to_tles, e, buff + 10);
         is_in_range = 0;
         id_matches = 1;
         }
      }
   fclose( ifile);
   if( ephem_lines_generated)
      printf( "%s", header_text);
   return( ephem_lines_generated);
}

static int set_location( ephem_t *e, const char *mpc_code, const char *obscode_file_name)
{
   FILE *ifile = fopen( obscode_file_name, "rb");
   char buff[200];
   int got_it = 0, planet;

   assert( ifile);
   while( !got_it && fgets_trimmed( buff, sizeof( buff), ifile))
      if( !memcmp( mpc_code, buff, 3))
         {
         mpc_code_t c;

         planet = get_mpc_code_info( &c, buff);
         assert( planet == 3);
         e->lat = c.lat;
         e->lon = c.lon;
         e->alt = c.alt;
         e->rho_cos_phi = c.rho_cos_phi;
         e->rho_sin_phi = c.rho_sin_phi;
         got_it = 1;
         }
   fclose( ifile);
   return( got_it ? 0 : -1);
}

#ifdef ON_LINE_VERSION
   #define OBSCODES_DOT_HTML_FILENAME  "/home/projectp/public_html/cgi-bin/fo"
   #define PATH_TO_TLES                "/home/projectp/public_html/tles"
#else
   #define OBSCODES_DOT_HTML_FILENAME  "/home/phred/.find_orb/ObsCodes.htm"
   #define PATH_TO_TLES                "/home/phred/tles"
#endif

static void error_help( void)
{
   printf( "Mandatory 'sat_eph' arguments:\n"
           "   -c(MPC code) : specify location\n"
           "   -t(date/time) : starting time of ephemeris\n"
           "   -n(#) : number of ephemeris steps\n"
           "   -s(#) : ephemeris step size in days\n"
           "   -o(#) : five digit NORAD number\n"
           "Optional arguments:\n"
           "   -v(#) : level of verbosity\n");
}

int main( const int argc, const char **argv)
{
   int i;
   ephem_t e;

   if( argc < 3)
      {
      error_help( );
      return( -1);
      }
   memset( &e, 0, sizeof( ephem_t));
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1])
         {
         const char *arg = (argv[i][2] || i == argc - 1 ?
                                             argv[i] + 2 : argv[i + 1]);

         switch( argv[i][1])
            {
            case 'c':
               set_location( &e, arg, OBSCODES_DOT_HTML_FILENAME);
               break;
            case 't':
               e.jd_start = get_time_from_string( 0., arg, FULL_CTIME_YMD, NULL);
               break;
            case 'n':
               e.n_steps = atoi( arg);
               break;
            case 's':
               e.step_size = atof( arg);
               break;
            case 'o':
               e.norad_number = atoi( arg);
               break;
            case 'v':
               verbose = 1 + atoi( arg);
               break;
            default:
               error_help( );
               return( -1);
            }
         }
   e.jd_end   = e.jd_start + (double)e.n_steps * e.step_size;
   generate_artsat_ephems( PATH_TO_TLES, &e);
   return( 0);
}
