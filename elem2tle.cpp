/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/* MOSTLY OBSOLETE.  See 'eph2tle.cpp' in the Find_Orb project
(https://github.com/Bill-Gray/find_orb) for a considerably better
approach to computing TLEs from orbital data.

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "comets.h"
#include "norad.h"
#include "norad_in.h"   /* for xke definition */
#include "date.h"

const double earth_mass_over_sun_mass = 2.98994e-6;
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define PI 3.141592653589793238462643383279502884197169399375105

int write_tle_from_vector( char *buff, const double *state_vect,
        const double epoch, const char *norad_desig, const char *intl_desig);

int verbose = 0;

static void set_tle_defaults( tle_t *tle)
{
   memset( tle, 0, sizeof( tle_t));
   strcpy( tle->intl_desig, "56999ZZ ");
   tle->classification = 'U';
   tle->ephemeris_type = '0';
}

#define centralize_angle(x) (fmod( (x) + PI * 10., PI + PI))

int vector_to_tle( tle_t *tle, const double *state_vect)
{
   ELEMENTS elem;
   int rval = 0, i;
   double tvect[6];
   const double max_ecc = .9999;

   for( i = 0; i < 6; i++)      /* cvt from km, km/min to AU, AU/day */
      tvect[i] = state_vect[i];
   elem.gm = xke * xke * earth_radius_in_km * earth_radius_in_km * earth_radius_in_km;
   calc_classical_elements( &elem, tvect, tle->epoch, 1);
   tle->xincl = centralize_angle( elem.incl);
   tle->xnodeo = centralize_angle( elem.asc_node);
   tle->omegao = centralize_angle( elem.arg_per);
   tle->xmo = centralize_angle( elem.mean_anomaly);

   if( elem.ecc > max_ecc || elem.major_axis <= 0.)
      rval = -1;
   else
      {
      tle->eo = elem.ecc;
      tle->xno = 1. / elem.t0;      /* xno is now in radians per minute */
      rval = 0;
      }
   if( tle->xincl < 0.)
      {
      tle->xincl = -tle->xincl;
      tle->xnodeo = centralize_angle( tle->xnodeo + PI);;
      tle->omegao = centralize_angle( tle->omegao + PI);;
      }
   return( rval);
}

static void show_results( const char *title, const tle_t *tle, const double *state_vect)
{
   if( title)
      printf( "%s\n", title);
   if( tle)
      {
      char buff[200];

      write_elements_in_tle_format( buff, tle);
      printf( "%s", buff);
      }
   printf("    %16.8f %16.8f %16.8f \n", state_vect[0], state_vect[1],
                                         state_vect[2]);
   printf("    %16.8f %16.8f %16.8f \n", state_vect[3] / 60.,
                                                state_vect[4] / 60.,
                                                state_vect[5] / 60.);
}

static int compute_new_state_vect( const tle_t *tle, double *state_vect,
                     const int ephem)
{
   double sat_params[N_SAT_PARAMS];
   int rval = 0;

   switch( ephem)
      {
      case 0:
         SGP_init( sat_params, tle);
         rval = SGP( 0., tle, sat_params, state_vect, state_vect + 3);
         break;
      case 1:
         SGP4_init( sat_params, tle);
         rval = SGP4( 0., tle, sat_params, state_vect, state_vect + 3);
         break;
      case 2:
         SGP8_init( sat_params, tle);
         rval = SGP8( 0., tle, sat_params, state_vect, state_vect + 3);
         break;
      case 3:
         SDP4_init( sat_params, tle);
         rval = SDP4( 0., tle, sat_params, state_vect, state_vect + 3);
         break;
      case 4:
         SDP8_init( sat_params, tle);
         rval = SDP8( 0., tle, sat_params, state_vect, state_vect + 3);
         break;
      default:
         printf( "??? ephem = %d\n", ephem);
         rval = -99;
         break;
      }
// if( rval)
//    printf( "??? rval = %d; ecc = %.6lf\n", rval, tle->eo);
   return( rval);
}

#define SIMPLEX_POINT struct simplex_point

SIMPLEX_POINT
   {
   double state_vect[6];
   double error;
   };

static double total_vector_diff( const double *vect1, const double *vect2)
{
   int i;
   double rval = 0.;

   for( i = 0; i < 6; i++)
      {
      double delta = vect1[i] - vect2[i];
      if( i >= 3)
         delta *= 1000.;
      rval += delta * delta;
      }
   return( rval);
}

static double compute_simplex_point_error( const double *state_vect, tle_t *tle,
            const double *start, const int ephem)
{
   double rval = 0., state_out[6];
   int compute_rval, vect_to_tle_rval;

   vect_to_tle_rval = vector_to_tle( tle, state_vect);
   if( vect_to_tle_rval == -1)
      return( 1.e+37);
   compute_rval = compute_new_state_vect( tle, state_out, ephem);
   if( compute_rval == SXPX_ERR_NEARLY_PARABOLIC
         || compute_rval == SXPX_ERR_NEGATIVE_MAJOR_AXIS
         || compute_rval == SXPX_ERR_NEGATIVE_XN
         || vect_to_tle_rval == -1)
      rval = 1.e+37;       /* invalid vector */
   else
      rval = total_vector_diff( state_out, start);
   return( rval);
}

static double try_simplex( SIMPLEX_POINT *simp, const double factor,
               tle_t *tle, const double *start, const int ephem)
{
   SIMPLEX_POINT new_point;
   int i, j;

   for( i = 0; i < 6; i++)
      {
      new_point.state_vect[i] = factor * simp->state_vect[i];
      for( j = 1; j < 7; j++)
         new_point.state_vect[i] += (1. - factor) * simp[j].state_vect[i] / 6.;
      }
   new_point.error = compute_simplex_point_error( new_point.state_vect, tle,
                                           start, ephem);
   if( new_point.error <= simp->error)
      *simp = new_point;
   return( new_point.error);
}

static void sort_simplexes( SIMPLEX_POINT *simp)
{
   int i;

   for( i = 0; i < 7; i++)          /* sort simplex points by error */
      if( simp[i].error < simp[i + 1].error)   /* highest to lowest */
         {
         SIMPLEX_POINT temp_elem = simp[i];

         simp[i] = simp[i + 1];
         simp[i + 1] = temp_elem;
         if( i)
            i -= 2;
         }
}

double dist_offset = 10000., vel_offset = 10.;

static void create_randomized_simplex( SIMPLEX_POINT *simp, const double *start_vect)
{
   int i;

   for( i = 0; i < 6; i++)
      {
      const double zval = (double)rand( ) / (double)RAND_MAX - .5;
      simp->state_vect[i] = start_vect[i]
                        + zval * (i < 3 ? dist_offset : vel_offset);
      }
}

static int initialize_simplexes( SIMPLEX_POINT *simp, const double *state_vect,
                                   const double *start_vect, const int ephem)
{
   int i, rval = 0;

   memcpy( simp[6].state_vect, start_vect, 6 * sizeof( double));
   assert( start_vect[0] && start_vect[1] && start_vect[2]);
   for( i = 0; i < 7 && !rval; i++)
      {
      tle_t tle;
      int iter = 0;

      set_tle_defaults( &tle);
      if( i != 6)
         create_randomized_simplex( simp + i, start_vect);
      while ( (simp[i].error = compute_simplex_point_error( simp[i].state_vect,
                          &tle, state_vect, ephem)) > 1e+36 && iter++ < 1000)
         create_randomized_simplex( simp, start_vect);
      if( iter >= 1000)
         rval = -1;
      }
   return( rval);
}

static int find_tle_via_simplex_method( tle_t *tle, const double *state_vect,
                     const double *start_vect, const int ephem)
{
   SIMPLEX_POINT simp[7];
   double best_rval_found = 1e+39, best_vect[6];
   int i, j, soln_found = 0, n_iterations = 0;
   int n_consecutive_contractions = 0;
   const int max_iterations = 43000;

   if( verbose)
      show_results( "Setting up:", NULL, start_vect);
   srand( 1);
   if( initialize_simplexes( simp, state_vect, start_vect, ephem))
      return( 0);       /* no solution found */
   while( !soln_found && n_iterations++ < max_iterations)
      {
      double ytry;

      sort_simplexes( simp);
      ytry = try_simplex( simp, -1., tle, state_vect, ephem);
      if( ytry <= simp[6].error)
         {
         if( verbose)
            {
            char buff[200];
            printf( "New record low: %f\n", ytry);

            write_elements_in_tle_format( buff, tle);
            printf( "%s", buff);
            }
         try_simplex( simp, 2., tle, state_vect, ephem);
         if( ytry < 1e-13)
            soln_found = true;
         if( ytry < best_rval_found)
            {
            best_rval_found = ytry;
            memcpy( best_vect, simp[0].state_vect, 6 * sizeof( double));
            }
         n_consecutive_contractions = 0;
         }
      else if( ytry > simp[1].error)
         {
         double ysave = simp[0].error;

         ytry = try_simplex( simp, .5, tle, state_vect, ephem);
         if( ytry > ysave)       /* still no success;  try contracting */
            {                    /* around lowest point: */
//          printf( "Contracting around best point\n");
            for( i = 0; i < 6; i++)
               {
               for( j = 0; j < 6; j++)
                  simp[i].state_vect[j] =
                           (simp[i].state_vect[j] + simp[6].state_vect[j]) / 2.;
               simp[i].error = compute_simplex_point_error( simp[i].state_vect, tle,
                                                               state_vect, ephem);
               }
            n_consecutive_contractions++;
            if( n_consecutive_contractions == 30)
               initialize_simplexes( simp, state_vect, best_vect, ephem);
            }
         else
            n_consecutive_contractions = 0;
         }
      if( n_iterations % 200 == 199)
         initialize_simplexes( simp, state_vect, best_vect, ephem);
      }
   sort_simplexes( simp);
   if( verbose)
      printf( "End err: %f\n", simp[6].error);
   vector_to_tle( tle, best_vect);
// vector_to_tle( tle, simp[6].state_vect);
   return( soln_found);
}

int compute_tle_from_state_vector( tle_t *tle, const double *state_vect, const int ephem,
                        double *trial_state)
{
   int n_failed_steps = 0, i;
   double state_out[6], best_vect[6], curr_err;
   const double thresh = 1e-12;

   memcpy( trial_state, state_vect, 6 * sizeof( double));
   if( vector_to_tle( tle, state_vect))
      {
      printf( "Immediate failure\n");
      return( -1);
      }
   memcpy( best_vect, state_vect, 6 * sizeof( double));
   compute_new_state_vect( tle, state_out, ephem);
   for( i = 0; i < 6; i++)
      trial_state[i] += state_vect[i] - state_out[i];
   curr_err = total_vector_diff( state_out, state_vect);
   if( verbose)
      show_results( "Initial guess", tle, state_out);
   if( curr_err < thresh)
      printf( "Got it right away\n");
   while( curr_err > thresh && n_failed_steps < 20)
      {
      double new_err = 0.;
      tle_t new_tle = *tle;

      if( vector_to_tle( &new_tle, trial_state))
         {
         memcpy( trial_state, best_vect, 6 * sizeof( double));
         show_results( "Simple failure:", tle, trial_state);
         return( -1);
         }
      compute_new_state_vect( &new_tle, state_out, ephem);
      new_err = total_vector_diff( state_out, state_vect);
      if( new_err > curr_err * .9)
          n_failed_steps++;      /* slow or no convergence */
      if( new_err < curr_err)
         {
         curr_err = new_err;
         *tle = new_tle;
         memcpy( best_vect, trial_state, 6 * sizeof( double));
         if( verbose)
            {
            printf( "New record %f\n", curr_err);
            show_results( NULL, tle, state_out);
            }
         }
      for( i = 0; i < 6; i++)
          trial_state[i] += state_vect[i] - state_out[i];
      }
   memcpy( trial_state, best_vect, 6 * sizeof( double));
   return( curr_err > thresh);
}

/* Main program */
int main( const int argc, const char **argv)
{
   const char *tle_filename = ((argc == 1) ? "test.tle" : argv[1]);
   FILE *ifile = fopen( tle_filename, "rb");
   tle_t tle; /* Pointer to two-line elements set for satellite */
   char line1[100], line2[100];
   int ephem = 1;       /* default to SGP4 */
   int i;               /* Index for loops etc */
   int n_failures = 0, n_simple = 0, n_simplex = 0;
   bool failures_only = false;

   for( i = 2; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'f':
               failures_only = true;
               break;
            case 'v':
               verbose = 1;
               break;
            case 'd':
               dist_offset = atof( argv[i] + 2);
               break;
            case 's':
               vel_offset = atof( argv[i] + 2);
               break;
            default:
               printf( "Option '%s' unrecognized\n", argv[i]);
               break;
            }
   if( !ifile)
      {
      printf( "Couldn't open input TLE file %s\n", tle_filename);
      exit( -1);
      }
   *line1 = '\0';
   while( fgets( line2, sizeof( line2), ifile))
      {
      int got_data = 0;
      double state_vect[6];

      set_tle_defaults( &tle);
      if( strlen( line2) > 110 && line2[7] == '.' && line2[18] == '.'
                     && line2[0] == '2' && line2[1] == '4')
         {
         got_data = 3;           /* Find_Orb state vector ephemeris */
         tle.epoch = atof( line2);
         sscanf( line2 + 13, "%lf %lf %lf %lf %lf %lf",
                    state_vect + 0, state_vect + 1, state_vect + 2,
                    state_vect + 3, state_vect + 4, state_vect + 5);
         }
      else if( strlen( line1) > 55 && !memcmp( line1 + 50, " (TDB)", 6))
         {                                  /* JPL Horizons vector */
         const double obliq_2000 = 23.4392911 * PI / 180.;

         tle.epoch = atof( line1);          /* get JD epoch from header... */
         strcpy( line1, line2);
         if( fgets( line2, sizeof( line2), ifile))
            got_data = 1;
         sscanf( line1, "%lf %lf %lf",
                    state_vect + 0, state_vect + 1, state_vect + 2);
         sscanf( line2, "%lf %lf %lf",
                    state_vect + 3, state_vect + 4, state_vect + 5);
                      /* Cvt ecliptic to equatorial 2000: */
         rotate_vector( state_vect    , obliq_2000, 0);
         rotate_vector( state_vect + 3, obliq_2000, 0);
         }
      else if( parse_elements( line1, line2, &tle) >= 0)
         got_data = 2;

      if( got_data == 1 || got_data == 3)
         tle.epoch -= 68.00 / 86400.;       /* rough convert TDT to UTC */

      if( got_data)     /* hey! we got a TLE! */
         {
         double sat_params[N_SAT_PARAMS],  trial_state[6];
         int simple_rval;
         bool failed = false;
         tle_t new_tle;

         if( got_data == 1 || got_data == 3)
            {
            ephem = 3;        /* Use SDP4 for JPL Horizons vectors */
            for( i = 0; i < 6 && fabs( state_vect[i]) < 1.; i++)
               ;
            if( i == 6)   /* all small quantities,  must be in AU & AU/day : */
               {
               for( i = 0; i < 6; i++)
                  state_vect[i] *= AU_IN_KM;
               for( i = 3; i < 6; i++)
                  state_vect[i] /= seconds_per_day;
               }
            for( i = 3; i < 6; i++)    /* cvt km/sec to km/min */
               state_vect[i] *= seconds_per_minute;
            if( !failures_only)
               show_results( "Before:", NULL, state_vect);
            }
         else
            {
            int is_deep = select_ephemeris( &tle);

            if( is_deep && (ephem == 1 || ephem == 2))
               ephem += 2;    /* switch to an SDx */
            if( !is_deep && (ephem == 3 || ephem == 4))
               ephem -= 2;    /* switch to an SGx */

            /* Calling of NORAD routines */
            /* Each NORAD routine (SGP, SGP4, SGP8, SDP4, SDP8)   */
            /* will be called in turn with the appropriate TLE set */
            switch( ephem)
               {
               case 0:
                  SGP_init( sat_params, &tle);
                  SGP( 0., &tle, sat_params, state_vect, state_vect + 3);
                  break;
               case 1:
                  SGP4_init( sat_params, &tle);
                  SGP4( 0., &tle, sat_params, state_vect, state_vect + 3);
                  break;
               case 2:
                  SGP8_init( sat_params, &tle);
                  SGP8( 0., &tle, sat_params, state_vect, state_vect + 3);
                  break;
               case 3:
                  SDP4_init( sat_params, &tle);
                  SDP4( 0., &tle, sat_params, state_vect, state_vect + 3);
                  break;
               case 4:
                  SDP8_init( sat_params, &tle);
                  SDP8( 0., &tle, sat_params, state_vect, state_vect + 3);
                  break;
               }
            if( !failures_only)
               show_results( "Before:", &tle, state_vect);
            }

         new_tle = tle;
         simple_rval = compute_tle_from_state_vector( &new_tle, state_vect, ephem, trial_state);
         if( simple_rval)
            {
            n_simplex++;
            find_tle_via_simplex_method( &new_tle, state_vect, trial_state, ephem);
            }
         else
            n_simple++;

         compute_new_state_vect( &new_tle, trial_state, ephem);
         for( i = 0; i < 6; i++)
            {
            trial_state[i] -= state_vect[i];
            if( fabs( trial_state[i]) > 1e-6)
               failed = true;
            }
         if( failed && failures_only)
            show_results( "Before:", &tle, state_vect);
         if( failed || !failures_only)
            show_results( (simple_rval ? "Simplex result:" : "Simplest method:"),
                                &new_tle, trial_state);
         if( failed)
            n_failures++;
         }
      strcpy( line1, line2);
      }
   fclose( ifile);
   printf( "%d solved with simple method; %d with simplex\n", n_simple, n_simplex);
   if( n_failures)
      printf( "%d failures\n", n_failures);
   return(0);
} /* End of main() */


