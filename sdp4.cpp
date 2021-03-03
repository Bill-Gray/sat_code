/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <stddef.h>
#include <math.h>
#include "norad.h"
#include "norad_in.h"

#include <stdio.h>

/* For high satellites,  we do a numerical integration that uses a
rather drastic set of simplifications.  We include the earth,
moon,  and sun,  but with low-precision approximations for the
positions of those last two.  References are to Meeus' _Astronomical
Algorithms_,  2nd edition.  Results are in meters from the center
of the earth,  in ecliptic coordinates of date.

The numerical integration is done using the 'classic' RK4 algorithm,
as described at (for example)

https://en.wikipedia.org/wiki/Runge–Kutta_methods

The solar and lunar positions are computed using Meeus' formulae,  which
are a little computationally intensive.  RK4 has the slight advantage of
requiring us to compute lunar/solar positions only for the steps
themselves and their midpoints.

You probably know that the 'elements' for traditional TLEs are fitted to
the SGP4 and SDP4 models : if you tried to numerically integrate TLEs
using a more sophisticated model,  you'd actually get _worse_ results.
Similarly,  the state vectors for my modified TLEs are integrated using
the following model,  which tries to balance accuracy and speed.  Just as
you shouldn't try to "improve" SGP4/SDP4,  you shouldn't try to "improve"
the following;  it'll only break backward compatibility.  */

static void raw_lunar_solar_position( const double jd, double *lunar_xyzr, double *solar_xyzr)
{
   const double j2000 = 2451545;    /* 1.5 Jan 2000 = JD 2451545 */
   const double t_cen = (jd - j2000) / 36525.;
            /* Mean lunar longitude, (47.1) */
   const double l_prime = 218.3164477 * pi / 180.
                        + (481267.88123421 * pi / 180.) * t_cen;
            /* Lunar mean anomaly,  (47.4) */
   const double m_prime = 134.9633964 * pi / 180.
                        + (477198.8675055 * pi / 180.) * t_cen;
            /* Solar mean longitude, (25.2) */
   const double l_solar = 280.46646 * pi / 180.
                        + (36000.76983 * pi / 180.) * t_cen;
            /* Solar mean anomaly, (47.3)  */
   const double m_solar = 357.5291092 * pi / 180.
                        + (35999.0502909 * pi / 180.) * t_cen;
            /* Lunar mean argument of latitude (47.5) */
   const double f = 93.2720950 * pi / 180.
                        + (483202.0175233 * pi / 180.) * t_cen;
   const double lunar_mean_elong = (297.8501921 * pi / 180.)
                        + (445267.1114034 * pi / 180.) * t_cen;
   const double term2 = 2. * lunar_mean_elong - m_prime;
   const double lunar_lon = l_prime   /* See table 47.A */
                        + (6.288774 * pi / 180.) * sin( m_prime)
                        + (1.274027 * pi / 180.) * sin( term2)
                        + (0.658314 * pi / 180.) * sin( 2. * lunar_mean_elong)
                        + (0.213618 * pi / 180.) * sin( 2. * m_prime)
                        - (0.185166 * pi / 180.) * sin( m_solar)
                        - (0.114332 * pi / 180.) * sin( 2. * f);
   const double lunar_lat = (5.128122 * pi / 180.) * sin( f)
                          + (0.280602 * pi / 180.) * sin( m_prime + f)
                          + (0.277693 * pi / 180.) * sin( m_prime - f)
                          + (0.173237 * pi / 180.) * sin( 2. * lunar_mean_elong - f);
   const double lunar_r = 385000560.       /* in meters */
                         - 20905355. * cos( m_prime)
                         -  3699111. * cos( term2)
                         -  2955968  * cos( 2. * lunar_mean_elong)
                         -   569925  * cos( 2. * m_solar);
   const double solar_ecc = 0.016708634;     /* (25.4) */
   const double solar_lon = l_solar          /* (above (25.5)) */
                        + (1.914602 * pi / 180.) * sin( m_solar);
   const double au_in_meters = 1.495978707e+11;
   const double solar_r = au_in_meters * (1. - solar_ecc * cos( m_solar));
   double tval;

   tval = lunar_r * cos( lunar_lat);
   *lunar_xyzr++ = tval * cos( lunar_lon);
   *lunar_xyzr++ = tval * sin( lunar_lon);
   *lunar_xyzr++ = lunar_r * sin( lunar_lat);
   *lunar_xyzr++ = lunar_r;

   *solar_xyzr++ = solar_r * cos( solar_lon);
   *solar_xyzr++ = solar_r * sin( solar_lon);
   *solar_xyzr++ = 0.;
   *solar_xyzr++ = solar_r;
}

/* For the RK4 integration,  we're frequently asking for the sun and
moon positions at the exact same time we needed for the preceding step.
Caching those positions saves recomputing them. */

void DLL_FUNC lunar_solar_position( const double jd,
                    double *lunar_xyzr, double *solar_xyzr)
{
   static double curr_jd = 0., lunar[4], solar[4];
   size_t i;

   if( curr_jd != jd)
      {
      curr_jd = jd;
      raw_lunar_solar_position( jd, lunar, solar);
      }
   for( i = 0; i < 4; i++)
      {
      if( lunar_xyzr)
         lunar_xyzr[i] = lunar[i];
      if( solar_xyzr)
         solar_xyzr[i] = solar[i];
      }
}

static const double sin_obliq_2000 = 0.397777155931913701597179975942380896684;
static const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;

static void equatorial_to_ecliptic( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 - vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 + vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

static void ecliptic_to_equatorial( double *vect)
{
   double temp;

   temp    = vect[2] * cos_obliq_2000 + vect[1] * sin_obliq_2000;
   vect[1] = vect[1] * cos_obliq_2000 - vect[2] * sin_obliq_2000;
   vect[2] = temp;
}

static void init_high_ephemeris( double *params, const tle_t *tle)
{
   const double *state_vect = &tle->xincl;   /* position at epoch,  in meters */
   size_t i;

   for( i = 0; i < 6; i++)
      params[i] = state_vect[i];
   equatorial_to_ecliptic( params);
   equatorial_to_ecliptic( params + 3);
}

#define c1           params[2]
#define c4           params[3]
#define xnodcf       params[4]
#define t2cof        params[5]
#define deep_arg     ((deep_arg_t *)( params + 10))

void DLL_FUNC SDP4_init( double *params, const tle_t *tle)
{
   init_t init;

   if( tle->ephemeris_type == 'H')
      {
      init_high_ephemeris( params, tle);
      return;
      }
   sxpx_common_init( params, tle, &init, deep_arg);
   deep_arg->sing = sin(tle->omegao);
   deep_arg->cosg = cos(tle->omegao);

   /* initialize Deep() */
   Deep_dpinit( tle, deep_arg);
#ifdef RETAIN_PERTURBATION_VALUES_AT_EPOCH
   /* initialize lunisolar perturbations: */
   deep_arg->t = 0.;                            /* added 30 Dec 2003 */
   deep_arg->solar_lunar_init_flag = 1;
   Deep_dpper( tle, deep_arg);
   deep_arg->solar_lunar_init_flag = 0;
#endif
} /*End of SDP4() initialization */

static inline double vector_len( const double *vect)
{
   double len2 = vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2];

   return( sqrt( len2));
}

/* Input position is in meters,  accel is in m/sec^2 */

static int calc_accel( const double jd, const double *pos, double *accel)
{
   size_t i;
   const double earth_gm = 3.9860044e+14;   /* in m^3/s^2 */
   const double solar_gm = 1.3271243994e+20;   /* m^3/s^2 */
   const double lunar_gm = 4.902798e+12;       /* m^3/s^2 */
   double r = vector_len( pos);
   double accel_factor = -earth_gm / (r * r * r);
   double lunar_xyzr[4], solar_xyzr[4];
   unsigned obj_idx;

   for( i = 0; i < 3; i++)
      accel[i] = accel_factor * pos[i];
   lunar_solar_position( jd, lunar_xyzr, solar_xyzr);
   for( obj_idx = 0; obj_idx < 2; obj_idx++)
      {
      double *opos = (obj_idx ? lunar_xyzr : solar_xyzr);
      double delta[3], d;
      const double gm = (obj_idx ? lunar_gm : solar_gm);
      double accel_factor2;

      accel_factor = gm / (opos[3] * opos[3] * opos[3]);
      for( i = 0; i < 3; i++)
         delta[i] = opos[i] - pos[i];
      d = vector_len( delta);
      accel_factor2 = gm / (d * d * d);
      for( i = 0; i < 3; i++)
         accel[i] -= accel_factor * opos[i] - accel_factor2 * delta[i];
      }
   return( 0);
}

static int calc_state_vector_deriv( const double jd,
                       const double state_vect[6], double deriv[6])
{
   deriv[0] = state_vect[3];
   deriv[1] = state_vect[4];
   deriv[2] = state_vect[5];
   return( calc_accel( jd, state_vect, deriv + 3));
}

/* NOTE: t_since is in minutes,  posn is in km, vel is in km/minutes.
State vector is in meters and m/s.  Hence some conversions...

'high_ephemeris()' does the actual RK4 numerical integration,  using a
simplified model of the earth and moon and a quite basic integration
step size adjustment so that it can take small steps when the object is
close to the earth or moon and larger steps when far away.  As described
above,  any temptation to "improve" the integration should be resisted. */

static int high_ephemeris( double tsince, const tle_t *tle, const double *params,
                                         double *pos, double *vel)
{
   const double meters_per_km = 1000.;
   const double seconds_per_minute = 60.;
   const double seconds_per_day =
               seconds_per_minute * minutes_per_day;     /* a.k.a. 86400 */
   size_t i, j;
   double jd = tle->epoch, state_vect[6];

   for( i = 0; i < 6; i++)
      state_vect[i] = params[i];
   tsince /= minutes_per_day;       /* input was in minutes;  days are */
   while( tsince)                   /* more convenient hereforth       */
      {
      double dt = tsince, dt_in_seconds;
      double max_step = 1.;
      double kvects[4][6];

      calc_state_vector_deriv( jd, state_vect, kvects[0]);
      for( j = 3; j < 6; j++)
         if( max_step > 1e-3 / fabs( kvects[0][j]))
            max_step = 1e-3 / fabs( kvects[0][j]);
      if( max_step < 1e-5)
         max_step = 1e-5;
      if( tsince > max_step)
         dt = max_step;
      else if( tsince < -max_step)
         dt = -max_step;
      dt_in_seconds = dt * seconds_per_day;
      for( j = 1; j < 4; j++)
         {
         const double step = (j == 3 ? dt_in_seconds : dt_in_seconds * .5);
         double tstate[6];

         for( i = 0; i < 6; i++)
            tstate[i] = state_vect[i] + step * kvects[j - 1][i];
         calc_state_vector_deriv( jd + (j == 3 ? dt : dt / 2.),
                        tstate, kvects[j]);
         }

      for( i = 0; i < 6; i++)
         state_vect[i] += (dt_in_seconds / 6.) *
               (kvects[0][i] + 2. * (kvects[1][i] + kvects[2][i]) + kvects[3][i]);
      jd += dt;
      tsince -= dt;
      }
   for( i = 0; i < 3; i++)
      {
      pos[i] = state_vect[i];
      vel[i] = state_vect[i + 3];
      }
   ecliptic_to_equatorial( vel);
   ecliptic_to_equatorial( pos);
                  /* Now,  cvt meters to km, meters/second to km/minute: */
   for( i = 0; i < 3; i++)
      {
      pos[i] /= meters_per_km;
      vel[i] *= seconds_per_minute / meters_per_km;
      }
   return( 0);
}

int DLL_FUNC SDP4( const double tsince, const tle_t *tle, const double *params,
                                         double *pos, double *vel)
{
  double
      a, tempa, tsince_squared,
      xl, xnoddf;

   if( tle->ephemeris_type == 'H')
      {
      double unused_vel[3];

      return( high_ephemeris( tsince, tle, params, pos, (vel ? vel : unused_vel)));
      }
  /* Update for secular gravity and atmospheric drag */
  deep_arg->omgadf = tle->omegao + deep_arg->omgdot * tsince;
  xnoddf = tle->xnodeo + deep_arg->xnodot * tsince;
  tsince_squared = tsince*tsince;
  deep_arg->xnode = xnoddf + xnodcf * tsince_squared;
  deep_arg->xn = deep_arg->xnodp;

  /* Update for deep-space secular effects */
  deep_arg->xll = tle->xmo + deep_arg->xmdot * tsince;
  deep_arg->t = tsince;

  Deep_dpsec( tle, deep_arg);

  tempa = 1-c1*tsince;
  if( deep_arg->xn < 0.)
     return( SXPX_ERR_NEGATIVE_XN);
  a = pow(xke/deep_arg->xn,two_thirds)*tempa*tempa;
  deep_arg->em -= tle->bstar*c4*tsince;

  /* Update for deep-space periodic effects */
  deep_arg->xll += deep_arg->xnodp * t2cof * tsince_squared;

  Deep_dpper( tle, deep_arg);

            /* Keeping xinc positive is not really necessary,  unless        */
            /* you're displaying elements and dislike negative inclinations. */
#ifdef KEEP_INCLINATION_POSITIVE
  if (deep_arg->xinc < 0.)       /* Begin April 1983 errata correction: */
     {
     deep_arg->xinc = -deep_arg->xinc;
     deep_arg->sinio = -deep_arg->sinio;
     deep_arg->xnode += pi;
     deep_arg->omgadf -= pi;
     }                          /* End April 1983 errata correction. */
#endif

  xl = deep_arg->xll + deep_arg->omgadf + deep_arg->xnode;
               /* Dundee change:  Reset cosio,  sinio for new xinc: */
  deep_arg->cosio = cos( deep_arg->xinc);
  deep_arg->sinio = sin( deep_arg->xinc);

  return( sxpx_posn_vel( deep_arg->xnode, a, deep_arg->em, deep_arg->cosio,
                deep_arg->sinio, deep_arg->xinc, deep_arg->omgadf,
                xl, pos, vel));
} /* SDP4 */
