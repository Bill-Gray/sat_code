/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <math.h>
#include "observe.h"

/*  Assorted functions useful in conjunction with the satellite code
   library for determining where the _observer_,  as well as the _target_,
   happens to be.  Combine the two positions,  and you can get the
   distance/RA/dec of the target as seen by the observer. */

#define PI 3.141592653589793238462643383279
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)

      /* function for Greenwich sidereal time,  ripped from 'deep.cpp' */

static inline double ThetaG( double jd)
{
  /* Reference:  The 1992 Astronomical Almanac, page B6. */
  const double omega_E = 1.00273790934;
                   /* Earth rotations per sidereal day (non-constant) */
  const double UT = fmod( jd + .5, 1.);
  const double seconds_per_day = 86400.;
  const double jd_2000 = 2451545.0;   /* 1.5 Jan 2000 = JD 2451545. */
  double t_cen, GMST, rval;

  t_cen = (jd - UT - jd_2000) / 36525.;
  GMST = 24110.54841 + t_cen * (8640184.812866 + t_cen *
                           (0.093104 - t_cen * 6.2E-6));
  GMST = fmod( GMST + seconds_per_day * omega_E * UT, seconds_per_day);
  if( GMST < 0.)
     GMST += seconds_per_day;
  rval = 2. * PI * GMST / seconds_per_day;

  return( rval);
} /*Function thetag*/

void DLL_FUNC observer_cartesian_coords( const double jd, const double lon,
              const double rho_cos_phi, const double rho_sin_phi,
              double *vect)
{
   const double angle = lon + ThetaG( jd);

   *vect++ = cos( angle) * rho_cos_phi * EARTH_MAJOR_AXIS / 1000.;
   *vect++ = sin( angle) * rho_cos_phi * EARTH_MAJOR_AXIS / 1000.;
   *vect++ = rho_sin_phi               * EARTH_MAJOR_AXIS / 1000.;
}

void DLL_FUNC earth_lat_alt_to_parallax( const double lat,
                    const double ht_in_meters,
                    double *rho_cos_phi, double *rho_sin_phi)
{
   const double u = atan( sin( lat) * EARTH_AXIS_RATIO / cos( lat));

   *rho_sin_phi = EARTH_AXIS_RATIO * sin( u) +
                           (ht_in_meters / EARTH_MAJOR_AXIS) * sin( lat);
   *rho_cos_phi = cos( u) + (ht_in_meters / EARTH_MAJOR_AXIS) * cos( lat);
}

void DLL_FUNC get_satellite_ra_dec_delta( const double *observer_loc,
                                 const double *satellite_loc, double *ra,
                                 double *dec, double *delta)
{
   double vect[3], dist2 = 0.;
   int i;

   for( i = 0; i < 3; i++)
      {
      vect[i] = satellite_loc[i] - observer_loc[i];
      dist2 += vect[i] * vect[i];
      }
   *delta = sqrt( dist2);
   *ra = atan2( vect[1], vect[0]);
   if( *ra < 0.)
      *ra += PI + PI;
   *dec = asin( vect[2] / *delta);
}

/* Formulae from Meeus' _Astronomical Algorithms_ for approximate precession.
More than accurate enough for our purposes.  */

static void precess( const double t_centuries, double *ra, double *dec)
{
   const double m = (3.07496 + .00186 * t_centuries / 2.) * (PI / 180.) / 240.;
   const double n = (1.33621 - .00057 * t_centuries / 2.) * (PI / 180.) / 240.;
   const double ra_rate  = m + n * sin( *ra) * tan( *dec);
   const double dec_rate = n * cos( *ra);

   *ra -= t_centuries * ra_rate * 100.;
   *dec -= t_centuries * dec_rate * 100.;
}

void DLL_FUNC epoch_of_date_to_j2000( const double jd, double *ra, double *dec)
{
   const double t_centuries = (jd - 2451545.) / 36525.;

   precess( t_centuries, ra, dec);
}

void DLL_FUNC j2000_to_epoch_of_date( const double jd, double *ra, double *dec)
{
   const double t_centuries = (jd - 2451545.) / 36525.;

   precess( -t_centuries, ra, dec);
}
