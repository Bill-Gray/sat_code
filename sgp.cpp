/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <math.h>
#include "norad.h"
#include "norad_in.h"

#define ao       params[0]
#define qo       params[1]
#define xlo      params[2]
#define d1o      params[3]
#define d2o      params[4]
#define d3o      params[5]
#define d4o      params[6]
#define omgdt    params[7]
#define xnodot   params[8]
#define c5       params[9]
#define c6       params[10]

void DLL_FUNC SGP_init( double *params, const tle_t *tle)
{
   double c1, c2, c3, c4, r1, cosio, sinio, a1, d1, po, po2no;

   c1 = ck2*1.5;
   c2 = ck2/4.;
   c3 = ck2/2.;
   r1 = ae;
   c4 = xj3*(r1*(r1*r1))/(ck2*4.);
   cosio = cos(tle->xincl);
   sinio = sin(tle->xincl);
   a1 = pow( xke / tle->xno, two_thirds);
   d1 = c1/a1/a1*(cosio*3.*cosio-1.)/pow( 1.-tle->eo*tle->eo, 1.5);
   ao = a1*(1.-d1*.33333333333333331-d1*d1-d1*
                            1.654320987654321*d1*d1);
   po = ao*(1.-tle->eo*tle->eo);
   qo = ao*(1.-tle->eo);
   xlo = tle->xmo+tle->omegao+tle->xnodeo;
   d1o = c3*sinio*sinio;
   d2o = c2*(cosio*7.*cosio-1.);
   d3o = c1*cosio;
   d4o = d3o*sinio;
   po2no = tle->xno/(po*po);
   omgdt = c1*po2no*(cosio*5.*cosio-1.);
   xnodot = d3o*-2.*po2no;
   c5 = c4*.5*sinio*(cosio*5.+3.)/(cosio+1.);
   c6 = c4*sinio;
}


int DLL_FUNC SGP( const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel)
{
  double
    temp, rdot, cosu, sinu, cos2u, sin2u, a, e,
    p, rr, u, ecose, esine, omgas, cosik, xinck,
    sinik, axnsl, aynsl,
    sinuk, rvdot, cosuk, coseo1, sineo1, pl,
    rk, uk, xl, su, ux, uy, uz, vx, vy, vz, pl2,
    xnodek, cosnok, xnodes, el2, eo1, r1, sinnok,
    xls, xmx, xmy, tem2, tem5;
   const double chicken_factor_on_eccentricity = 1.e-6;

   int i, rval = 0;

    /* Update for secular gravity and atmospheric drag */
   a = tle->xno+(tle->xndt2o*2.+tle->xndd6o*3.*tsince)*tsince;
   if( a < 0.)
      rval = SXPX_ERR_NEGATIVE_MAJOR_AXIS;
   e = e6a;
   if( e > 1. - chicken_factor_on_eccentricity)
      rval = SXPX_ERR_NEARLY_PARABOLIC;
   if( rval)
      {
      for( i = 0; i < 3; i++)
         {
         pos[i] = 0.;
         if( vel)
            vel[i] = 0.;
         }
      return( rval);
      }
   a = ao * pow( tle->xno / a, two_thirds);
   if( a * (1. - e) < 1. && a * (1. + e) < 1.)   /* entirely within earth */
      rval = SXPX_WARN_ORBIT_WITHIN_EARTH;     /* remember, e can be negative */
   if( a * (1. - e) < 1. || a * (1. + e) < 1.)   /* perigee within earth */
      rval = SXPX_WARN_PERIGEE_WITHIN_EARTH;
   if (a > qo) e = 1.-qo/a;
   p = a*(1.-e*e);
   xnodes = tle->xnodeo+xnodot*tsince;
   omgas = tle->omegao+omgdt*tsince;
   r1 = xlo+(tle->xno+omgdt+xnodot+
       (tle->xndt2o+tle->xndd6o*tsince)*tsince)*tsince;
   xls = FMod2p(r1);

   /* Long period periodics */
   axnsl = e*cos(omgas);
   aynsl = e*sin(omgas)-c6/p;
   r1 = xls-c5/p*axnsl;
   xl = FMod2p(r1);

   /* Solve Kepler's equation */
   r1 = xl-xnodes;
   u = FMod2p(r1);
   eo1 = u;
   tem5 = 1.;

   i = 0;
   do
     {
       sineo1 = sin(eo1);
       coseo1 = cos(eo1);
       if (fabs(tem5) < e6a) break;
       tem5 = 1.-coseo1*axnsl-sineo1*aynsl;
       tem5 = (u-aynsl*coseo1+axnsl*sineo1-eo1)/tem5;
       tem2 = fabs(tem5);
       if (tem2 > 1.) tem5 = tem2/tem5;
       eo1 += tem5;
     }
   while(i++ < 10);

   /* Short period preliminary quantities */
   ecose = axnsl*coseo1+aynsl*sineo1;
   esine = axnsl*sineo1-aynsl*coseo1;
   el2 = axnsl*axnsl+aynsl*aynsl;
   pl = a*(1.-el2);
   pl2 = pl*pl;
   rr = a*(1.-ecose);
   rdot = xke*sqrt(a)/rr*esine;
   rvdot = xke*sqrt(pl)/rr;
   temp = esine/(sqrt(1.-el2)+1.);
   sinu = a/rr*(sineo1-aynsl-axnsl*temp);
   cosu = a/rr*(coseo1-axnsl+aynsl*temp);
   su = atan2(sinu, cosu);

   /* Update for short periodics */
   sin2u = (cosu+cosu)*sinu;
   cos2u = 1.-2.*sinu*sinu;
   rk = rr+d1o/pl*cos2u;
   uk = su-d2o/pl2*sin2u;
   xnodek = xnodes+d3o*sin2u/pl2;
   xinck = tle->xincl+d4o/pl2*cos2u;

   /* Orientation vectors */
   sinuk = sin(uk);
   cosuk = cos(uk);
   sinnok = sin(xnodek);
   cosnok = cos(xnodek);
   sinik = sin(xinck);
   cosik = cos(xinck);
   xmx = -sinnok*cosik;
   xmy = cosnok*cosik;
   ux = xmx*sinuk+cosnok*cosuk;
   uy = xmy*sinuk+sinnok*cosuk;
   uz = sinik*sinuk;
   vx = xmx*cosuk-cosnok*sinuk;
   vy = xmy*cosuk-sinnok*sinuk;
   vz = sinik*cosuk;

   /* Position and velocity */
   pos[0] = rk*ux*earth_radius_in_km;
   pos[1] = rk*uy*earth_radius_in_km;
   pos[2] = rk*uz*earth_radius_in_km;
   if( vel)
      {
      vel[0] = (rdot*ux + rvdot * vx)*earth_radius_in_km;
      vel[1] = (rdot*uy + rvdot * vy)*earth_radius_in_km;
      vel[2] = (rdot*uz + rvdot * vz)*earth_radius_in_km;
      }
   return( rval);
} /* SGP */
