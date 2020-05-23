/* Copyright (C) 2018, Project Pluto.  See LICENSE. */

#include <math.h>
#include <assert.h>
#include "norad.h"
#include "norad_in.h"

/* params[1] and [6]-[9] were used in earlier implementations,  but are
   now unused */

#define c2           params[0]
#define c1           params[2]
#define c4           params[3]
#define xnodcf       params[4]
#define t2cof        params[5]

void sxpall_common_init( const tle_t *tle, deep_arg_t *deep_arg)
{
   const double a1 = pow(xke / tle->xno, two_thirds);  /* in Earth radii */
   double del1, ao, delo, tval;

   /* Recover original mean motion (xnodp) and   */
   /* semimajor axis (aodp) from input elements. */
   deep_arg->cosio = cos( tle->xincl);
   deep_arg->cosio2 = deep_arg->cosio * deep_arg->cosio;
   deep_arg->eosq = tle->eo*tle->eo;
   deep_arg->betao2 = 1-deep_arg->eosq;
   deep_arg->betao = sqrt(deep_arg->betao2);
   tval = 1.5 * ck2 * (3. * deep_arg->cosio2 - 1.) / (deep_arg->betao * deep_arg->betao2);
   del1 = tval / (a1 * a1);
   ao = a1 * (1. - del1 * (1. / 3. + del1 * ( 1. + 134./81. * del1)));
   delo = tval / (ao * ao);
   deep_arg->xnodp = tle->xno / (1+delo);   /* in radians/minute */
   deep_arg->aodp = ao / (1-delo);
}

void sxpx_common_init( double *params, const tle_t *tle,
                                  init_t *init, deep_arg_t *deep_arg)
{
   double
         eeta, etasq, perige, pinv, pinvsq,
         psisq, qoms24, temp1, temp2, temp3,
         cosio4, tsi_squared, x3thm1, xhdot1;

   sxpall_common_init( tle, deep_arg);
   x3thm1 = 3. * deep_arg->cosio2 - 1.;
   /* For perigee below 156 km, the values */
   /* of s and qoms2t are altered.         */
   init->s4 = s_const;
   qoms24 = qoms2t;
   perige = (deep_arg->aodp * (1-tle->eo) - ae) * earth_radius_in_km;
   if( perige < 156.)
      {
      double temp_val, temp_val_squared;

      if(perige <= 98.)
         init->s4 = 20;
      else
         init->s4 = perige-78.;
      temp_val = (120. - init->s4) * ae / earth_radius_in_km;
      temp_val_squared = temp_val * temp_val;
      qoms24 = temp_val_squared * temp_val_squared;
      init->s4 = init->s4 / earth_radius_in_km + ae;
      }  /* End of if(perige <= 156) */

   pinv = 1. / (deep_arg->aodp * deep_arg->betao2);
   pinvsq = pinv * pinv;
   init->tsi = 1. / (deep_arg->aodp - init->s4);
   init->eta = deep_arg->aodp*tle->eo*init->tsi;
   etasq = init->eta*init->eta;
   eeta = tle->eo*init->eta;
   psisq = fabs(1-etasq);
   tsi_squared = init->tsi * init->tsi;
   init->coef = qoms24 * tsi_squared * tsi_squared;
   init->coef1 = init->coef / pow(psisq,3.5);
   c2 = init->coef1 * deep_arg->xnodp * (deep_arg->aodp*(1+1.5*etasq+eeta*
   (4+etasq))+0.75*ck2*init->tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
   c1 = tle->bstar*c2;
   deep_arg->sinio = sin(tle->xincl);
   c4 = 2*deep_arg->xnodp*init->coef1*deep_arg->aodp*deep_arg->betao2*
        (init->eta*(2+0.5*etasq)+tle->eo*(0.5+2*etasq)-2*ck2*init->tsi/
        (deep_arg->aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*
        (1.5-0.5*eeta))+0.75*(1. - deep_arg->cosio2) *(2*etasq-eeta*(1+etasq))*
        cos(2*tle->omegao)));
   cosio4 = deep_arg->cosio2 * deep_arg->cosio2;
   temp1 = 3*ck2*pinvsq*deep_arg->xnodp;
   temp2 = temp1 * ck2 * pinvsq;
   temp3 = 1.25 * ck4 * pinvsq * pinvsq * deep_arg->xnodp;
   deep_arg->xmdot = deep_arg->xnodp
            + temp1 * deep_arg->betao* x3thm1 / 2.
            + temp2 * deep_arg->betao*
                    (13-78*deep_arg->cosio2+137*cosio4) / 16.;
   deep_arg->omgdot = -temp1 * (1. - 5 * deep_arg->cosio2) / 2.
              + temp2 * (7-114*deep_arg->cosio2+395*cosio4) / 16.
              + temp3 * (3-36*deep_arg->cosio2+49*cosio4);
   xhdot1 = -temp1*deep_arg->cosio;
   deep_arg->xnodot = xhdot1+(temp2*(4-19*deep_arg->cosio2) / 2.
           + 2*temp3*(3-7*deep_arg->cosio2))*deep_arg->cosio;
   xnodcf = 3.5*deep_arg->betao2*xhdot1*c1;
   t2cof = 1.5*c1;
}

inline double centralize_angle( const double ival)
{
   double rval = fmod( ival, twopi);

   if( rval > pi)
      rval -= twopi;
   else if( rval < - pi)
      rval += twopi;
   return( rval);
}

#define MAX_KEPLER_ITER 10

int sxpx_posn_vel( const double xnode, const double a, const double ecc,
      const double cosio, const double sinio,
      const double xincl, const double omega,
      const double xl, double *pos, double *vel)
{
  /* Long period periodics */
   const double axn = ecc*cos(omega);
   double temp = 1/(a*(1.-ecc*ecc));
   const double xlcof = .125 * a3ovk2 * sinio * (3+5*cosio)/ (1. + cosio);
   const double aycof = 0.25 * a3ovk2 * sinio;
   const double xll = temp*xlcof*axn;
   const double aynl = temp*aycof;
   const double xlt = xl+xll;
   const double ayn = ecc*sin(omega)+aynl;
   const double elsq = axn*axn+ayn*ayn;
   const double capu = centralize_angle( xlt - xnode);
   const double chicken_factor_on_eccentricity = 1.e-6;
   double epw = capu;
   double temp1, temp2;
   double ecosE, esinE, pl, r;
   double betal;
   double u, sinu, cosu, sin2u, cos2u;
   double rk, uk, xnodek, xinck;
   double sinuk, cosuk, sinik, cosik, sinnok, cosnok, xmx, xmy;
   double sinEPW, cosEPW;
   double ux, uy, uz;
   int i, rval = 0;

/* Dundee changes:  items dependent on cosio get recomputed: */
   const double cosio_squared = cosio * cosio;
   const double x3thm1 = 3.0 * cosio_squared - 1.0;
   const double sinio2 = 1.0 - cosio_squared;
   const double x7thm1 = 7.0 * cosio_squared - 1.0;

                /* Added 29 Mar 2003,  modified 26 Sep 2006:  extremely    */
                /* decayed satellites can end up "orbiting" within the     */
                /* earth.  Eventually,  the semimajor axis becomes zero,   */
                /* then negative.  In that case,  or if the orbit is near  */
                /* to parabolic,  we zero the posn/vel and quit.  If the   */
                /* object has a perigee or apogee indicating a crash,  we  */
                /* just flag it.  Revised 28 Oct 2006.                     */

   if( a < 0.)
      rval = SXPX_ERR_NEGATIVE_MAJOR_AXIS;
   if( elsq > 1. - chicken_factor_on_eccentricity)
      rval = SXPX_ERR_NEARLY_PARABOLIC;
   for( i = 0; i < 3; i++)
      {
      pos[i] = 0.;
      if( vel)
         vel[i] = 0.;
      }
   if( rval)
      return( rval);
   if( a * (1. - ecc) < 1. && a * (1. + ecc) < 1.)   /* entirely within earth */
      rval = SXPX_WARN_ORBIT_WITHIN_EARTH;     /* remember, e can be negative */
   if( a * (1. - ecc) < 1. || a * (1. + ecc) < 1.)   /* perigee within earth */
      rval = SXPX_WARN_PERIGEE_WITHIN_EARTH;
  /* Solve Kepler's' Equation */
   for( i = 0; i < MAX_KEPLER_ITER; i++)
      {
      const double newton_raphson_epsilon = 1e-12;
      double f, fdot, delta_epw;
      int do_second_order_newton_raphson = 1;

      sinEPW = sin( epw);
      cosEPW = cos( epw);
      ecosE = axn * cosEPW + ayn * sinEPW;
      esinE = axn * sinEPW - ayn * cosEPW;
      f = capu - epw + esinE;
      if (fabs(f) < newton_raphson_epsilon) break;
      fdot = 1. - ecosE;
      delta_epw = f / fdot;
      if( !i)
         {
         const double max_newton_raphson = 1.25 * fabs( ecc);

         do_second_order_newton_raphson = 0;
         if( delta_epw > max_newton_raphson)
            delta_epw = max_newton_raphson;
         else if( delta_epw < -max_newton_raphson)
            delta_epw = -max_newton_raphson;
         else
            do_second_order_newton_raphson = 1;
         }
      if( do_second_order_newton_raphson)
         delta_epw = f / (fdot + 0.5*esinE*delta_epw);
                             /* f/(fdot - 0.5*fdotdot * f / fdot) */
      epw += delta_epw;
      }

   if( i == MAX_KEPLER_ITER)
      return( SXPX_ERR_CONVERGENCE_FAIL);

  /* Short period preliminary quantities */
   temp = 1-elsq;
   pl = a*temp;
   r = a*(1-ecosE);
   temp2 = a / r;
   betal = sqrt(temp);
   temp = esinE/(1+betal);
   cosu = temp2 * (cosEPW - axn + ayn * temp);
   sinu = temp2 * (sinEPW - ayn - axn * temp);
   u = atan2( sinu, cosu);
   sin2u = 2*sinu*cosu;
   cos2u = 2*cosu*cosu-1;
   temp1 = ck2 / pl;
   temp2 = temp1 / pl;

  /* Update for short periodics */
   rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*sinio2*cos2u;
   uk = u-0.25*temp2*x7thm1*sin2u;
   xnodek = xnode+1.5*temp2*cosio*sin2u;
   xinck = xincl+1.5*temp2*cosio*sinio*cos2u;

  /* Orientation vectors */
   sinuk = sin(uk);
   cosuk = cos(uk);
   sinik = sin(xinck);
   cosik = cos(xinck);
   sinnok = sin(xnodek);
   cosnok = cos(xnodek);
   xmx = -sinnok*cosik;
   xmy = cosnok*cosik;
   ux = xmx*sinuk+cosnok*cosuk;
   uy = xmy*sinuk+sinnok*cosuk;
   uz = sinik*sinuk;

  /* Position and velocity */
   pos[0] = rk * ux * earth_radius_in_km;
   pos[1] = rk * uy * earth_radius_in_km;
   pos[2] = rk * uz * earth_radius_in_km;
   if( vel)
      {
      const double rdot = xke * sqrt(a) * esinE / r;
      const double rfdot = xke * sqrt(pl) / r;
      const double xn = xke / (a * sqrt(a));
      const double rdotk = rdot - xn * temp1 * sinio2 * sin2u;
      const double rfdotk = rfdot + xn * temp1 * (sinio2 * cos2u + 1.5 * x3thm1);
      const double vx = xmx * cosuk - cosnok * sinuk;
      const double vy = xmy * cosuk - sinnok * sinuk;
      const double vz = sinik*cosuk;

      vel[0] = (rdotk * ux + rfdotk * vx) * earth_radius_in_km;
      vel[1] = (rdotk * uy + rfdotk * vy) * earth_radius_in_km;
      vel[2] = (rdotk * uz + rfdotk * vz) * earth_radius_in_km;
      }
   return( rval);
} /*SGP4*/
