/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <math.h>
#include <stdio.h>
#include "norad.h"
#include "norad_in.h"

#define c1           params[2]
#define c4           params[3]
#define xnodcf       params[4]
#define t2cof        params[5]
#define p_aodp         params[10]
#define p_cosio        params[11]
#define p_sinio        params[12]
#define p_omgdot       params[13]
#define p_xmdot        params[14]
#define p_xnodot       params[15]
#define p_xnodp        params[16]
#define c5             params[17]
#define d2             params[18]
#define d3             params[19]
#define d4             params[20]
#define delmo          params[21]
#define p_eta          params[22]
#define omgcof         params[23]
#define sinmo          params[24]
#define t3cof          params[25]
#define t4cof          params[26]
#define t5cof          params[27]
#define xmcof          params[28]
#define simple_flag *((int *)( params + 29))
#define MINIMAL_E    1.e-4
#define ECC_EPS      1.e-6     /* Too low for computing further drops. */

void DLL_FUNC SGP4_init( double *params, const tle_t *tle)
{
   deep_arg_t deep_arg;
   init_t init;
   double eeta, etasq;

   sxpx_common_init( params, tle, &init, &deep_arg);
   p_aodp =   deep_arg.aodp;
   p_cosio =  deep_arg.cosio;
   p_sinio =  deep_arg.sinio;
   p_omgdot = deep_arg.omgdot;
   p_xmdot =  deep_arg.xmdot;
   p_xnodot = deep_arg.xnodot;
   p_xnodp =  deep_arg.xnodp;
   p_eta = deep_arg.aodp*tle->eo*init.tsi;
// p_eta = init.eta;

   eeta = tle->eo*p_eta;
   /* For perigee less than 220 kilometers, the "simple" flag is set */
   /* and the equations are truncated to linear variation in sqrt a  */
   /* and quadratic variation in mean anomaly.  Also, the c3 term,   */
   /* the delta omega term, and the delta m term are dropped.        */
   simple_flag = ((p_aodp*(1-tle->eo)/ae) < (220./earth_radius_in_km+ae));
   if( !simple_flag)
      {
      const double c1sq = c1*c1;
      double temp;

      simple_flag = 0;
      delmo = 1. + p_eta * cos(tle->xmo);
      delmo *= delmo * delmo;
      d2 = 4*p_aodp*init.tsi*c1sq;
      temp = d2*init.tsi*c1/3;
      d3 = (17*p_aodp+init.s4)*temp;
      d4 = 0.5*temp*p_aodp*init.tsi*(221*p_aodp+31*init.s4)*c1;
      t3cof = d2+2*c1sq;
      t4cof = 0.25*(3*d3+c1*(12*d2+10*c1sq));
      t5cof = 0.2*(3*d4+12*c1*d3+6*d2*d2+15*c1sq*(2*d2+c1sq));
      sinmo = sin(tle->xmo);
      if( tle->eo < MINIMAL_E)
         omgcof = xmcof = 0.;
      else
         {
         const double c3 =
              init.coef * init.tsi * a3ovk2 * p_xnodp * ae * p_sinio / tle->eo;

         xmcof = -two_thirds * init.coef * tle->bstar * ae / eeta;
         omgcof = tle->bstar*c3*cos(tle->omegao);
         }
      } /* End of if (isFlagClear(SIMPLE_FLAG)) */
   etasq = p_eta * p_eta;
   c5 = 2*init.coef1*p_aodp * deep_arg.betao2*(1+2.75*(etasq+eeta)+eeta*etasq);
} /* End of SGP4() initialization */

int DLL_FUNC SGP4( const double tsince, const tle_t *tle, const double *params,
                                                    double *pos, double *vel)
{
  double
        a, e, omega, omgadf,
        temp, tempa, tempe, templ, tsq,
        xl, xmdf, xmp, xnoddf, xnode;

  /* Update for secular gravity and atmospheric drag. */
  xmdf = tle->xmo+p_xmdot*tsince;
  omgadf = tle->omegao+p_omgdot*tsince;
  xnoddf = tle->xnodeo+p_xnodot*tsince;
  omega = omgadf;
  xmp = xmdf;
  tsq = tsince*tsince;
  xnode = xnoddf+xnodcf*tsq;
  tempa = 1-c1*tsince;
  tempe = tle->bstar*c4*tsince;
  templ = t2cof*tsq;
  if( !simple_flag)
    {
      const double delomg = omgcof*tsince;
      double delm = 1. + p_eta * cos(xmdf);
      double tcube, tfour;

      delm = xmcof * (delm * delm * delm - delmo);
      temp = delomg+delm;
      xmp = xmdf+temp;
      omega = omgadf-temp;
      tcube = tsq*tsince;
      tfour = tsince*tcube;
      tempa = tempa-d2*tsq-d3*tcube-d4*tfour;
      tempe = tempe+tle->bstar*c5*(sin(xmp)-sinmo);
      templ = templ+t3cof*tcube+tfour*(t4cof+tsince*t5cof);
    }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

  a = p_aodp*tempa*tempa;
  e = tle->eo-tempe;
         /* A highly arbitrary lower limit on e,  of 1e-6: */
  if( e < ECC_EPS)
     e = ECC_EPS;
  xl = xmp+omega+xnode+p_xnodp*templ;
  if( tempa < 0.)       /* force negative a,  to indicate error condition */
     a = -a;
  return( sxpx_posn_vel( xnode, a, e, p_cosio, p_sinio, tle->xincl,
                                          omega, xl, pos, vel));
} /*SGP4*/
