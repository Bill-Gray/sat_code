/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <math.h>
#include <string.h>
#include "norad.h"
#include "norad_in.h"

#define tthmun              params[0]
#define sini2               params[1]
#define cosi2               params[2]
#define unm5th              params[3]
#define unmth2              params[4]
#define xmdt1               params[5]
#define xgdt1               params[6]
#define xhdt1               params[7]
#define xndt                params[8]
#define edot                params[9]
#define ed                  params[10]
#define gamma               params[11]
#define omgdt               params[12]
#define ovgpp               params[13]
#define pp                  params[14]
#define qq                  params[15]
#define sini                params[16]
#define cosi                params[17]
#define cosio_2             params[18]
#define xlldot              params[19]
#define xnd                 params[20]
#define xnodot_             params[21]
#define xnodp_              params[22]
#define simple_flag       *((int *)( params + 23))

void sxpall_common_init( const tle_t *tle, deep_arg_t *deep_arg);

void sxp8_common_init( double *params, const tle_t *tle, deep_arg_t *deep_arg)
{
   const double half_inclination = tle->xincl*.5;
   const double cosio4 = deep_arg->cosio2 * deep_arg->cosio2;
   double po, pom2, pardt1, pardt2, pardt4;

   deep_arg->sing = sin( tle->omegao);
   deep_arg->cosg = cos( tle->omegao);
   sini2 = sin( half_inclination);
   cosi2 = cos( half_inclination);
   tthmun = deep_arg->cosio2 * 3. - 1.;
   unm5th = 1.-deep_arg->cosio2 * 5.;
   unmth2 = 1.-deep_arg->cosio2;
   po = deep_arg->aodp * deep_arg->betao2;
   pom2 = 1./(po*po);
   pardt1 = 3. * ck2 * pom2 * deep_arg->xnodp;
   pardt2 = pardt1 * ck2 * pom2;
   pardt4 = ck4 * 1.25 * pom2 * pom2 * deep_arg->xnodp;
   xmdt1 = .5 * pardt1 * deep_arg->betao * tthmun;
   xgdt1 = -.5 * pardt1 * unm5th;
   xhdt1 = -pardt1 * deep_arg->cosio;
   deep_arg->xmdot = deep_arg->xnodp+xmdt1+pardt2*.0625*deep_arg->betao*
               (13.-deep_arg->cosio2*78.+cosio4*137.);
   deep_arg->omgdot = xgdt1+pardt2*.0625*(7.-deep_arg->cosio2*
                     114.+cosio4*395.)+pardt4*(3.-deep_arg->cosio2*
                     36.+cosio4*49.);
   deep_arg->xnodot = xhdt1+(pardt2*.5*(4.-deep_arg->cosio2*19.)+pardt4*
      2.*(3.-deep_arg->cosio2*7.))*deep_arg->cosio;
}

void DLL_FUNC SGP8_init( double *params, const tle_t *tle)
{
   const double rho = .15696615;
   const double b = tle->bstar*2./rho;
   double
         alpha2, b1, b2, b3, c0,
         c1, c4, c5, cos2g, d1, d2, d3, d4,
         d5, eeta, eta, eta2, eddot,  etdt,
         po, psim2, r1, tsi, xndtn;
   deep_arg_t deep_arg;

   sxpall_common_init( tle, &deep_arg);
   sxp8_common_init( params, tle, &deep_arg);
   sini = sin( tle->xincl);
   cosi = deep_arg.cosio;
   cosio_2 = deep_arg.cosio2;

   /* Initialization */
   xnodp_ = deep_arg.xnodp;
   xlldot = deep_arg.xmdot;
   omgdt = deep_arg.omgdot;
   xnodot_ = deep_arg.xnodot;
   po = deep_arg.aodp * deep_arg.betao2;
   tsi = 1./(po-s_const);
   eta = tle->eo*s_const*tsi;
   eta2 = eta * eta;
   psim2 = (r1 = 1./(1.-eta2), fabs(r1));
   alpha2 = deep_arg.eosq+1.;
   eeta = tle->eo*eta;
   cos2g = deep_arg.cosg * deep_arg.cosg * 2. - 1.;
   d5 = tsi*psim2;
   d1 = d5/po;
   d2 = eta2*(eta2*4.5+36.)+12.;
   d3 = eta2*(eta2*2.5+15.);
   d4 = eta*(eta2*3.75+5.);
   b1 = ck2*tthmun;
   b2 = -ck2*unmth2;
   b3 = a3ovk2*sini;
   r1 = tsi, r1 *= r1;
   c0 = b*.5*rho*qoms2t*xnodp_*deep_arg.aodp*(r1*r1)*
                pow(psim2, 3.5)/sqrt(alpha2);
   r1 = alpha2;
   c1 = xnodp_*1.5*(r1*r1)*c0;
   c4 = d1*d3*b2;
   c5 = d5*d4*b3;
   xndt = c1*(eta2*(deep_arg.eosq*34.+3.)+2.+eeta*5.*(eta2+4.)
     +deep_arg.eosq*8.5+d1*d2*b1+c4*cos2g+c5*deep_arg.sing);
   xndtn = xndt/xnodp_;

   /* If drag is very small, the isimp flag is set and the */
   /* equations are truncated to linear variation in mean  */
   /* motion and quadratic variation in mean anomaly       */
   r1 = xndtn * minutes_per_day;
   if( fabs(r1) > .00216)
      {
      const double d6 = eta*(eta2*22.5+30.);
      const double d7 = eta*(eta2*12.5+5.);
      const double d8 = eta2*(eta2+6.75)+1.;
      const double d9 = eta*(deep_arg.eosq*68.+6.)+tle->eo*(eta2*15.+20.);
      const double d10 = eta*5.*(eta2+4.)+tle->eo*(eta2*68.+17.);
      const double d11 = eta*(eta2*18.+72.);
      const double d12 = eta*(eta2*10.+30.);
      const double d13 = eta2*11.25+5.;
      const double d20 = two_thirds*.5*xndtn;
      const double c8 = d1*d7*b2;
      const double c9 = d5*d8*b3;
      const double sin2g = deep_arg.sing*2.*deep_arg.cosg;
      double d1dt, d2dt, d3dt, d4dt, d5dt, temp;
      double d14, d15, d16, d17, d18, d19, d23, d25, aldtal, psdtps;
      double c4dt, c5dt, c0dtc0, c1dtc1, rr2;
      double d1ddt, etddt, tmnddt, tsdtts, tsddts, xnddt, xntrdt;

      simple_flag = 0;
      edot = -c0*(eta*(eta2+4.+deep_arg.eosq*(eta2*7.+15.5))+tle->eo*
                 (eta2*15.+5.)+d1*d6*b1+c8*cos2g+c9*deep_arg.sing);
      tsdtts = deep_arg.aodp*2.*tsi*(d20*deep_arg.betao2+tle->eo*edot);
      aldtal = tle->eo*edot/alpha2;
      etdt = (edot+tle->eo*tsdtts)*tsi*s_const;
      psdtps = -eta*etdt*psim2;
      c0dtc0 = d20+tsdtts*4.-aldtal-psdtps*7.;
      c1dtc1 = xndtn+aldtal*4.+c0dtc0;
      d14 = tsdtts-psdtps*2.;
      d15 = (d20+tle->eo*edot/deep_arg.betao2)*2.;
      d1dt = d1*(d14+d15);
      d2dt = etdt*d11;
      d3dt = etdt*d12;
      d4dt = etdt*d13;
      d5dt = d5*d14;
      c4dt = b2*(d1dt*d3+d1*d3dt);
      c5dt = b3*(d5dt*d4+d5*d4dt);
      d16 = d9*etdt+d10*edot+b1*(d1dt*d2+d1*d2dt)+c4dt*
            cos2g+c5dt*deep_arg.sing+xgdt1*(c5*deep_arg.cosg-c4*2.*sin2g);
      xnddt = c1dtc1*xndt+c1*d16;
      eddot = c0dtc0*edot-c0*((eta2*3.+4.+eeta*30.+deep_arg.eosq*
        (eta2*21.+15.5))*etdt+(eta2*15.+5.+eeta*
         (eta2*14.+31.))*edot+b1*(d1dt*d6+d1*etdt*
             (eta2*67.5+30.))+b2*(d1dt*d7+d1*etdt*
        (eta2*37.5+5.))*cos2g+b3*(d5dt*d8+d5*etdt*eta*
        (eta2*4.+13.5))*deep_arg.sing+xgdt1*(c9*deep_arg.cosg-c8*2.*sin2g));
      r1 = edot;
      d25 = r1*r1;
      r1 = xndtn;
      d17 = xnddt/xnodp_-r1*r1;
      tsddts = tsdtts*2.*(tsdtts-d20)+deep_arg.aodp*tsi*
               (two_thirds*deep_arg.betao2*d17-d20*4.*tle->eo*edot+
               (d25+tle->eo*eddot)*2.);
      etddt = (eddot+edot*2.*tsdtts)*tsi*s_const+tsddts*eta;
      r1 = tsdtts;
      d18 = tsddts-r1*r1;
      r1 = psdtps;
      rr2 = psdtps;
      d19 = -(r1*r1)/eta2-eta*etddt*psim2-rr2*rr2;
      d23 = etdt*etdt;
      d1ddt = d1dt*(d14+d15)+d1*(d18-d19*2.+two_thirds*d17+
         (alpha2*d25/deep_arg.betao2+tle->eo*eddot)*2./deep_arg.betao2);
      r1  = aldtal;
      xntrdt = xndt*(two_thirds*2.*d17+(d25+tle->eo*eddot)*3./
               alpha2-r1*r1*6.+d18*4.-d19*7.)+
                    c1dtc1*xnddt+c1*(c1dtc1*d16+d9*etddt+d10*
          eddot+d23*(eeta*30.+6.+deep_arg.eosq*68.)+etdt*edot*
          (eta2*30.+40.+eeta*272.)+d25*(eta2*68.+17.)+
          b1*(d1ddt*d2+d1dt*2.*d2dt+d1*(etddt*d11+d23*
          (eta2*54.+72.)))+b2*(d1ddt*d3+d1dt*2.*d3dt+d1
          *(etddt*d12+d23*(eta2*30.+30.)))*cos2g+b3*
          ((d5dt*d14+d5*(d18-d19*2.))*d4+d4dt*2.*d5dt+
          d5*(etddt*d13+eta*22.5*d23))*deep_arg.sing+xgdt1*((d20*
          7.+tle->eo*4.*edot/deep_arg.betao2)*(c5*deep_arg.cosg-c4*2.*
          sin2g)+(c5dt*2.*deep_arg.cosg-c4dt*4.*sin2g-xgdt1*
          (c5*deep_arg.sing+c4*4.*cos2g))));
      tmnddt = xnddt*1e9;
      r1 = tmnddt;
      temp = r1*r1-xndt*1e18*xntrdt;
      r1 = tmnddt;
      pp = (temp+r1*r1)/temp;
      gamma = -xntrdt/(xnddt*(pp-2.));
      xnd = xndt/(pp*gamma);
      qq = 1.-eddot/(edot*gamma);
      ed = edot/(qq*gamma);
      ovgpp = 1./(gamma*(pp+1.));
      }
   else
      {
      simple_flag = 1;
      edot = -two_thirds*xndtn*(1.-tle->eo);
      }  /* End of if (fabs(r1) > .00216) */
} /* End of SGP8() initialization */

int DLL_FUNC SGP8( const double tsince, const tle_t *tle, const double *params,
                                           double *pos, double *vel)
{
   int i;
   double
        am, aovr, axnm, aynm, beta, beta2m,
        cose, cosos, cs2f2g, csf, csfg,
        cslamb, di, diwc, dr, ecosf, em, fm,
        g1, g10, g13, g14, g2, g3, g4, g5,
        omgasm, pm, r1, rdot, rm, rr, rvdot,
        sine, sinos, sn2f2g, snf, snfg,
        sni2du, snlamb, temp, ux, uy,
        uz, vx, vy, vz, xlamb, xmam, xn,
        xnodes, y4, y5, z1, z7, zc2, zc5;

  /* Update for secular gravity and atmospheric drag */
  r1 = tle->xmo+xlldot*tsince;
  xmam = FMod2p(r1);
  omgasm = tle->omegao+omgdt*tsince;
  xnodes = tle->xnodeo+xnodot_*tsince;
  if( !simple_flag)
    {
      double temp1;

      temp = 1.-gamma*tsince;
      temp1 = pow(temp, pp);
      xn = xnodp_+xnd*(1.-temp1);
      em = tle->eo+ed*(1.-pow(temp, qq));
      z1 = xnd*(tsince+ovgpp*(temp*temp1-1.));
    }
  else
    {
      xn = xnodp_+xndt*tsince;
      em = tle->eo+edot*tsince;
      z1 = xndt*.5*tsince*tsince;
    }  /* if(isFlagClear(SIMPLE_FLAG)) */

  z7 = two_thirds*3.5*z1/xnodp_;
  r1 = xmam+z1+z7*xmdt1;
  xmam = FMod2p(r1);
  omgasm += z7*xgdt1;
  xnodes += z7*xhdt1;

  /* Solve Kepler's equation */
  zc2 = xmam+em*sin(xmam)*(em*cos(xmam)+1.);

  i = 0;
  do
    {
      double cape;

      sine = sin(zc2);
      cose = cos(zc2);
      zc5 = 1./(1.-em*cose);
      cape = (xmam+em*sine-zc2)*zc5+zc2;
      r1 = cape-zc2;
      if(fabs(r1) <= e6a) break;
      zc2 = cape;
    }
  while(i++ < 10 );

  /* Short period preliminary quantities */
  am = pow( xke / xn, two_thirds);
  beta2m = 1.-em*em;
  sinos = sin(omgasm);
  cosos = cos(omgasm);
  axnm = em*cosos;
  aynm = em*sinos;
  pm = am*beta2m;
  g1 = 1./pm;
  g2 = ck2*.5*g1;
  g3 = g2*g1;
  beta = sqrt(beta2m);
  g4 = a3ovk2 * .25 * sini;
  g5 = a3ovk2 * .25 * g1;
  snf = beta*sine*zc5;
  csf = (cose-em)*zc5;
  fm = atan2(snf, csf);
  if( fm < 0.)
     fm += pi + pi;
  snfg = snf*cosos+csf*sinos;
  csfg = csf*cosos-snf*sinos;
  sn2f2g = snfg*2.*csfg;
  r1 = csfg;
  cs2f2g = r1*r1*2.-1.;
  ecosf = em*csf;
  g10 = fm-xmam+em*snf;
  rm = pm/(ecosf+1.);
  aovr = am/rm;
  g13 = xn*aovr;
  g14 = -g13*aovr;
  dr = g2*(unmth2*cs2f2g-tthmun*3.)-g4*snfg;
  diwc = g3*3.*sini*cs2f2g-g5*aynm;
  di = diwc*cosi;

  /* Update for short period periodics */
  sni2du = sini2*(g3*((1.-cosio_2*7.)*.5*sn2f2g-unm5th*3.*g10)-
      g5*sini*csfg*(ecosf+2.))-g5*.5f*cosio_2*axnm/cosi2;
  xlamb = fm+omgasm+xnodes+g3*((cosi*6.+1.-cosio_2*7.)*
     .5*sn2f2g-(unm5th+cosi*2.)*3.*g10)+g5*sini*
          (cosi*axnm/(cosi+1.)-(ecosf+2.)*csfg);
  y4 = sini2*snfg+csfg*sni2du+snfg*.5*cosi2*di;
  y5 = sini2*csfg-snfg*sni2du+csfg*.5*cosi2*di;
  rr = rm+dr;
  rdot = xn*am*em*snf/beta+g14*(g2*2.*unmth2*sn2f2g+g4*csfg);
  r1 = am;
  rvdot = xn*(r1*r1)*beta/rm+g14 *
          dr+am*g13*sini*diwc;

  /* Orientation vectors */
  snlamb = sin(xlamb);
  cslamb = cos(xlamb);
  temp = (y5*snlamb-y4*cslamb)*2.;
  ux = y4*temp+cslamb;
  vx = y5*temp-snlamb;
  temp = (y5*cslamb+y4*snlamb)*2.;
  uy = -y4*temp+snlamb;
  vy = -y5*temp+cslamb;
  temp = sqrt(1.-y4*y4-y5*y5)*2.;
  uz = y4*temp;
  vz = y5*temp;

  /* Position and velocity */
  pos[0] = rr*ux*earth_radius_in_km;
  pos[1] = rr*uy*earth_radius_in_km;
  pos[2] = rr*uz*earth_radius_in_km;
  if( vel)
     {
     vel[0] = (rdot*ux+rvdot*vx)*earth_radius_in_km;
     vel[1] = (rdot*uy+rvdot*vy)*earth_radius_in_km;
     vel[2] = (rdot*uz+rvdot*vz)*earth_radius_in_km;
     }
   return( 0);
} /* SGP8 */
