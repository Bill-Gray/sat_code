/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
 *  main.c          April 10  2001
 *
 *  A skeleton main() function to demonstrate the use of
 *  the various NORAD ephimerides in the norad.c file.
 *  The TLE sets used are those provided in NORAD's spacetrack
 *  report #3 so that a comparison with their own results
 *  can be made. The results produced by these software agree
 *  with NORAD's to the 5th or 6th figure, whatever this means!
 *  (But please note that NORAD uses mostly 'float' types with
 *  only a few 'doubles' for the critical variables).
 */

#include <stdio.h>
#include <string.h>
#include "norad.h"

#define INTENTIONALLY_UNUSED_PARAMETER( param) (void)(param)

static double test_data[5 * 6 * 5] = {

                        /* SGP: */
           2328.96594238,  -5995.21600342,   1719.97894287,
                       2.91110113,     -0.98164053,      -7.09049922,
           2456.00610352,  -6071.94232177,   1222.95977784,
                       2.67852119,     -0.44705850,      -7.22800565,
           2567.39477539,  -6112.49725342,    713.97710419,
                       2.43952477,      0.09884824,      -7.31899641,
           2663.03179932,  -6115.37414551,    195.73919105,
                       2.19531813,      0.65333930,      -7.36169147,
           2742.85470581,  -6079.13580322,   -328.86091614,
                       1.94707947,      1.21346101,      -7.35499924,

                        /* SGP4: */
           2328.97048951,  -5995.22076416,   1719.97067261,
                       2.91207230,     -0.98341546,      -7.09081703,
           2456.10705566,  -6071.93853760,   1222.89727783,
                       2.67938992,     -0.44829041,      -7.22879231,
           2567.56195068,  -6112.50384522,    713.96397400,
                       2.44024599,      0.09810869,      -7.31995916,
           2663.09078980,  -6115.48229980,    196.39640427,
                       2.19611958,      0.65241995,      -7.36282432,
           2742.55133057,  -6079.67144775,   -326.38095856,
                       1.94850229,      1.21106251,      -7.35619372,

                        /* SGP8: */
           2328.87265015,  -5995.21289063,   1720.04884338,
                       2.91210661,     -0.98353850,      -7.09081554,
           2456.04577637,  -6071.90490722,   1222.84086609,
                       2.67936245,     -0.44820847,      -7.22888553,
           2567.68383789,  -6112.40881348,    713.29282379,
                       2.43992555,      0.09893919,      -7.32018769,
           2663.49508667,  -6115.18182373,    194.62816810,
                       2.19525236,      0.65453661,      -7.36308974,
           2743.29238892,  -6078.90783691,   -329.73434067,
                       1.94680957,      1.21500109,      -7.35625595,

                        /* SDP4: */
           7473.37066650,    428.95261765,   5828.74786377,
                       5.10715413,      6.44468284,      -0.18613096,
          -3305.22537232,  32410.86328125, -24697.1767581,
                      -1.30113538,     -1.15131518,      -0.28333528,
          14271.28759766,  24110.46411133,  -4725.76837158,
                      -0.32050445,      2.67984074,      -2.08405289,
          -9990.05883789,  22717.35522461, -23616.89062501,
                      -1.01667246,     -2.29026759,       0.72892364,
           9787.86975097,  33753.34667969, -15030.81176758,
                      -1.09425066,      0.92358845,      -1.52230928,

                        /* SDP8: (gotta fix) */
           7469.47631836,    415.99390792,   5829.64318848,
                       5.11402285,      6.44403201,      -0.18296110,
          -3337.38992310,  32351.39086914, -24658.63037109,
                      -1.30200730,     -1.15603013,      -0.28164955,
          14226.54333496,  24236.08740234,  -4856.19744873,
                      -0.33951668,      2.65315416,      -2.08114153,
         -10151.59838867,  22223.69848633, -23392.39770508,
                      -1.00112480,     -2.33532837,       0.76987664,
           9420.08203125,  33847.21875000, -15391.06469727,
                      -1.11986055,      0.85410149,      -1.49506933 };

/* Main program */
int main( const int argc, const char **unused_argv)
{
  double vel[3], pos[3]; /* Satellite position and velocity vectors */
  double *test_ptr = test_data;

  tle_t tle; /* Pointer to two-line elements set for satellite */

  /* Data for the prediction type and time period */
  double   ts = 0.;    /* Time since TLE epoch to start predictions */
  double   tf = 1440.; /* Time over which predictions are required  */
  double delt = 360.;  /* Time interval between predictions         */

  double tsince; /* Time since epoch (in minutes) */

  int i; /* Index for loops etc */
  const char *tle_data[16] = {
      "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87",
      "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058",
      "1 11801U          80230.29629788  .01431103  00000-0  14311-1       2",
      "2 11801U 46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     2",
          /* GOES 9 */
      "1 23581U 95025A   01311.43599209 -.00000094  00000-0  00000+0 0  8214",
      "2 23581   1.1236  93.7945 0005741 214.4722 151.5103  1.00270260 23672",
          /* Cosmos 1191 */
      "1 11871U 80057A   01309.36911127 -.00000499 +00000-0 +10000-3 0 08380",
      "2 11871 067.5731 001.8936 6344778 181.9632 173.2224 02.00993562062886",
          /* ESA-GEOS 1 */
      "1 09931U 77029A   01309.17453186 -.00000329 +00000-0 +10000-3 0 05967",
      "2 09931 026.4846 264.1300 6609654 082.2734 342.9061 01.96179522175451",
          /* Cosmos 1217 */
      "1 12032U 80085A   01309.42683181  .00000182  00000-0  10000-3 0  3499",
      "2 12032  65.2329  86.7607 7086222 172.0967 212.4632  2.00879501101699",
          /* Molniya 3-19Rk */
      "1 13446U 82083E   01283.10818257  .00098407  45745-7  54864-3 0  6240",
      "2 13446  62.1717  83.8458 7498877 273.9677 320.2568  2.06357523137203",
          /* Ariane Deb */
      "1 23246U 91015G   01311.70347086  .00004957  00000-0  43218-2 0  8190",
      "2 23246   7.1648 263.6949 5661268 241.8299  50.5793  4.44333001129208" };

   INTENTIONALLY_UNUSED_PARAMETER( unused_argv);
   for( i = 1; i <= 17; i++)  /* Loop for each type of ephemeris */
      {
      int tle_idx = ((i - 2) / 2) * 2, err_code;
      int ephem, is_deep;
      const char *ephem_names[6] = { NULL, "SGP ", "SGP4", "SGP8", "SDP4", "SDP8" };
      double sat_params[N_SAT_PARAMS];

      /* Select the sgp or sdp TLE set for use below */
      if( tle_idx < 0) tle_idx = 0;
      printf( "\n%s\n%s", tle_data[tle_idx], tle_data[tle_idx + 1]);
      err_code = parse_elements( tle_data[tle_idx], tle_data[tle_idx + 1], &tle);
      if( err_code)
         printf( "\nError parsing elements: %d\n", err_code);

      if( i <= 5)
         ephem = i;
      else
         ephem = 4 + ((i - 5) % 2);

      /* Select ephemeris type */
      /* Will select a "deep" (SDPx) or "general" (SGPx) ephemeris  */
      /* depending on the TLE parameters of the satellite:          */
      is_deep = select_ephemeris( &tle);

/*    printf( "BStar: %.8lf\n", tle.bstar);  */
      if( is_deep)
         printf("\nDeep-Space type Ephemeris (SDP*) selected:");
      else
         printf("\nNear-Earth type Ephemeris (SGP*) selected:");

      /* Print some titles for the results */
      printf("\nEphem:%s   Tsince         "
             "X/Xdot           Y/Ydot           Z/Zdot\n", ephem_names[ephem]);

      /* Calling of NORAD routines */
      /* Each NORAD routine (SGP, SGP4, SGP8, SDP4, SDP8)   */
      /* will be called in turn with the appropriate TLE set */
      switch( ephem)
         {
         case 1:
            SGP_init( sat_params, &tle);
            break;
         case 2:
            SGP4_init( sat_params, &tle);
            break;
         case 3:
            SGP8_init( sat_params, &tle);
            break;
         case 4:
            SDP4_init( sat_params, &tle);
            break;
         case 5:
            SDP8_init( sat_params, &tle);
            break;
         }

      for( tsince = ts; tsince <= tf; tsince += delt)
         {
         switch( ephem)
            {
            case 1:
               SGP(tsince, &tle, sat_params, pos, vel);
               break;
            case 2:
               SGP4(tsince, &tle, sat_params, pos, vel);
               break;
            case 3:
               SGP8(tsince, &tle, sat_params, pos, vel);
               break;
            case 4:
               SDP4(tsince, &tle, sat_params, pos, vel);
               break;
            case 5:
               SDP8(tsince, &tle, sat_params, pos, vel);
               break;
            }

         /* Calculate and print results */
         vel[0] /= 60.;    /* cvt km/minute to km/second */
         vel[1] /= 60.;
         vel[2] /= 60.;

         if( argc > 1 && i <= 5)   /* wanna show _difference from test data_ */
            {
            pos[0] -= *test_ptr++;
            pos[1] -= *test_ptr++;
            pos[2] -= *test_ptr++;
            vel[0] -= *test_ptr++;
            vel[1] -= *test_ptr++;
            vel[2] -= *test_ptr++;
            }
         printf("       %12.4f   %16.8f %16.8f %16.8f \n",
                            tsince,pos[0],pos[1],pos[2]);
         printf("                      %16.8f %16.8f %16.8f \n",
                                   vel[0],vel[1],vel[2]);

         } /* End of for(tsince = ts; tsince <= tf; tsince += delt) */
      } /* End of for (i=1; i<=17; i++) */

  return(0);
} /* End of main() */

/*------------------------------------------------------------------*/
