/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/* tle_out.cpp: code to create ASCII TLEs (Two-Line Elements)

/* Code to convert the in-memory artificial satellite elements into
the "standard" TLE (Two-Line Element) form described at

https://en.wikipedia.org/wiki/Two-line_elements       */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "norad.h"

      /* Useful constants to define,  in case the value of PI or the number
         of minutes in a day should change: */
#define PI 3.141592653589793238462643383279502884197169399375105
#define MINUTES_PER_DAY 1440.
#define MINUTES_PER_DAY_SQUARED (MINUTES_PER_DAY * MINUTES_PER_DAY)
#define MINUTES_PER_DAY_CUBED (MINUTES_PER_DAY_SQUARED * MINUTES_PER_DAY)

#define AE 1.0
#define J1900 (2451545.5 - 36525. - 1.)

static int add_tle_checksum_data( char *buff)
{
   int count = 69, rval = 0;

   if( (*buff != '1' && *buff != '2') || buff[1] != ' ')
      return( 0);    /* not a .TLE */
   while( --count)
      {
      if( *buff < ' ' || *buff > 'Z')
         return( 0);             /* wups!  those shouldn't happen in .TLEs */
      if( *buff > '0' && *buff <= '9')
         rval += *buff - '0';
      else if( *buff == '-')
         rval++;
      buff++;
      }
   if( *buff != 10 && buff[1] != 10 && buff[2] != 10)
      return( 0);                              /* _still_ not a .TLE */
   *buff++ = (char)( '0' + (rval % 10));
   *buff++ = 13;
   *buff++ = 10;
   *buff++ = '\0';
   return( 1);
}

static double zero_to_two_pi( double angle_in_radians)
{
   angle_in_radians = fmod( angle_in_radians, PI + PI);
   if( angle_in_radians < 0.)
      angle_in_radians += PI + PI;
   return( angle_in_radians);
}

/* See comments for get_high_value() in 'get_el.cpp'.  Essentially,  we are
   writing out a state vector in a convoluted base-36 form.  */

static void set_high_value( char *obuff, const double value)
{
   int64_t oval = (int64_t)fabs( value);
   int i;

   *obuff++ = (value >= 0. ? '+' : '-');
   for( i = 7; i >= 0; i--, oval /= (int64_t)36)
      {
      obuff[i] = (char)( oval % (int64_t)36);
      if( obuff[i] < 10)
         obuff[i] += '0';
      else
         obuff[i] += 'A' - 10;
      }
   obuff[8] = ' ';
}


/* The second derivative of the mean motion,  'xnddo6',  and the 'bstar'
drag term,  are stored in a simplified scientific notation.  To save
valuable bytes,  the decimal point and 'E' are assumed.     */

static void put_sci( char *obuff, double ival)
{
   if( !ival)
      memcpy( obuff, " 00000-0", 7);
   else
      {
      int oval, exponent = 0;

      if( ival > 0.)
         *obuff++ = ' ';
      else
         {
         *obuff++ = '-';
         ival = -ival;
         }
      while( true)
         {
         if( ival > 1.)    /* avoid integer overflow */
            oval = 100000;
         else
            oval = (int)( ival * 100000. + .5);
         if( oval > 99999)
            {
            ival /= 10;
            exponent++;
            }
         else if( oval < 10000)
            {
            ival *= 10;
            exponent--;
            }
         else
            break;
         }
      sprintf( obuff, "%5d", oval);
      if( exponent > 0)
         {
         obuff[5] = '+';
         obuff[6] = (char)( '0' + exponent);
         }
      else
         {
         obuff[5] = '-';
         obuff[6] = (char)( '0' - exponent);
         }
      }
}

/* SpaceTrack TLEs have,  on the second line,  leading zeroes in front of the
inclination,  ascending node,  argument of perigee,  and mean motion.  Which
is why I've used this format string :

      sprintf( line2 + 8, "%08.4f %08.4f %07ld %08.4f %08.4f %011.8f", ...)

   'classfd.tle' and some other sources don't use leading zeroes.  For them,
use the following format string for those four quantities :

      sprintf( line2 + 8, "%8.4f %8.4f %07ld %8.4f %8.4f %11.8f", ...)  */

void DLL_FUNC write_elements_in_tle_format( char *buff, const tle_t *tle)
{
   long year = (long)( tle->epoch - J1900) / 365 + 1;
   double day_of_year;
   char *line2;

   do
      {
      double start_of_year;

      year--;
      start_of_year = J1900 + (double)year * 365. + (double)((year - 1) / 4);
      day_of_year = tle->epoch - start_of_year;
      }
      while( day_of_year < 1.);
   if( year < 0 || year > 200)      /* bogus input */
      {
      year = 56;
      day_of_year = 0.;
      }
   sprintf( buff,
/*                                     xndt2o    xndd6o   bstar  eph bull */
           "1 %05d%c %-8s %02ld%12.8f -.000hit00 +00000-0 +00000-0 %c %4dZ\n",
           tle->norad_number, tle->classification, tle->intl_desig,
           year % 100L, day_of_year,
           tle->ephemeris_type, tle->bulletin_number);
   if( buff[20] == ' ')       /* fill in leading zeroes for day of year */
      buff[20] = '0';
   if( buff[21] == ' ')
      buff[21] = '0';
   if( tle->ephemeris_type != 'H')     /* "normal",  standard TLEs */
      {
      double deriv_mean_motion = tle->xndt2o * MINUTES_PER_DAY_SQUARED / (2. * PI);
      if( deriv_mean_motion >= 0)
         buff[33] = ' ';
      deriv_mean_motion = fabs( deriv_mean_motion * 100000000.) + .5;
      sprintf( buff + 35, "%08ld", (long)deriv_mean_motion);
      buff[43] = ' ';
      put_sci( buff + 44, tle->xndd6o * MINUTES_PER_DAY_CUBED / (2. * PI));
      put_sci( buff + 53, tle->bstar / AE);
      }
   else
      {
      size_t i;
      const double *posn = &tle->xincl;

      for( i = 0; i < 3; i++)
         set_high_value( buff + 33 + i * 10, posn[i]);
      buff[62] = 'H';
      }
   add_tle_checksum_data( buff);
   line2 = buff + strlen( buff);
   sprintf( line2, "2 %05d ", tle->norad_number);
   if( tle->ephemeris_type != 'H')     /* "normal",  standard TLEs */
      sprintf( line2 + 8, "%08.4f %08.4f %07ld %08.4f %08.4f %011.8f",
           zero_to_two_pi( tle->xincl) * 180. / PI,
           zero_to_two_pi( tle->xnodeo) * 180. / PI,
           (long)( tle->eo * 10000000. + .5),
           zero_to_two_pi( tle->omegao) * 180. / PI,
           zero_to_two_pi( tle->xmo) * 180. / PI,
           tle->xno * MINUTES_PER_DAY / (2. * PI));
   else
      {
      size_t i;
      const double *vel = &tle->xincl + 3;

      memset( line2 + 8, ' ', 25);     /* reserved for future use */
      for( i = 0; i < 3; i++)
         set_high_value( line2 + 33 + i * 10, vel[i] * 1e+4);
      }
   sprintf( line2 + 63, "%5dZ\n", tle->revolution_number);
   add_tle_checksum_data( line2);
}
