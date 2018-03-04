/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "norad.h"

#define PI 3.141592653589793238462643383279502884197
#define TWOPI (2. * PI)
#define MINUTES_PER_DAY 1440.
#define MINUTES_PER_DAY_SQUARED (MINUTES_PER_DAY * MINUTES_PER_DAY)
#define MINUTES_PER_DAY_CUBED (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED)
#define AE 1.0
                             /* distance units, earth radii */

/* TLEs have four angles on line 2,  given in the form DDD.DDDD.  This
can be parsed more quickly as an integer,  then cast to double and
converted to radians,  all in one step.    */

static int get_angle( const char *buff)
{
   int rval = 0;

   while( *buff == ' ')
      buff++;
   while( *buff != ' ')
      {
      if( *buff != '.')
         rval = rval * 10 + (int)( *buff - '0');
      buff++;
      }
   return( rval);
}

/* Converts the quasi scientific notation of the "Motion Dot Dot/6" or
"BSTAR" field to double.  The input will always be of the form

sdddddSe

   ....where s is blank or + or -;  ddddd is a five-digit mantissa;
S is + or - or blank;  and e is a single-digit exponent.  A decimal
point is assumed before the five-digit mantissa.  */

static double sci( const char *string)
{
   double rval = 0.;

   if( string[1] != ' ')
      {
      const int ival = atoi( string);

      if( ival)
         {
         rval = (double)ival * 1.e-5;
         if( string[7] != '0')
            {
            int exponent = string[7] - '0';

            if( string[6] == '-')
               while( exponent--)
                  rval *= .1;
            else
               while( exponent--)
                  rval *= 10.;
            }
         }
      }
   return( rval);
}


/* Does a checksum modulo 10 on the given line.  Digits = their
value, '-' = 1, all other chars = 0.  Returns 0 if ok, a negative
value if it's definitely not a TLE line,  positive if it's all OK
except the checksum.  This last was added because people sometimes
want to use TLEs without worrying about the checksum. */

int DLL_FUNC tle_checksum( const char *buff)
{
   int rval = 0;
   int count = 69;

   if( (*buff != '1' && *buff != '2') || buff[1] != ' ')
      return( -1);
   while( --count)
      {
      if( *buff > '0' && *buff <= '9')
         rval += *buff - '0';
      else if( *buff == '-')
         rval++;
      if( *buff < ' ' || *buff > 'z')           /* invalid character */
         return( -2);
      buff++;
      }
   rval -= *buff++ - '0';
   if( *buff > ' ')                 /* line unterminated */
      rval = -3;
   else
      {
      rval %= 10;
      if( rval < 0)
         rval += 10;
      }
   return( rval % 10);
}

/* The "standard" SDP4 model fails badly for very high-flying satellites.
As a non-standard extension,  I'm simply storing state vectors for them,
using the following somewhat odd scheme :

1 40391U 15007B   15091.99922241 sxxxxxxxx syyyyyyyy szzzzzzzzh  9997
2 49391 [valid range, accuracy]  saaaaaaaa sbbbbbbbb scccccccc    0 8

   Epoch,  int'l & NORAD IDs are stored in the standard manner.  The
'ephemeris type' is h (rather than the otherwise universal 0).  The
xyz position and vx, vy, vz velocity are stored as 8-digit signed
hexadecimal values,  hence a range of +/- 2^32.

  x, y, z are in meters,  and hence cover a range +/- 4.29 million km.
vx, vy, vz are in 10^-4 m/s,  range +/- 429 km/s.  I do have ideas in
mind if we have to go beyond that... but let's hope that doesn't happen. */

static double get_high_value( const char *iptr)
{
   int bytes_scanned;
   unsigned value = 0;
   double rval = 0.;

   if( *iptr == '+' || *iptr == '-')
      if( sscanf( iptr + 1, "%x%n", &value, &bytes_scanned) == 1
                     && bytes_scanned == 8)
         rval = (*iptr == '-' ? -(double)value : (double)value);
   return( rval);
}

static inline double get_eight_places( const char *ptr)
{
   return( (double)atoi( ptr) + (double)atoi(ptr + 4) * 1e-8);
}

/* Meteor 2-08                                                           */
/* 1 13113U          88245.60005115 0.00000076           63463-4 0  5998 */
/* 2 13113  82.5386 288.0994 0015973 147.1294 213.0868 13.83869004325321 */

#define J2000 2451545.5
#define J1900 (J2000 - 36525. - 1.)

/* parse_elements returns:
         0 if the elements are parsed without error;
         1 if they're OK except the first line has a checksum error;
         2 if they're OK except the second line has a checksum error;
         3 if they're OK except both lines have checksum errors;
         a negative value if the lines aren't at all parseable */

int DLL_FUNC parse_elements( const char *line1, const char *line2, tle_t *sat)
{
   int rval, checksum_problem = 0;

   if( *line1 != '1' || *line2 != '2')
      rval = -4;
   else
      {
      rval = tle_checksum( line1);
      if( rval > 0)
         {
         checksum_problem = 1;  /* there's a checksum problem,  but it's */
         rval = 0;              /* not fatal; continue processing the TLE */
         }
      }

   if( rval)
      rval -= 100;
   else
      {
      rval = tle_checksum( line2);
      if( rval > 0)
         {
         checksum_problem |= 2;  /* there's a checksum problem,  but it's */
         rval = 0;               /* not fatal; continue processing the TLE */
         }
      }

   if( !rval)
      {
      char tbuff[13];
      int year = line1[19] - '0';

      if( line1[18] >= '0')
         year += (line1[18] - '0') * 10;
      if( year < 57)          /* cycle around Y2K */
         year += 100;
      sat->epoch = get_eight_places( line1 + 20) + J1900
             + (double)( year * 365 + (year - 1) / 4);
      sat->norad_number = atoi( line1 + 2);
      memcpy( tbuff, line1 + 64, 4);
      tbuff[4] = '\0';
      sat->bulletin_number = atoi( tbuff);
      sat->classification = line1[7];       /* almost always 'U' */
      memcpy( sat->intl_desig, line1 + 9, 8);
      sat->intl_desig[8] = '\0';
      memcpy( tbuff, line2 + 63, 5);
      tbuff[5] = '\0';
      sat->revolution_number = atoi( tbuff);
      sat->ephemeris_type = line1[62];
      if( sat->ephemeris_type == 'H')
         {
         size_t i;
         double *state_vect = &sat->xincl;

         for( i = 0; i < 3; i++)
            {
            state_vect[i]     = get_high_value( line1 + 33 + i * 10);
            state_vect[i + 3] = get_high_value( line2 + 33 + i * 10) * 1e-4;
            }
         return( 0);
         }

      sat->xmo = (double)get_angle( line2 + 43) * (PI / 180e+4);
      sat->xnodeo = (double)get_angle( line2 + 17) * (PI / 180e+4);
      sat->omegao = (double)get_angle( line2 + 34) * (PI / 180e+4);
      sat->xincl = (double)get_angle( line2 + 8) * (PI / 180e+4);
      sat->eo = atoi( line2 + 26) * 1.e-7;

      /* Make sure mean motion is null-terminated, since rev. no.
          may immediately follow. */
      memcpy( tbuff, line2 + 51, 12);
      tbuff[12] = '\0';
            /* Input mean motion,  derivative of mean motion and second  */
            /* deriv of mean motion,  are all in revolutions and days.   */
            /* Convert them here to radians and minutes:                 */
      sat->xno = get_eight_places( tbuff) * TWOPI / MINUTES_PER_DAY;
      sat->xndt2o = (double)atoi( line1 + 35)
                        * 1.e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
      if( line1[33] == '-')
         sat->xndt2o *= -1.;
      sat->xndd6o = sci( line1 + 44) * TWOPI / MINUTES_PER_DAY_CUBED;

      sat->bstar = sci( line1 + 53) * AE;
      }
   return( rval ? rval : checksum_problem);
}
