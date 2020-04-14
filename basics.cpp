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

#include <math.h>
#include "norad.h"
#include "norad_in.h"

/*------------------------------------------------------------------*/

/* FMOD2P */
double FMod2p( const double x)
{
   double rval = fmod( x, twopi);

   if( rval < 0.)
      rval += twopi;
   return( rval);
} /* fmod2p */

#define EPHEM_TYPE_DEFAULT    '0'
#define EPHEM_TYPE_SGP        '1'
#define EPHEM_TYPE_SGP4       '2'
#define EPHEM_TYPE_SDP4       '3'
#define EPHEM_TYPE_SGP8       '4'
#define EPHEM_TYPE_SDP8       '5'
#define EPHEM_TYPE_HIGH       'H'

/*------------------------------------------------------------------*/

void sxpall_common_init( const tle_t *tle, deep_arg_t *deep_arg);
                                                            /* common.c */

/* Selects the type of ephemeris to be used (SGP*-SDP*) */
int DLL_FUNC select_ephemeris( const tle_t *tle)
{
   int rval;

   if( tle->ephemeris_type == EPHEM_TYPE_HIGH)
      rval = 1;         /* force high-orbit state vector model */
   else if( tle->xno <= 0. || tle->eo > 1. || tle->eo < 0.)
      rval = -1;                 /* error in input data */
   else if( tle->ephemeris_type == EPHEM_TYPE_SGP4
         || tle->ephemeris_type == EPHEM_TYPE_SGP8)
      rval = 0;         /* specifically marked non-deep */
   else if( tle->ephemeris_type == EPHEM_TYPE_SDP4
         || tle->ephemeris_type == EPHEM_TYPE_SDP8)
      rval = 1;         /* specifically marked deep */
   else
      {
      deep_arg_t deep_arg;

      sxpall_common_init( tle, &deep_arg);
      /* Select a deep-space/near-earth ephemeris */
      /* If the orbital period is greater than 225 minutes... */
      if (twopi / deep_arg.xnodp >= 225.)
         rval = 1;      /* yes,  it should be a deep-space (SDPx) ephemeris */
      else
         rval = 0;      /* no,  you can go with an SGPx ephemeris */
      }
   return( rval);
} /* End of select_ephemeris() */

/*------------------------------------------------------------------*/

long DLL_FUNC sxpx_library_version( void)
{
   return( 0x100);
}
