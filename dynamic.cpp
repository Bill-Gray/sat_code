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

#include <stdio.h>
#include <stdlib.h>
#include "windows.h"
#include "norad.h"

typedef void (__stdcall *sxpx_init_fn)( double *params, const tle_t *tle);
typedef int (__stdcall *sxpx_fn)( const double tsince, const tle_t *tle,
                     const double *params, double *pos, double *vel);
typedef long (__stdcall *sxpx_version_fn)( void);

static HINSTANCE load_sat_code_lib( const int unload)
{
   static HINSTANCE h_sat_code_lib = (HINSTANCE)0;
   static int first_time = 1;

   if( unload)
      {
      if( h_sat_code_lib)
         FreeLibrary( h_sat_code_lib);
      h_sat_code_lib = NULL;
      first_time = 1;
      }
   else if( first_time)
      {
      h_sat_code_lib = LoadLibrary( "sat_code.dll");
      first_time = 0;
      }
   return( h_sat_code_lib);
}

/* 26 Nov 2002:  revised following two functions slightly so that the
   return values distinguish between "didn't get the function" and
   "didn't get the library" */

int SXPX_init( double *params, const tle_t *tle, const int sxpx_num)
{
   static sxpx_init_fn func[5];
   static char already_done[5];
   int rval = 0;
   HINSTANCE h_sat_code_lib;

   if( !params)         /* flag to unload library */
      {
      int i;

      load_sat_code_lib( -1);
      for( i = 0; i < 5; i++)
         already_done[i] = 0;
      return( 0);
      }
   h_sat_code_lib = load_sat_code_lib( 0);
   if( !already_done[sxpx_num])
      {
      if( h_sat_code_lib)
         func[sxpx_num] = (sxpx_init_fn)GetProcAddress( h_sat_code_lib,
                                (LPCSTR)( sxpx_num + 1));
      already_done[sxpx_num] = 1;
      }
   if( func[sxpx_num])
      (*func[sxpx_num])( params, tle);
   else
      rval = -1;
   if( !h_sat_code_lib)
      rval = -2;
   return( rval);
}

int SXPX( const double tsince, const tle_t *tle, const double *params,
                               double *pos, double *vel, const int sxpx_num)
{
   static sxpx_fn func[5];
   static char already_done[5];
   int rval = 0;
   HINSTANCE h_sat_code_lib = load_sat_code_lib( 0);

   if( !already_done[sxpx_num])
      {
      if( h_sat_code_lib)
         func[sxpx_num] = (sxpx_fn)GetProcAddress( h_sat_code_lib,
                                (LPCSTR)( sxpx_num + 6));
      already_done[sxpx_num] = 1;
      }
   if( func[sxpx_num])
      rval = (*func[sxpx_num])( tsince, tle, params, pos, vel);
   else
      rval = -1;
   if( !h_sat_code_lib)
      rval = -2;
   return( rval);
}

long get_sat_code_lib_version( void)
{
   HINSTANCE h_sat_code_lib = load_sat_code_lib( 0);
   long rval;

   if( !h_sat_code_lib)
      rval = -2;
   else
      {
      sxpx_version_fn func =
               (sxpx_version_fn)GetProcAddress( h_sat_code_lib, (LPCSTR)20);

      if( !func)
         rval = -1;
      else
         rval = (*func)( );
      }
   return( rval);
}
