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

#ifdef _WIN32
#define DLL_FUNC __stdcall
#else
#define DLL_FUNC
#endif

#ifdef __cplusplus
extern "C" {
#endif

void DLL_FUNC earth_lat_alt_to_parallax( const double lat,
                    const double ht_in_meters,
                    double *rho_cos_phi, double *rho_sin_phi);
void DLL_FUNC observer_cartesian_coords( const double jd, const double lon,
              const double rho_cos_phi, const double rho_sin_phi,
              double *vect);
void DLL_FUNC get_satellite_ra_dec_delta( const double *observer_loc,
                                 const double *satellite_loc, double *ra,
                                 double *dec, double *delta);
void DLL_FUNC epoch_of_date_to_j2000( const double jd, double *ra, double *dec);
void DLL_FUNC j2000_to_epoch_of_date( const double jd, double *ra, double *dec);

#ifdef __cplusplus
}                       /* end of 'extern "C"' section */
#endif
