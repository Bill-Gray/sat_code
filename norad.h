/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

/*
 *  norad.h v. 01.beta 03/17/2001
 *
 *  Header file for norad.c
 */

#ifndef NORAD_H
#define NORAD_H 1

/* #define RETAIN_PERTURBATION_VALUES_AT_EPOCH 1 */

/* Two-line-element satellite orbital data */
typedef struct
{
  double epoch, xndt2o, xndd6o, bstar;
  double xincl, xnodeo, eo, omegao, xmo, xno;
  int norad_number, bulletin_number, revolution_number;
  char classification;    /* "U" = unclassified;  only type I've seen */
  char ephemeris_type;
  char intl_desig[9];
} tle_t;

   /* NOTE: xndt2o and xndt6o are used only in the "classic" SGP, */
   /* not in SxP4 or SxP8. */
   /* epoch is a Julian Day,  UTC */
   /* xmo = mean anomaly at epoch,  radians */
   /* xno = mean motion at epoch,  radians/minute*/

#ifdef RETAIN_PERTURBATION_VALUES_AT_EPOCH
   #define DEEP_ARG_T_PARAMS     87
#else
   #define DEEP_ARG_T_PARAMS     81
#endif

#define N_SGP_PARAMS          11
#define N_SGP4_PARAMS         30
#define N_SGP8_PARAMS         25
#define N_SDP4_PARAMS        (10 + DEEP_ARG_T_PARAMS)
#define N_SDP8_PARAMS        (11 + DEEP_ARG_T_PARAMS)

/* 94 = maximum possible size of the 'deep_arg_t' structure,  in 8-byte units */
/* You can use the above constants to minimize the amount of memory used,
   but if you use the following constant,  you can be assured of having
   enough memory for any of the five models: */

#define N_SAT_PARAMS         (11 + DEEP_ARG_T_PARAMS)

/* Byte 63 of the first line of a TLE contains the ephemeris type.  The */
/* following five values are recommended,  but it seems the non-zero    */
/* values are only used internally;  "published" TLEs all have type 0.  */
/* However,  I've had occasion to produce SGP4 TLEs for high-fliers, in */
/* cases where I couldn't get an SDP4 fit.                              */

#define TLE_EPHEMERIS_TYPE_DEFAULT           0
#define TLE_EPHEMERIS_TYPE_SGP               1
#define TLE_EPHEMERIS_TYPE_SGP4              2
#define TLE_EPHEMERIS_TYPE_SDP4              3
#define TLE_EPHEMERIS_TYPE_SGP8              4
#define TLE_EPHEMERIS_TYPE_SDP8              5

#define SXPX_DPSEC_INTEGRATION_ORDER         0
#define SXPX_DUNDEE_COMPLIANCE               1
#define SXPX_ZERO_PERTURBATIONS_AT_EPOCH     2


/* SDP4 and SGP4 can return zero,  or any of the following error/warning codes.
The 'warnings' result in a mathematically reasonable value being returned,
and perigee within the earth is completely reasonable for an object that's
just left the earth or is about to hit it.  The 'errors' mean that no
reasonable position/velocity was determined.       */

#define SXPX_ERR_NEARLY_PARABOLIC         -1
#define SXPX_ERR_NEGATIVE_MAJOR_AXIS      -2
#define SXPX_WARN_ORBIT_WITHIN_EARTH      -3
#define SXPX_WARN_PERIGEE_WITHIN_EARTH    -4
#define SXPX_ERR_NEGATIVE_XN              -5
#define SXPX_ERR_CONVERGENCE_FAIL         -6

/* Function prototypes */
/* norad.c */

         /* The Win32 version can be compiled to make a .DLL,  if the     */
         /* functions are declared to be of type __stdcall... _and_ the   */
         /* functions must be declared to be extern "C",  something I     */
         /* overlooked and added 24 Sep 2002.  The DLL_FUNC macro lets    */
         /* this coexist peacefully with other OSes.                      */

#ifdef _WIN32
#define DLL_FUNC __stdcall
#else
#define DLL_FUNC
#endif

#ifdef __cplusplus
extern "C" {
#endif


void DLL_FUNC SGP_init( double *params, const tle_t *tle);
int  DLL_FUNC SGP(  const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel);

void DLL_FUNC SGP4_init( double *params, const tle_t *tle);
int  DLL_FUNC SGP4( const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel);

void DLL_FUNC SGP8_init( double *params, const tle_t *tle);
int  DLL_FUNC SGP8( const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel);

void DLL_FUNC SDP4_init( double *params, const tle_t *tle);
int  DLL_FUNC SDP4( const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel);

void DLL_FUNC SDP8_init( double *params, const tle_t *tle);
int  DLL_FUNC SDP8( const double tsince, const tle_t *tle, const double *params,
                                     double *pos, double *vel);

int DLL_FUNC select_ephemeris( const tle_t *tle);
int DLL_FUNC parse_elements( const char *line1, const char *line2, tle_t *sat);
int DLL_FUNC tle_checksum( const char *buff);
void DLL_FUNC write_elements_in_tle_format( char *buff, const tle_t *tle);

void DLL_FUNC sxpx_set_implementation_param( const int param_index,
                                              const int new_param);
void DLL_FUNC sxpx_set_dpsec_integration_step( const double new_step_size);
void DLL_FUNC lunar_solar_position( const double jd,
                    double *lunar_xyzr, double *solar_xyzr);

#ifdef __cplusplus
}                       /* end of 'extern "C"' section */
#endif

         /* Following are in 'dynamic.cpp',  for C/C++ programs that want  */
         /* to load 'sat_code.dll' and use its functions at runtime.  They */
         /* only make sense in the Win32 world: */
#ifdef _WIN32
int SXPX_init( double *params, const tle_t *tle, const int sxpx_num);
int SXPX( const double tsince, const tle_t *tle, const double *params,
                               double *pos, double *vel, const int sxpx_num);
#endif
#endif
