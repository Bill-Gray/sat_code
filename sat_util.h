#ifndef SAT_UTIL_H_INCLUDED
#define SAT_UTIL_H_INCLUDED

/* A few functions that are used in common by sat_id, sat_id2, and sat_eph. */

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

FILE *local_then_config_fopen( const char *filename, const char *permits);
char *fgets_trimmed( char *buff, const int buffsize, FILE *ifile);

#if !defined( _WIN32)
void make_config_dir_name( char *oname, const char *iname);
#endif

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */
#endif  /* #ifndef SAT_UTIL_H_INCLUDED */
