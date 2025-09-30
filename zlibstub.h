#define gzFile FILE *
#define gzopen fopen
#define gzgets( ifile, buff, buffsize)    fgets( buff, buffsize, ifile)
#define gzclose fclose
