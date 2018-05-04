# sat_code

C/C++ code for the SGP4/SDP4 satellite motion model,  and for manipulating TLEs
(Two-Line Elements). Full details at
[http://www.projectpluto.com/sat_code.htm](http://www.projectpluto.com/sat_code.htm).

The only dependency is on the [`lunar`](https://github.com/Bill-Gray/lunar) (basic
astronomical ephemeris/time functions) library.  Make and `make install` that
library before attempting to build this code.

On Linux,  run `make` to build the library and various test executables.
(You can also do this with MinGW under Windows.)  In Linux,  you
can then run `make install` to put libraries in `/usr/local/lib` and some
include files in `/usr/local/include`.  (You will probably have to make that
`sudo make install`.)  For BSD,  and probably OS/X,  run `gmake CLANG=Y`
(GNU make,  with the clang compiler),  then `sudo gmake install`.

On Windows,  run `nmake -f msvc.mak` with MSVC++.  Optionally,  add
`-BITS_32=Y` for 32-bit code.
