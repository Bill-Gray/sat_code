(Translation/additions to ftp://ftp.kiam1.rssi.ru/pub/gps/spectr-rg/nu/00_read.me .
This also applies to the files for Spektr-R in the
ftp://ftp.kiam1.rssi.ru/pub/gps/spectr-r/nu/ directory.)

Format of .nu files :


Line 1:  0.000000000000E+00   // 0 = coordinates are in J2000, geocentric, equatorial
Line 2:  1.000000000000E+00   // orbit number (always zero for Spektr-RG)
Line 3:  2.011111500000E+07   // Epoch date,  YYYYMMDD (example is 2011 11 15)
Line 4:  1.492674329226E+04   // Epoch time, HH:MM:SS (example is 01:49:26.74329226)
Line 5:  5.174588077311E+00   // x-coordinate,  in thousands of kilometers
Line 6: -2.277098498781E+00   // y-coordinate
Line 7: -3.461034587339E+00   // z-coordinate
Line 8:  4.762109598043E+00   // x-velocity,  in kilometers/second
Line 9:  4.109503042673E+00   // y-velocity
Line 10: 4.560902129453E+00   // z-velocity
Line 11: 3.000000000000E-02   // Ballistic coefficient,  m^3/seconds^2/kilograms (?)
Line 12: 1.248000000000E-05   // Solar radiation coefficient,  unitless

   For Spektr-R,  the orbit number is non-zero and gives the number
of completed orbits.  It is not defined for Spektr-RG.

   Note that the file name gives the epoch in UTC.  The epoch given
within the file is three hours ahead of that (Moscow time).  The
positions match those determined by optical astrometry to within about
20 kilometers.  It is possible that the epoch is in TT (plus three
hours);  the difference is below what I can detect.

   The above solar radiation coefficient means that the sun's gravity
is counteracted by a force 1.248e-5 times as great.  It corresponds
approximately to an area/mass ratio of 0.015 m^2/kg,  which is
quite close to what we've been getting from the optical astrometry
orbit solution.
