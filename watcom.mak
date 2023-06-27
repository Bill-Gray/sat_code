# Makefile for OpenWATCOM

all: test2.exe test_sat.exe obs_test.exe obs_tes2.exe sat_id.exe test_out.exe out_comp.exe

out_comp.exe: out_comp.cpp
   wcl386 -zq -W4 -Ox out_comp.cpp

test2.exe: test2.obj wsatlib.lib
   wcl386 -zq -k10000 test2.obj wsatlib.lib

test_sat.exe: test_sat.obj wsatlib.lib
   wcl386 -zq -k10000 test_sat.obj wsatlib.lib

obs_test.exe: obs_test.obj wsatlib.lib
   wcl386 -zq -k10000 obs_test.obj wsatlib.lib

obs_tes2.exe: obs_tes2.obj wsatlib.lib
   wcl386 -zq -k10000 obs_tes2.obj wsatlib.lib

WAT_LIB=../watlib

sat_id.exe: sat_id.obj sat_util.obj wsatlib.lib $(WAT_LIB)/wafuncs.lib
   wcl386 -zq -k10000 sat_id.obj sat_util.obj wsatlib.lib $(WAT_LIB)/wafuncs.lib

test_out.exe: test_out.obj wsatlib.lib
   wcl386 -zq -k10000 test_out.obj wsatlib.lib

#CFLAGS=-W4 -Ox -j -zq -DRETAIN_PERTURBATION_VALUES_AT_EPOCH
CFLAGS=-W4 -Ox -j -zq -i=..\include

wsatlib.lib: sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj deep.obj &
     basics.obj get_el.obj observe.obj common.obj tle_out.obj
   wlib -q wsatlib.lib  +sgp.obj +sgp4.obj +sgp8.obj +sdp4.obj +sdp8.obj
   wlib -q wsatlib.lib  +deep.obj +basics.obj +get_el.obj +observe.obj
   wlib -q wsatlib.lib  +common.obj +tle_out.obj

.cpp.obj:
   wcc386 $(CFLAGS) $<

.c.obj:
   wcc386 $(CFLAGS) $<

basics.obj:

deep.obj:

get_el.obj:

observe.obj:

common.obj:

obs_test.obj:

obs_tes2.obj:

sat_id.obj:

sat_util.obj:

tle_out.obj:

test_sat.obj:

test2.obj:

sgp.obj:

sgp4.obj:

sgp8.obj:

sdp4.obj:

sdp8.obj:

test_out.obj:

clean:
   -del *.obj
   -del *.exe
   -del *.o
   -del *.map
   -del *.lib
   -del *.exp
   -del *.dll
   -del *.a
