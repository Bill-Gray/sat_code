# MSVC makefile for a DLL version
# NOTE: hasn't been used or updated in years;  some work required
# before it would be fit for purpose.

all:  test2.exe test_sat.exe obs_test.exe obs_tes2.exe sat_id.exe test_out.exe sm_sat.dll out_comp.exe

out_comp.exe: out_comp.cpp
   cl -W3 -Ox -nologo out_comp.cpp

test2.exe: test2.obj sat_code.lib
   link test2.obj sat_code.lib

test_sat.exe: test_sat.obj sat_code.lib
   cl -nologo test_sat.obj sat_code.lib

test_out.exe: test_out.obj tle_out.obj sat_code.lib
   cl -nologo test_out.obj tle_out.obj sat_code.lib

obs_test.exe: obs_test.obj observe.obj sat_code.lib
   cl -nologo obs_test.obj observe.obj sat_code.lib

obs_tes2.exe: obs_tes2.obj observe.obj sat_code.lib
   cl -nologo obs_tes2.obj observe.obj sat_code.lib

sat_id.exe:   sat_id.obj sat_util.obj sat_code.lib
   cl -nologo sat_id.obj sat_util.obj sat_code.lib

sat_code.lib: sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj deep.obj \
   basics.obj get_el.obj observe.obj common.obj
   del sat_code.lib
   del sat_code.dll
   link /DLL /IMPLIB:sat_code.lib /DEF:sat_code.def /MAP:sat_code.map \
            sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj deep.obj \
            basics.obj get_el.obj observe.obj common.obj

sm_sat.dll: sgp4.obj basics.obj get_el.obj common.obj
   del sm_sat.lib
   del sm_sat.dll
   link /DLL /IMPLIB:sm_sat.lib /DEF:sm_sat.def /MAP:sm_sat.map \
            sgp4.obj basics.obj get_el.obj common.obj

#CFLAGS=-W3 -c -LD -Ox -DRETAIN_PERTURBATION_VALUES_AT_EPOCH
CFLAGS=-W3 -c -LD -Ox -nologo -D_CRT_SECURE_NO_WARNINGS

.cpp.obj:
   cl $(CFLAGS) $<

common.obj:

sgp.obj:

sgp4.obj:

sgp8.obj:

sdp4.obj:

sdp8.obj:

deep.obj:

basics.obj:

get_el.obj:

observe.obj:

test2.obj:

test_out.obj:

test_sat.obj:

tle_out.obj:

obs_test.obj:

obs_tes2.obj:

sat_id.obj:

clean:
   del sat_code.dll
   del *.obj
   del *.exe
   del *.idb
   del sat_code.exp
   del sat_code.map
   del sat_code$(BITS).lib
