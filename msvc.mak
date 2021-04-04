# Makefile for MSVC
all:  dropouts.exe fix_tles.exe line2.exe mergetle.exe obs_test.exe \
   obs_tes2.exe out_comp.exe sat_eph.exe sat_id.exe  \
   test2.exe test_out.exe test_sat.exe tle2mpc.exe

COMMON_FLAGS=-nologo -W3 -EHsc -c -FD -D_CRT_SECURE_NO_WARNINGS
RM=del

!ifdef BITS_32
BITS=32
!else
BITS=64
!endif

CFLAGS=-MT -O1 -D "NDEBUG" $(COMMON_FLAGS)
LINK=link /nologo /stack:0x8800

OBJS= sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj deep.obj \
     basics.obj get_el.obj common.obj tle_out.obj

dropouts.exe: dropouts.obj
   $(LINK) dropouts.obj

fix_tles.exe: fix_tles.obj sat_code$(BITS).lib
   $(LINK)    fix_tles.obj sat_code$(BITS).lib

line2.exe: line2.obj observe.obj sat_code$(BITS).lib
   $(LINK) line2.obj observe.obj sat_code$(BITS).lib lunar$(BITS).lib

mergetle.exe: mergetle.obj
   $(LINK) mergetle.obj

obs_test.exe: obs_test.obj observe.obj sat_code$(BITS).lib
   $(LINK)    obs_test.obj observe.obj sat_code$(BITS).lib

obs_tes2.exe: obs_tes2.obj observe.obj sat_code$(BITS).lib
   $(LINK)    obs_tes2.obj observe.obj sat_code$(BITS).lib

out_comp.exe: out_comp.obj
   $(LINK)    out_comp.obj

sat_code$(BITS).lib: $(OBJS)
   del sat_code$(BITS).lib
   lib /OUT:sat_code$(BITS).lib $(OBJS)

sat_id.exe: sat_id.obj sat_util.obj observe.obj sat_code$(BITS).lib
   $(LINK)  sat_id.obj sat_util.obj observe.obj sat_code$(BITS).lib lunar$(BITS).lib

sat_eph.exe: sat_eph.obj sat_util.obj observe.obj sat_code$(BITS).lib
    $(LINK)  sat_eph.obj sat_util.obj observe.obj sat_code$(BITS).lib lunar$(BITS).lib

test2.exe: test2.obj sat_code$(BITS).lib
   $(LINK) test2.obj sat_code$(BITS).lib

test_out.exe: test_out.obj sat_code$(BITS).lib
   $(LINK)    test_out.obj sat_code$(BITS).lib

test_sat.exe: test_sat.obj sat_code$(BITS).lib
   $(LINK)    test_sat.obj sat_code$(BITS).lib

tle2mpc.exe: tle2mpc.obj observe.obj sat_code$(BITS).lib
   $(LINK)   tle2mpc.obj observe.obj sat_code$(BITS).lib lunar$(BITS).lib

.cpp.obj:
   cl $(CFLAGS) $<

clean:
   del *.obj
   del *.exe
   del *.idb
   del sat_code$(BITS).lib

install:
   copy norad.h ..\myincl
   copy sat_code$(BITS).lib ..\lib
