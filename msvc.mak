# Makefile for MSVC
all:  test2.exe test_sat.exe obs_test.exe obs_tes2.exe sat_id.exe out_comp.exe

!ifdef BITS_32
COMMON_FLAGS=-nologo -W3 -EHsc -c -FD
LUNAR=lunar
RM=rm
!else
COMMON_FLAGS=-nologo -W3 -EHsc -c -FD -D_CRT_SECURE_NO_WARNINGS
LUNAR=lunar64
RM=del
!endif

CFLAGS=-MT -O1 -D "NDEBUG" $(COMMON_FLAGS)
LINK=link /nologo /stack:0x8800

OBJS= sgp.obj sgp4.obj sgp8.obj sdp4.obj sdp8.obj deep.obj \
     basics.obj get_el.obj common.obj tle_out.obj

sat_code.lib: $(OBJS)
   del sat_code.lib
   lib /OUT:sat_code.lib $(OBJS)

test2.exe: test2.obj sat_code.lib
   $(LINK) test2.obj sat_code.lib

test_sat.exe: test_sat.obj sat_code.lib
   $(LINK)    test_sat.obj sat_code.lib

obs_test.exe: obs_test.obj observe.obj sat_code.lib
   $(LINK)    obs_test.obj observe.obj sat_code.lib

obs_tes2.exe: obs_tes2.obj observe.obj sat_code.lib
   $(LINK)    obs_tes2.obj observe.obj sat_code.lib

sat_id.exe: sat_id.obj observe.obj sat_code.lib
   $(LINK)  sat_id.obj observe.obj sat_code.lib $(LUNAR).lib

out_comp.exe: out_comp.obj
   $(LINK)    out_comp.obj

.cpp.obj:
   cl $(CFLAGS) $<

clean:
   del *.obj
   del *.exe
   del *.idb
   del sat_code.lib

install:
   copy norad.h ..\myincl
