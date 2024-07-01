# GNU MAKE Makefile for artsat code and utilities
#
# Usage: make [W64=Y] [W32=Y] [MSWIN=Y] [tgt]
#
#	where tgt can be any of:
# [all|get_high|mergetle|obs_tes2|...]
#
# ...see below for complete list.  Note that you have to 'make tle_date' and/or
# 'make tle_date.cgi' separately;  they have a dependency on the 'lunar'
# library (also on my GitHub site).
#
#	'W64'/'W32' = cross-compile for 64- or 32-bit Windows,  using MinGW,
#     on a Linux box
#	'MSWIN' = compile for Windows,  using MinGW and PDCurses,  on a Windows machine
#	'CC=clang' or 'CC=g++-4.8' = use clang or older GCC
#    (I've used gmake CLANG=Y on PC-BSD;  probably works on OS/X too)
# None of these: compile using default g++ on Linux,  for Linux
#

# As CC is an implicit variable, a simple CC?=g++ doesn't work.
# We have to use this trick from https://stackoverflow.com/a/42958970
ifeq ($(origin CC),default)
	CC=gcc
	CXX=g++
endif

ifeq ($(shell uname -s),FreeBSD)
	CC=cc
	CXX=c++
endif

EXE=
RM=rm -f

# I'm using 'mkdir -p' to avoid error messages if the directory exists.
# It may fail on very old systems,  and will probably fail on non-POSIX
# systems.  If so,  change to '-mkdir' and ignore errors.

ifdef MSWIN
	EXE=.exe
	MKDIR=-mkdir
else
	MKDIR=mkdir -p
endif

LIB_DIR=$(INSTALL_DIR)/lib

ifdef W64
	CC=x86_64-w64-mingw32-gcc
	CXX=x86_64-w64-mingw32-g++
	EXE=.exe
	LIB_DIR=$(INSTALL_DIR)/win_lib
endif

ifdef W32
	CC=i686-w64-mingw32-gcc
	CXX=i686-w64-mingw32-g++
	EXE=.exe
	LIB_DIR=$(INSTALL_DIR)/win_lib32
endif

# You can have your include files in ~/include and libraries in
# ~/lib,  in which case only the current user can use them;  or
# (with root privileges) you can install them to /usr/local/include
# and /usr/local/lib for all to enjoy.

PREFIX?=~
ifdef GLOBAL
	INSTALL_DIR=/usr/local
else
	INSTALL_DIR=$(PREFIX)
endif

INCL=$(INSTALL_DIR)/include

all: dropouts$(EXE) fake_ast$(EXE) fix_tles$(EXE) get_high$(EXE) \
	line2$(EXE) mergetle$(EXE) obs_tes2$(EXE) obs_test$(EXE) \
	out_comp$(EXE) sat_cgi$(EXE) sat_eph$(EXE) sat_id$(EXE) \
	sat_id2$(EXE) sat_id3$(EXE) summarize$(EXE) \
	test_des$(EXE) test_out$(EXE) test_sat$(EXE) test2$(EXE) tle2mpc$(EXE)

CFLAGS+=-Wextra -Wall -O3 -pedantic -Wshadow

ifdef UCHAR
	CFLAGS += -funsigned-char
endif

ifdef DEBUG
	CFLAGS += -g
endif

ifndef NO_ERRORS
	CFLAGS += -Werror
endif


clean:
	$(RM) *.o
	$(RM) dropouts$(EXE)
	$(RM) fake_ast$(EXE)
	$(RM) fix_tles$(EXE)
	$(RM) get_high$(EXE)
	$(RM) get_vect$(EXE)
	$(RM) libsatell.a
	$(RM) line2$(EXE)
	$(RM) mergetle$(EXE)
	$(RM) obs_tes2$(EXE)
	$(RM) obs_test$(EXE)
	$(RM) out_comp$(EXE)
	$(RM) sat_cgi$(EXE)
	$(RM) sat_eph$(EXE)
	$(RM) sat_id$(EXE)
	$(RM) sat_id2$(EXE)
	$(RM) sat_id3$(EXE)
	$(RM) summarize$(EXE)
	$(RM) test2$(EXE)
	$(RM) test_des$(EXE)
	$(RM) test_out$(EXE)
	$(RM) test_sat$(EXE)
	$(RM) tle2mpc$(EXE)
	$(RM) tle_date$(EXE)
	$(RM) tle_date.cgi

install:
	$(MKDIR) $(LIB_DIR)
	cp libsatell.a $(LIB_DIR)
	cp norad.h     $(INSTALL_DIR)/include
	$(MKDIR) $(INSTALL_DIR)/bin
	cp sat_id$(EXE)  $(INSTALL_DIR)/bin

install_lib:
	$(MKDIR) $(LIB_DIR)
	cp libsatell.a $(LIB_DIR)
	cp norad.h     $(INSTALL_DIR)/include

uninstall:
	rm $(INSTALL_DIR)/lib/libsatell.a
	rm $(INSTALL_DIR)/include/norad.h
	rm $(INSTALL_DIR)/bin/sat_id

uninstall_lib:
	rm $(INSTALL_DIR)/lib/libsatell.a
	rm $(INSTALL_DIR)/include/norad.h

OBJS= sgp.o sgp4.o sgp8.o sdp4.o sdp8.o deep.o basics.o get_el.o common.o tle_out.o

get_high$(EXE):	 get_high.o get_el.o
	$(CC) $(CFLAGS) -o get_high$(EXE) get_high.o get_el.o

mergetle$(EXE):	 mergetle.o
	$(CC) $(CFLAGS) -o mergetle$(EXE) mergetle.o -lm

dropouts$(EXE):	 dropouts.o
	$(CC) $(CFLAGS) -o dropouts$(EXE) dropouts.o

obs_tes2$(EXE):	 obs_tes2.o observe.o libsatell.a
	$(CC) $(CFLAGS) -o obs_tes2$(EXE) obs_tes2.o observe.o libsatell.a -lm

obs_test$(EXE):	 obs_test.o observe.o libsatell.a
	$(CC) $(CFLAGS) -o obs_test$(EXE) obs_test.o observe.o libsatell.a -lm

fake_ast$(EXE):	 fake_ast.o observe.o libsatell.a
	$(CC) $(CFLAGS) -o fake_ast$(EXE) fake_ast.o observe.o libsatell.a -lm

fix_tles$(EXE):	 fix_tles.o libsatell.a
	$(CC) $(CFLAGS) -o fix_tles$(EXE) fix_tles.o libsatell.a -lm

get_vect$(EXE):	 	get_vect.cpp	observe.o libsatell.a
	$(CXX) $(CFLAGS) -o get_vect$(EXE) -I $(INCL) get_vect.cpp observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

line2$(EXE):	 	line2.cpp libsatell.a
	$(CXX) $(CFLAGS) -o line2$(EXE) -I $(INCL) line2.cpp libsatell.a -lm -L $(LIB_DIR) -llunar

out_comp$(EXE):	 out_comp.o
	$(CC) $(CFLAGS) -o out_comp$(EXE) out_comp.o -lm

libsatell.a: $(OBJS)
	rm -f libsatell.a
	ar rv libsatell.a $(OBJS)

sat_eph$(EXE):	 	sat_eph.c	observe.o libsatell.a
	$(CC) $(CFLAGS) -o sat_eph$(EXE) -I $(INCL) sat_eph.c observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

sat_cgi$(EXE):	 	sat_eph.c	observe.o libsatell.a
	$(CC) $(CFLAGS) -o sat_cgi$(EXE) -I $(INCL) sat_eph.c observe.o -DON_LINE_VERSION libsatell.a -lm -L $(LIB_DIR) -llunar

sat_id$(EXE):	 	sat_id.cpp sat_util.o	observe.o libsatell.a
	$(CXX) $(CFLAGS) -o sat_id$(EXE) -I $(INCL) sat_id.cpp sat_util.o observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

sat_id2$(EXE):	 	sat_id2.cpp sat_id.cpp sat_util.o observe.o libsatell.a
	$(CXX) $(CFLAGS) -o sat_id2$(EXE) -I $(INCL) -DON_LINE_VERSION sat_id2.cpp sat_id.cpp sat_util.o observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

sat_id3$(EXE):	 	sat_id3.cpp sat_id.cpp sat_util.o observe.o libsatell.a
	$(CXX) $(CFLAGS) -o sat_id3$(EXE) -I $(INCL) -DON_LINE_VERSION sat_id3.cpp sat_id.cpp sat_util.o observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

summarize$(EXE):	 	summarize.c	observe.o libsatell.a
	$(CC) $(CFLAGS) -o summarize$(EXE) -I $(INCL) summarize.c observe.o libsatell.a -lm -L $(LIB_DIR) -llunar

test2$(EXE):	 	test2.o sgp.o libsatell.a
	$(CC) $(CFLAGS) -o test2$(EXE) test2.o sgp.o libsatell.a -lm

tle_date$(EXE):	 	tle_date.o
	$(CC) $(CFLAGS) -o tle_date$(EXE) tle_date.o -L $(LIB_DIR) -llunar -lm

tle_date.o: tle_date.c
	$(CC) $(CFLAGS) -o tle_date.o -c -I../include tle_date.c

tle_date.cgi:	 	tle_date.c
	$(CC) $(CFLAGS) -o tle_date.cgi  -I../include -DON_LINE_VERSION tle_date.c -L $(LIB_DIR) -llunar

tle2mpc$(EXE):	 	tle2mpc.cpp libsatell.a
	$(CXX) $(CFLAGS) -o tle2mpc$(EXE) -I $(INCL) tle2mpc.cpp libsatell.a -lm -L $(LIB_DIR) -llunar

test_des$(EXE):	 test_des.o libsatell.a
	$(CC) $(CFLAGS) -o test_des$(EXE) test_des.o libsatell.a -lm

test_out$(EXE):	 test_out.o tle_out.o get_el.o sgp4.o common.o
	$(CC) $(CFLAGS) -o test_out$(EXE) test_out.o tle_out.o get_el.o sgp4.o common.o -lm

test_sat$(EXE):	 test_sat.o libsatell.a
	$(CC) $(CFLAGS) -o test_sat$(EXE) test_sat.o libsatell.a -lm

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<
