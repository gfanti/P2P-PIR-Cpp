##  Percy++
##  Copyright 2007,2013 Ian Goldberg <iang@cs.uwaterloo.ca>
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of version 2 of the GNU General Public License as
##  published by the Free Software Foundation.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  There is a copy of the GNU General Public License in the COPYING file
##  packaged with this plugin; if you cannot find it, write to the Free
##  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
##  02111-1307  USA

#Set to `true` if the build should support SPIR
SPIR_SUPPORT=false

CXXFLAGS=-Wall -g -O2 -pedantic -I/usr/local/include/NTL -std=c++11
LDLIBS=-lntl -lgmp -pthread -lsocket++ -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lgcrypt

TESTS=findpolys_test rr_test time_findpolys time_findpolys_gf28 time_findpolys_w8 time_findpolys_w16 time_findpolys_w32 testdistserver

TARGETS=pirserver pirserver_master pirclient splitdatabase

RUNFILES=database database.* out.client out.real

CLIENT_O=percyclient.o percyparams.o recover.o rsdecoder.o percyio.o \
		FXY.o gf2e.o subset_iter.o subset.o portfolio.o
SERVER_O=percyserver.o percyparams.o datastore.o percyio.o gf2e.o distserver.o \
		threadedserver.o cmdtools.o
SRCS=$(subst .o,.cc,$(CLIENT_O) $(SERVER_O) pirclient.o pirserver.o splitdatabase.o percyio.o)
LIBS=libpercyserver.a libpercyclient.a

ifeq ($(SPIR_SUPPORT),true)
	CXXFLAGS+= -I../PolyCommit -I../PBCWrapper -DSPIR_SUPPORT
	LDLIBS+= -L../PBCWrapper -lPBC -lpbc
	CLIENT_O+= ../PolyCommit/PolyCommitCommon.o pspir_crypt.o spirclient.o
	SERVER_O+= ../PolyCommit/PolyCommitCommon.o pspir_crypt.o spirserver.o
endif

all: $(TARGETS)

libs: $(LIBS)

libpercyserver.a: $(SERVER_O)
	ar rcs $@ $^

libpercyclient.a: $(CLIENT_O)
	ar rcs $@ $^

tests: $(TESTS)

pirserver: pirserver.o libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

pirserver_master: pirserver.cc libpercyserver.a
	g++ $(CXXFLAGS) -DDIST_MASTER -o $@ $^ $(LDLIBS)

pirserver_mpi: pirserver.cc mpicomm.o mpiserver.o libpercyserver.a
	mpic++ $(CXXFLAGS) -DMPI_DIST_SERVER -o $@ $^ $(LDLIBS)

pirclient: pirclient.o libpercyclient.a
	g++ -o $@ $^ $(LDLIBS)

splitdatabase: splitdatabase.o percyio.o
	g++ -o $@ $^ $(LDLIBS)

time_findpolys: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -o $@ $^ $(LDLIBS) # w128

time_findpolys_gf28: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ -static $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_GF28 -o $@ $^ $(LDLIBS)

time_findpolys_gf24: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ -static $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_GF24 -o $@ $^ $(LDLIBS)

time_findpolys_w8: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W8 -o $@ $^ $(LDLIBS)

time_findpolys_w16: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W16 -o $@ $^ $(LDLIBS)

time_findpolys_w32: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W32 -o $@ $^ $(LDLIBS)

rr_test: rsdecoder.cc FXY.o gf2e.o
	g++ $(CXXFLAGS) -DTEST_RR -o $@ $^ $(LDLIBS)

findpolys_test: rsdecoder.cc FXY.o gf2e.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTEST_FINDPOLYS -o $@ $^ $(LDLIBS)

testdistserver: testdistserver.cc cmdtools.o
	g++ -o $@ $^ -std=c++11

clean:
	-rm -f *.o

veryclean: clean
	-rm -f $(TARGETS) $(TESTS) $(LIBS)

distclean: veryclean
	-rm -f $(RUNFILES)

depend:
	makedepend -Y -- $(CXXFLAGS) -- $(SRCS) 2>/dev/null

# DO NOT DELETE

percyclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h
percyclient.o: /usr/local/include/NTL/ZZ_pX.h percyclient.h
percyclient.o: /usr/local/include/NTL/vec_GF2E.h percyresult.h percytypes.h
percyclient.o: percyparams.h gf2e.h rsdecoder.h FXY.h portfolio.h recover.h
percyclient.o: percyclient_impl.h rsdecoder_impl.h subset.h subset_iter.h
percyparams.o: percyparams.h percytypes.h percyio.h
percyparams.o: /usr/local/include/NTL/ZZ.h
recover.o: subset.h percytypes.h subset_iter.h recover.h rsdecoder.h
recover.o: percyresult.h FXY.h portfolio.h gf2e.h
rsdecoder.o: rsdecoder.h percyresult.h percytypes.h FXY.h portfolio.h
rsdecoder.o: recover.h gf2e.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
FXY.o: FXY.h
gf2e.o: gf2e.h
subset_iter.o: subset_iter.h percytypes.h
subset.o: subset.h percytypes.h
portfolio.o: portfolio.h
percyserver.o: /usr/local/include/NTL/vec_ZZ_p.h percytypes.h percyserver.h
percyserver.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h gf2e.h
percyserver.o: percyserver_impl.h
percyparams.o: percyparams.h percytypes.h percyio.h
percyparams.o: /usr/local/include/NTL/ZZ.h
datastore.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h
datastore.o: percytypes.h percyio.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
gf2e.o: gf2e.h
distserver.o: percyio.h /usr/local/include/NTL/ZZ.h distserver.h datastore.h
distserver.o: percyparams.h percytypes.h percyserver.h
distserver.o: /usr/local/include/NTL/vec_ZZ_p.h gf2e.h
threadedserver.o: percyio.h /usr/local/include/NTL/ZZ.h threadedserver.h
threadedserver.o: datastore.h percyparams.h percytypes.h percyserver.h
threadedserver.o: /usr/local/include/NTL/vec_ZZ_p.h gf2e.h
cmdtools.o: cmdtools.h
pirclient.o: /usr/local/include/NTL/ZZ_p.h percyclient.h
pirclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h
pirclient.o: /usr/local/include/NTL/vec_GF2E.h percyresult.h percytypes.h
pirclient.o: percyparams.h gf2e.h rsdecoder.h FXY.h portfolio.h recover.h
pirclient.o: percyclient_impl.h config.h version.h
pirserver.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h
pirserver.o: percytypes.h percyserver.h /usr/local/include/NTL/vec_ZZ_p.h
pirserver.o: gf2e.h config.h version.h distserver.h threadedserver.h
pirserver.o: cmdtools.h
splitdatabase.o: percyio.h /usr/local/include/NTL/ZZ.h config.h version.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
