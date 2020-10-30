CC=g++
#CFLAGS=-I /share/apps/boost_1_49_0  
#CFLAGS=-I /apps/gcc-4.8.2/openmpi-1.8.3/boost-mpi/1.56.0/include
CFLAGS1= -I plugins/smhasher-master/src/
CFLAGS= -O3 -std=c++0x
DEBUG= -g -Wall
MPFLAG= -fopenmp

all: utilities/scaffoldExtend  utilities/combineSrAssemblyLr \
     utilities/extension2Cluster utilities/kCluster \
	 utilities/annotate  utilities/srConfirm \
	 utilities/srQuantify utilities/uniq_spliceTx \
	 utilities/updatedExtendGpd utilities/getConsensus \
	 utilities/lrCount utilities/srAbundance \
	 utilities/lrTag utilities/pathtest \
     bin/idpdenovo_core \
	 bin/pre-test 

utilities/lrTag : src/lrTag.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/scaffoldExtend : src/scaffoldExtend.cpp
	$(CC) $(CFLAGS) $(MPFLAG) -o $@  $^

utilities/combineSrAssemblyLr : src/combineSrAssemblyLr.cpp
	$(CC) $(CFLAGS) $(MPFLAG)  -o $@  $^

utilities/extension2Cluster : src/extension2Cluster.cpp
	$(CC) $(CFLAGS) $(MPFLAG) $(CFLAGS1) -o $@  $^

utilities/kCluster : src/kCluster.cpp
	$(CC) $(CFLAGS) $(MPFLAG) $(CFLAGS1)  -o $@  $^ -lrt

utilities/annotate : src/annotate.cpp
	$(CC) $(CFLAGS) $(MPFLAG) -o $@  $^ -lrt

utilities/srConfirm : src/srConfirm.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/uniq_spliceTx : src/uniq_spliceTx.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/updatedExtendGpd : src/updatedExtendGpd.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/getConsensus : src/getConsensus.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/lrCount : src/lrCount.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/srQuantify : src/srQuantify.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/srAbundance : src/srAbundance.cpp
	$(CC) $(CFLAGS) -o $@  $^

utilities/pathtest : src/pathtest.cpp
	$(CC) $(CFLAGS) -o $@  $^

bin/idpdenovo_core : src/idpdenovo.cpp
	$(CC) $(CFLAGS) -o $@  $^

bin/pre-test : src/pre-test.cpp
	$(CC) $(CFLAGS) -o $@  $^

clean :
	rm -f utilities/*
	rm -f bin/pre-test
	rm -f bin/idpdenovo_core
