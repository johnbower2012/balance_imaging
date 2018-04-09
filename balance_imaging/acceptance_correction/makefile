include ../makefile_defs.mk
BALANCE_LOCAL = ../local
ROOT_ROOT = /Applications/root_v6.04.00/lib
RHIC_HOME = ${MADAI_HOME}/rhic

ROOT_LIBS = $(shell root-config --libs)
ROOT_INCDIR = $(shell root-config --incdir)

INCLUDE = -I ${BALANCE_LOCAL}/include -I ${MADAI_INSTALLDIR}/include/madai -I ${MADAI_GSLPATH}/include -I ${ROOT_INCDIR}

LIBRARY = -L ${BALANCE_LOCAL}/lib -L ${MADAI_GSLPATH}/lib -L ${MADAI_INSTALLDIR}/lib/madai -L ${ROOT_ROOT}

corallibs :
	make -C ${RHIC_HOME}/coral install
	
b3dlibs :
	make -C ${RHIC_HOME}/b3d2 install

ballibs : 
	make -C ../newbalance install

statlibs :
	make -C ~/git/stat/rhicstat install

balance : ballibs corallibs b3dlibs
	${MADAI_CPP} ${MADAI_CFLAGS} balmain.cc ${INCLUDE} ${LIBRARY} -lbalance -lb3d -lcoralutils -lgsl -lgslcblas -o balance

mixed : ballibs corallibs b3dlibs
	${MADAI_CPP} ${MADAI_CFLAGS} mixedmain.cc ${INCLUDE} ${LIBRARY} -lbalance -lb3d -lcoralutils -lgsl -lgslcblas -o mixed

analyze : analyze.cc corallibs
	${MADAI_CPP} ${OPT} analyze.cc -o analyze

rhicstat : rhicstatmain.cc statlibs corallibs
	${MADAI_CPP} rhicstatmain.cc ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lcoralutils -lgsl -lgslcblas -lrhicstat -o rhicstat

getparsamples : getparsamplesmain.cc statlibs corallibs
	${MADAI_CPP} getparsamplesmain.cc ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lcoralutils -lgsl -lgslcblas -lrhicstat -o getparsamples