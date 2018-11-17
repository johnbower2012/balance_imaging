## output file name
MAIN=parameters

## suffixes used:
#.SUFFIXES: .cpp .o .x .h

## software directories
#### include: *.h  obj: *.o  lib: *.a  src: *.cpp
SOFT=../software
IDIR=$(SOFT)/include
ODIR=$(SOFT)/build
LDIR=$(SOFT)/lib
SDIR=$(SOFT)/src
EDIR=$(SOFT)/error

## compiler and compiler flags
CC=g++
CFLAGS=-g
AFLAGS=-I$(IDIR) -fast -W -Wall -WShadow -Wconversion

## additional libraries to link
LIBS=-larmadillo

## dependecies
#### pattern substition adds folder location to *.h files
_DEPS=coshfunc.h defs.h file.h main.h parameters.h parametermap.h pricomana.h randy.h sampling.h observables_class.h
DEPS=$(patsubst %, $(IDIR)/%, $(_DEPS))
#DEPS=$(wildcard $(IDIR)/*.h)

## object files
#### pattern substition adds folder location to *.o files
_OBJ=coshfunc.o file.o main.o observables_class.o parameters.o parametermap.o pricomana.o randy.o
OBJ=$(patsubst %, $(ODIR)/%, $(_OBJ))
