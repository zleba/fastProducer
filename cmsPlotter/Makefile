###########################################################################################
#
#  simple Makefile for fastNLO for Diffraction
#
###########################################################################################

# please set your fastNLO prefix path here, or use fnlo-tk-config
#FastNLOInstallDir=$(fnlo-tk-config --prefix)
FastNLOInstallDir=../fastnloBuild/

ROOT_LIBS=$(shell root-config --libs)
ROOT_INCLUDE=$(shell root-config --cflags)

LHA_LIBS=$(shell lhapdf-config --libdir)
LHA_INCLUDE=$(shell lhapdf-config --incdir)

ARMA_DIR=../arma/install



# flags and variables
OBJS 	= fnlo-tk-h1diffpdf.o 
LDFLAGS	= -lm -lpthread  -lgfortran  -lstdc++   # `gtk-config --cflags` `gtk-config --libs` 
CC 	= g++
CFLAGS 	= -Wall -g -O2 -I$(FastNLOInstallDir)/include 
EXEC	= fastNLO



runFastNLO: runFastNLO.cc
	$(CC) -g  runFastNLO.cc $(LDFLAGS) \
	$(ROOT_INCLUDE)  -I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-I$(LHA_INCLUDE)     \
	-L$(LHA_LIBS) -lLHAPDF \
	-Wl,-rpath $(LHA_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@


checkTables: checkTables.cc
	$(CC) -g  checkTables.cc $(LDFLAGS) -I../PlottingHelper/ \
	$(ROOT_INCLUDE)  -I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-I$(LHA_INCLUDE)     \
	-L$(LHA_LIBS) -lLHAPDF \
	-Wl,-rpath $(LHA_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@


plotJets: plotJets.cc
	$(CC) -g  plotJets.cc $(LDFLAGS) -I../PlottingHelper/ \
	$(ROOT_INCLUDE)  -I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-I$(LHA_INCLUDE)     \
	-L$(LHA_LIBS) -lLHAPDF \
	-Wl,-rpath $(LHA_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@

calcTheory: calcTheory.cc
	$(CC) -g -O2  calcTheory.cc $(LDFLAGS) -I../PlottingHelper/ \
	$(ROOT_INCLUDE)  -I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-I$(LHA_INCLUDE)     \
	-L$(LHA_LIBS) -lLHAPDF \
	-Wl,-rpath $(LHA_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@



fitTheory: fitTheory.cc
	$(CC) -g -O2  $< $(LDFLAGS) -I../PlottingHelper/ \
	$(ROOT_INCLUDE)  -I$(FastNLOInstallDir)/include -L$(FastNLOInstallDir)/lib -lfastnlotoolkit \
	-L../PlottingHelper/ -lPlottingHelper -Wl,-rpath,../PlottingHelper   \
	-I$(LHA_INCLUDE)     \
	-L$(LHA_LIBS) -lLHAPDF \
	-Wl,-rpath $(LHA_LIBS) \
	$(ROOT_LIBS)   -Wl,-rpath $(FastNLOInstallDir)/lib \
	-o $@
