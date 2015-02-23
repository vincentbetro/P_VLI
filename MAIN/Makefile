
include ../Makefile.options

LIBRARY_PATH += -L../IO -L../REFINE -L../SMOOTH

LIBS += -lgeometry -lrefine -lsmooth -loctree -lio -lUtility -lcgns -lhdf5 
#LIBS += -lrefine -lgeometry -lsmooth -loctree -lio -lUtility -lcgns /usr/local/hdf5-64-intel/lib/libhdf5.a  

LDFLAGS = $(LIBRARY_PATH) $(LIBS)

VLIC_SRCS = P_HUGG_VLI.cpp \

VLISRCS = $(VLIC_SRCS)
VLIOBJECTS = $(VLIC_SRCS:.cpp=.o)

MoveC_SRCS = P_Move.cpp \

MoveSRCS = $(MoveC_SRCS)
MoveOBJECTS = $(MoveC_SRCS:.cpp=.o)

all: P_VLI P_Move

P_VLI: io refine smooth $(VLIOBJECTS)
	$(LINK) -o $@ $(VLIOBJECTS) $(LDFLAGS)
	cp $@ ../$@.$(EXE_SUFFIX)
	cp $@ ../../$@.$(EXE_SUFFIX)

P_Move: io refine smooth $(MoveOBJECTS)
	$(LINK) -o $@ $(MoveOBJECTS) $(LDFLAGS)
	cp $@ ../$@.$(EXE_SUFFIX)
	cp $@ ../../$@.$(EXE_SUFFIX)

io:
	cd ../IO; $(MAKE)

refine:
	cd ../REFINE; $(MAKE)

smooth:
	cd ../SMOOTH; $(MAKE)

clean:
	/bin/rm -f *.o
	/bin/rm -f P_VLI

include ../suffix.rules
