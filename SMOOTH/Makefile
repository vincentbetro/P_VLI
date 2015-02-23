
include ../Makefile.options

SRCS=\
	displaced.cpp \
	output_GMRES.cpp \
	output_K.cpp \
	smooth.cpp \
	mesh_move.cpp \

OBJECTS = $(SRCS:.cpp=.o)

smooth.a:	$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

smooth.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libsmooth.a
	/bin/rm -f libsmooth.so

include ../suffix.rules
