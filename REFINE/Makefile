
include ../Makefile.options

SRCS=\
	viscous_insert.cpp \
	mesh_deform.cpp \

OBJECTS = $(SRCS:.cpp=.o)

refine.a:	$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

refine.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f librefine.a
	/bin/rm -f librefine.so

include ../suffix.rules
