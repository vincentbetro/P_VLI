
SRCDIRS = REFINE SMOOTH IO MAIN

all:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) ; \
	  cd .. ; \
	done

clean:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) clean ; \
	  cd .. ; \
	done
	/bin/rm -f P_VLI.MACOSX

