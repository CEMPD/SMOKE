# Version @(#)$Id$ $Source$ $Date$ 

lib:
	f77 $(FFLAGS) -c *.F
	ar rv $(DISTDIR)/libemlib.a *.o
	rm *.o

debug:
	f77 $(DBGFLAGS) -c *.F
	ar rv $(DISTDIR)/libemlib.debug.a *.o
	rm *.o


