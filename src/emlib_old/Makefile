# Version @(#)$Id$ $Source$ $Date$ 

lib:
	f77 $(FFLAGS) -c *.F
	cc -c *.c
	ar rv $(DISTDIR)/libemstuff.a *.o
	rm *.o

debug:
	f77 $(DBGFLAGS) -c *.F
	cc -g -C -c *.c
	ar rv $(DISTDIR)/libemstuff.debug.a *.o
	rm *.o


