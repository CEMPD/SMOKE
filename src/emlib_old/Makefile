# Version @(#)$Id$ $Source$ $Date$ 

lib:
	f77 $(FFLAGS) -c *.F *.f
	ar rv $(EDSS_LIB)/libemlib.a *.o
	rm *.o

debug:
	f77 $(DBGFLAGS) -c *.F *.f
	ar rv $(EDSS_LIB)/libemlib.debug.a *.o
	rm *.o


