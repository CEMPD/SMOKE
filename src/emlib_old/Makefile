# Version @(#)$Id$ $Source$ $Date$ 

lib:
#	f77 $(FFLAGS) -c *.F *.f
	f77 $(FFLAGS) -c *.F
	ar rv $(EDSS_LIB)/libemlib.a *.o
	rm *.o

debug:
#	f77 $(DBGFLAGS) -c *.F *.f
	f77 $(DBGFLAGS) -c *.F
	ar rv $(EDSS_LIB)/libemlib.debug.a *.o
	rm *.o


