# Version "@(#)$Id$ $Source$ $Date$ 

lib:
	f77 $(FFLAGS) -Nl300 -c *.F
	ar rv $(DISTDIR)/libmo5av1.a *.o
	rm *.o

debug:
	f77 $(DBGFLAGS) -Nl300 -c *.F
	ar rv $(DISTDIR)/libmo5av1.debug.a *.o
	rm *.o
