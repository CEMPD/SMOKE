# Version "@(#)$Id$ $Source$ $Date$ 

#
lib:
	f77 $(FFLAGS) -Nl300 -c *.F
	ar rv $(DISTDIR)/libmo5b.a *.o
	rm *.o
 
debug:
	f77 $(DBGFLAGS) -Nl300 -c *.F
	ar rv $(DISTDIR)/libmo5b.debug.a *.o
	rm *.o
