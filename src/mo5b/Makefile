# Version "@(#)$Id$ $Source$ $Date$ 

#
lib:
	f77 $(FFLAGS) -c *.F
	ar rv $(DISTDIR)/libmo5b.a *.o
	rm *.o
 
debug:
	f77 $(DBGFLAGS) -c *.F
	ar rv $(DISTDIR)/libmo5b.debug.a *.o
	rm *.o
