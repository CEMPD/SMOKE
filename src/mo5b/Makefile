# Version "@(#)$Id$ $Source$ $Date$ 

#
lib:
	f77 $(FFLAGSM5) -c *.F
	ar rv $(EDSS_LIB)/libmo5b.a *.o
	rm *.o
 
debug:
	f77 $(DBGFLAGSM5) -c *.F
	ar rv $(EDSS_LIB)/libmo5b.debug.a *.o
	rm *.o
