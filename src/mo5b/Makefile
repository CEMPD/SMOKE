# Version "@(#)$Id$ $Source$ $Date$ 

#
lib:
	f77 $(FFLAGSM5) -c a*.F b*.F
	f77 $(FFLAGSM5) -c c*.F
	f77 $(FFLAGSM5) -c d*.F e*.F f*.F
	f77 $(FFLAGSM5) -c g*.F h*.F
	f77 $(FFLAGSM5) -c i*.F l*.F m*.F n*.F
	f77 $(FFLAGSM5) -c o*.F
	f77 $(FFLAGSM5) -c p*.F q*.F r*.F s*.F t*.F u*.F v.F w*.F y*.F z*.F

	ar rv $(EDSS_LIB)/libmo5b.a a*.o b*.o
	ar rv $(EDSS_LIB)/libmo5b.a c*.o
	ar rv $(EDSS_LIB)/libmo5b.a d*.o e*.o f*.o
	ar rv $(EDSS_LIB)/libmo5b.a g*.o h*.o
	ar rv $(EDSS_LIB)/libmo5b.a i*.o l*.o m*.F n*.F
	ar rv $(EDSS_LIB)/libmo5b.a o*.o
	ar rv $(EDSS_LIB)/libmo5b.a p*.o q*.o r*.o s*.o t*.o u*.o \
                                          v*.o w*.o y*.o z*.o
	rm *.o
 
debug:
	f77 $(DBGFLAGSM5) -c a*.F b*.F
	f77 $(DBGFLAGSM5) -c c*.F
	f77 $(DBGFLAGSM5) -c d*.F e*.F f*.F
	f77 $(DBGFLAGSM5) -c g*.F h*.F
	f77 $(DBGFLAGSM5) -c i*.F l*.F m*.F n*.F
	f77 $(DBGFLAGSM5) -c o*.F
	f77 $(DBGFLAGSM5) -c p*.F q*.F r*.F s*.F t*.F u*.F v.F w*.F y*.F z*.F
	ar rv $(EDSS_LIB)/libmo5b.debug.a a*.o b*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a c*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a d*.o e*.o f*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a g*.o h*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a i*.o l*.o m*.F n*.F
	ar rv $(EDSS_LIB)/libmo5b.debug.a o*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a p*.o q*.o r*.o s*.o t*.o u*.o \
                                          v*.o w*.o y*.o z*.o

	rm *.o
