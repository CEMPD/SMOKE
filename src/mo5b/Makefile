# Version "@(#)$Id$ $Source$ $Date$ 

#
lib:
	f77 $(FFLAGSM5) -c a*.F b*.F
	f77 $(FFLAGSM5) -c c*.F
	f77 $(FFLAGSM5) -c d*.F e*.F f*.F
	f77 $(FFLAGSM5) -c g*.F h*.F
	f77 $(FFLAGSM5) -c i*.F l*.F m*.F n*.F
	f77 $(FFLAGSM5) -c ob*.F on*.F ot*.F outa*.F outb*.F outc*.F outd*.F \
                           oute*.F outg*.F 
	f77             -c outhd1.F 
	f77 $(FFLAGSM5) -c outhd2.F outhd3.F outhd4.F outhdr.F outi*.F outm*.F \
                           outn*.F outp*.F outs*.F outt*.F outz*.F ov*.F 
	f77 $(FFLAGSM5) -c p*.F q*.F r*.F s*.F t*.F u*.F v*.F w*.F y*.F z*.F

	ar rv $(EDSS_LIB)/libmo5b.a a*.o b*.o
	ar rv $(EDSS_LIB)/libmo5b.a c*.o
	ar rv $(EDSS_LIB)/libmo5b.a d*.o e*.o f*.o
	ar rv $(EDSS_LIB)/libmo5b.a g*.o h*.o
	ar rv $(EDSS_LIB)/libmo5b.a i*.o l*.o m*.o n*.o
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
	f77 $(DBGFLAGSM5) -c p*.F q*.F r*.F s*.F t*.F u*.F v*.F w*.F y*.F z*.F
	ar rv $(EDSS_LIB)/libmo5b.debug.a a*.o b*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a c*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a d*.o e*.o f*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a g*.o h*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a i*.o l*.o m*.o n*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a o*.o
	ar rv $(EDSS_LIB)/libmo5b.debug.a p*.o q*.o r*.o s*.o t*.o u*.o \
                                          v*.o w*.o y*.o z*.o

	rm *.o
