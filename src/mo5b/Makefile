# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 

# Makefile for Mobile5b library

all: lib

lib:
	$(FC) $(FFLAGS) -c a*.F b*.F
	$(FC) $(FFLAGS) -c c*.F
	$(FC) $(FFLAGS) -c d*.F e*.F f*.F
	$(FC) $(FFLAGS) -c g*.F h*.F
	$(FC) $(FFLAGS) -c i*.F l*.F m*.F n*.F
	$(FC) $(FFLAGS) -c ob*.F on*.F ot*.F outa*.F outb*.F outc*.F outd*.F \
                           oute*.F outg*.F 
	$(FC)           -c outhd1.F 
	$(FC) $(FFLAGS) -c outhd2.F outhd3.F outhd4.F outhdr.F outi*.F \
                           outm*.F outn*.F outp*.F outs*.F outt*.F outz*.F ov*.F
	$(FC) $(FFLAGS) -c p*.F q*.F r*.F s*.F t*.F u*.F v*.F w*.F y*.F z*.F

	ar rv $(SMK_BIN)/libmo5b.a a*.o b*.o
	ar rv $(SMK_BIN)/libmo5b.a c*.o
	ar rv $(SMK_BIN)/libmo5b.a d*.o e*.o f*.o
	ar rv $(SMK_BIN)/libmo5b.a g*.o h*.o
	ar rv $(SMK_BIN)/libmo5b.a i*.o l*.o m*.o n*.o
	ar rv $(SMK_BIN)/libmo5b.a o*.o
	ar rv $(SMK_BIN)/libmo5b.a p*.o q*.o r*.o s*.o t*.o u*.o \
                                   v*.o w*.o y*.o z*.o
#	rm *.o
 
debug:
	$(FC) $(DBGFLAGS) -c a*.F b*.F
	$(FC) $(DBGFLAGS) -c c*.F
	$(FC) $(DBGFLAGS) -c d*.F e*.F f*.F
	$(FC) $(DBGFLAGS) -c g*.F h*.F
	$(FC) $(DBGFLAGS) -c i*.F l*.F m*.F n*.F
	$(FC) $(DBGFLAGS) -c o*.F
	$(FC) $(DBGFLAGS) -c p*.F q*.F r*.F s*.F t*.F u*.F v*.F w*.F y*.F z*.F
	ar rv $(SMK_BIN)/libmo5b.debug.a a*.o b*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a c*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a d*.o e*.o f*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a g*.o h*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a i*.o l*.o m*.o n*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a o*.o
	ar rv $(SMK_BIN)/libmo5b.debug.a p*.o q*.o r*.o s*.o t*.o u*.o \
                                         v*.o w*.o y*.o z*.o

clean:
	/bin/rm *.o
