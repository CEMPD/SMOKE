# Version @(#)$Id$ $Source$ $Date$ 

lib:
	f77 $(FFLAGS) -c *.F
	cc -c *.c
	ar rv $(em)/libemstuff.a *.o
	rm *.o

debug:
	f77 -g -C -c *.F
	cc -g -C -c *.c
	ar rv $(em)/libemstuff.debug.a *.o
	rm *.o


