# Version @(#)$Id$
# Path    $Source$
# Date    $Date$
#
FC        = f90

LIBS      = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBPATH   = $(SMK_BIN)/libsmoke.a   $(SMK_BIN)/libedsstools.a

##

all:   normbeis3v0.9 beis3v0.9 

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

NORM=normbeis3v0.9.o rdb3fac.o

normbeis3v0.9: $(NORM) $(ALLI) $(LIBPATH)
	echo NORMBEISV30.9:
	$(FC) $(FFLAGS) -o $(@) $(NORM) $(LIBS)


TMPB = tmpbeis3v0.9.o czangle.o getpar.o prebmet.o rdbpro.o hrbeis3.o 

beis3v0.9: $(TMPB) $(ALLI) $(LIBPATH)
	echo BEIS3V09:
	$(FC) $(FFLAGS) -o $(@) $(TMPB) $(LIBS)

#
## Module dependencies
#
include modmake.inc

# Modules
#

normbeis3v0.9.o: $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ $<

rdb3fac.o: $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ $<

tmpbeis3.o: $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ $<
hrbeis3.o: $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ $<
rdbpro.o:$(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $<

#
clean:
	/bin/rm *.o
