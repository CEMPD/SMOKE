# Version @(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: spcmat.debug 
#
FC=f90
FFLAGS = $(DBGFLAGS)
LDFLAGS = -g
LDIRS = -L/env/proj/ioapi/IRIX6n32f90 \
        -L$(EDSS_LIB) \
        -L$(EL_SRC)/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH= $(EL_SRC)/$(EDSS_OS)/libsmoke.debug.a
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
SPC= spcmat.o opensmat.o asgnspro.o
     
#
spcmat.debug: $(SPC) $(ALLI) $(LIBPATH)
	echo SPCMAT:
	$(FC) $(FFLAGS) -o $(@) $(SPC) $(LIBS)

include $(MD_SRC)/modmake.inc

asgnspro.o: $(MODSOURC) $(MODXREF) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
opensmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
spcmat.o: $(MODSOURC) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<


#
clean:
	rm *.o
