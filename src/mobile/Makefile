# Version "@(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: mbtest.debug
#
FC=f90
FFLAGS = $(DBGFLAGS)
LDFLAGS = -g
LDIRS = -L$(EDSS_LIB)  \
        -L$(EL_SRC)/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH= $(EL_SRC)/$(EDSS_OS)/libsmoke.debug.a
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
RAW=mbtest.o rdidamb.o
#
mbtest.debug: $(RAW) $(ALLI) $(LIBPATH)
	echo MBTEST:
	$(FC) $(FFLAGS) -o $(@) $(RAW) $(LIBS)
#
# Modules
#
MODSOURC = $(MD_SRC)/modsourc.o
#
## Module dependencies
#
$(MODSOURC):
	$(FC) $(FFLAGS) -c -o $@ $<
	mv modsourc.M $(MD_SRC)

mbtest.o  : $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<

rdidamb.o  : $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm *.o
