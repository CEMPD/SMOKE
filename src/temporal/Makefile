# Version @(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: temporal.debug 
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
TMP= temporal.o asgntpro.o genhemis.o opentmp.o
     
#
temporal.debug: $(TMP) $(ALLI) $(LIBPATH)
	echo TEMPORAL:
	$(FC) $(FFLAGS) -o $(@) $(TMP) $(LIBS)

include $(MD_SRC)/modmake.inc

asgntpro.o: $(MODSOURC) $(MODXREF) $(MODTPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genhemis.o: $(MODXREF) $(MODTPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
opentmp.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
temporal.o: $(MODSOURC) $(MODXREF) $(MODTPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	rm *.o
