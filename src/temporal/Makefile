# Version @(#)$Id$ 
# Path    $Source$
#
LIBS          = $(LDIRS) $(LIBS_SYS) 
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) 
LIBPATH       = $(SMK_BIN)/libsmoke.a
#LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.a
#

all: temporal

debug: temporal.debug

.f.o:
	
	$(FC) $(FFLAGS) -c -o $@ $< 

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
TMP= asgntpro.o genhemis.o mktmat.o opentmp.o opentmpin.o \
     setdaylt.o temporal.o tmnamunt.o updtmat.o wrtsup.o

#
temporal: $(TMP) $(ALLI) $(LIBPATH)
	echo TEMPORAL:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(TMP) $(LIBS)
#
#
temporal.debug: $(TMP) $(ALLI) $(LIBPATH_DEBUG)
	echo TEMPORAL \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(TMP) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc

asgntpro.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ asgntpro.f
genhemis.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODEMFAC) $(MODMET) $(MODINFO) \
            $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ genhemis.f
mktmat.o: $(MODTMPRL) 
	$(FC) $(FFLAGS) -c -o $@ mktmat.f
opentmp.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opentmp.f
opentmpin.o: $(MODMET) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ opentmpin.f
setdaylt.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ setdaylt.f
temporal.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODEMFAC) $(MODMET) $(MODINFO) \
            $(MODDAYHR) $(MODLISTS) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ temporal.f
tmnamunt.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ tmnamunt.f
updtmat.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ updtmat.f
#
clean:
	/bin/rm *.o
