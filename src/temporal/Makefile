# Version @(#)$Id$ 
# Path    $Source$
#
LIBS          = $(LDIRS) $(LIBS_SYS)  $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS)  $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#

all: temporal

debug: temporal.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $< 

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
TMP= asgntpro.o chkemfac.o genhemis.o mktmat.o opentmp.o opentmpin.o \
     rdmetprof.o rdemfac.o setdaylt.o temporal.o tmnamunt.o updtmat.o wrtsup.o 
     

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
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgntpro.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODINFO)
chkemfac.o: $(MODEMFAC) $(MODINFO) $(MODLISTS) $(MODMBSET)
genhemis.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODEMFAC) $(MODMET) $(MODINFO) \
            $(MODDAYHR)
mktmat.o: $(MODTMPRL) 
opentmp.o: $(MODEMFAC) $(MODINFO) $(MODFILESET)
opentmpin.o: $(MODMET) $(MODINFO) $(MODGRID)
rdemfac.o: $(MODSOURC) $(MODINFO) $(MODMBSEP)
rdmetpro.o: $(MODXREF)
setdaylt.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
temporal.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODEMFAC) $(MODMET) $(MODINFO) \
            $(MODDAYHR) $(MODLISTS) $(MODGRID) $(MODMBSET)
tmnamunt.o: $(MODEMFAC) $(MODINFO)
updtmat.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODDAYHR)
wrtsup.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODDAYHR)
#
clean:
	/bin/rm -f *.o
