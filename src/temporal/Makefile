# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS)
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#

all: temporal

debug: temporal.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $< 

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
TMP= asgntpro.o genhemis.o mktmat.o opentmp.o opentmpin.o pretmpr.o rdefacs.o \
     setdaylt.o temporal.o tmnamunt.o tmpintrp.o updtmat.o wrtsup.o opentmp_u.o \
     rdustati.o getrnums.o
#qatmpr.o

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
opentmp_u.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opentmp_u.f
opentmpin.o: $(MODMET) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opentmpin.f
rdefacs.o: $(MODMOBIL) $(MODEMFAC) $(MODINFO) $(MODMET) $(MODXREF)
	$(FC) $(FFLAGS) -c -o $@ rdefacs.f
setdaylt.o: $(MODSOURC) $(MODSTCY) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ setdaylt.f
temporal.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODEMFAC) $(MODMET) $(MODINFO) \
            $(MODDAYHR) $(MODLISTS)
	$(FC) $(FFLAGS) -c -o $@ temporal.f
tmnamunt.o: $(MODEMFAC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ tmnamunt.f
tmpintrp.o: $(MODSOURC) $(MODMET) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ tmpintrp.f
qatmpr.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ qatmpr.f
updtmat.o: $(MODSOURC) $(MODXREF) $(MODTMPRL) $(MODDAYHR)
	$(FC) $(FFLAGS) -c -o $@ updtmat.f
rdustati.o: $(MODINFO) $(MODUNCERT)
	$(FC) $(FFLAGS) -c -o $@ rdustati.f
getrnums.o:  $(MODUNCERT)
	$(FC) $(FFLAGS) -c -o $@ getrnums.f
#
clean:
	/bin/rm *.o
