# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG)
LIBS2         = $(LDIRS) $(LIBS) -lmo5b
LIBS2_DEBUG   = $(LDIRS) $(LIBS_DEBUG) -lmo5b.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
FFLAGS_NOOPT  = $(DBGFLAGS)
##

all: premobl emisfac

debug: premobl.debug emisfac.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

PRE=premobl.o adjsmet.o chkfulldy.o dyminmax.o openpsiout.o opensmet.o \
    wrsmet.o genefmet.o setoutdt.o

# NOTE: the optimized version of premobl does not work

premobl: $(PRE) $(ALLI) $(LIBPATH)
	echo PREMOBL:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(PRE) $(LIBS)

premobl.debug: $(PRE) $(ALLI) $(LIBPATH_DEBUG)
	echo 'PREMOBL (debug):'
	$(FC) $(FFLAGS) -o $(@) $(PRE) $(LIBS_DEBUG)
#
#
EMS=emisfac.o chkefhdr.o chkemfac.o efcombo.o emfacio.o expndeft.o \
    expndpsi.o mobil5b.o rdefdat.o rdeftmpr.o rdefupd.o seteftmpr.o setefupd.o \
    wremfacs.o

emisfac: $(EMS) $(ALLI) $(LIBPATH)
	echo EMISFAC:
	$(FC) $(FFLAGS) -I$(M5_SRC) -o $(SMK_BIN)/$(@) $(EMS) $(LIBS2)

emisfac.debug: $(EMS) $(ALLI) $(LIBPATH_DEBUG)
	echo 'EMISFAC (debug):'
	$(FC) $(FFLAGS) -I$(M5_SRC) -o $(@) $(EMS) $(LIBS2_DEBUG)
#
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc

adjsmet.o  : $(MODSOURC) $(MODINFO)          # This should not get optimization
	$(FC) $(FFLAGS) -c -o $@ adjsmet.f
chkefhdr.o: $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ chkefhdr.f
chkemfac.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ chkemfac.f
efcombo.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ efcombo.f
emfacio.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ emfacio.f
emisfac.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ emisfac.f
expndeft.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ expndeft.f
expndpsi.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ expndpsi.f
genefmet.o: $(MODXREF) $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ genefmet.f
mobil5b.o : $(MODEMFACS) $(MODINFO) $(MODMET) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -I$(M5_SRC) -c -o $@ mobil5b.f
opensmet.o  : $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ opensmet.f
premobl.o  : $(MODSOURC) $(MODXREF) $(MODEMFAC) $(MODINFO) $(MODMET) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ premobl.f
rdefdat.o : $(MODEMFAC) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ rdefdat.f
rdeftmpr.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ rdeftmpr.f
seteftmpr.o : $(MODEMFAC) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ seteftmpr.f
setefupd.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ setefupd.f
wremfacs.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ wremfacs.f

clean:
	/bin/rm *.o
