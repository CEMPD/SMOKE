# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBS2         = $(LDIRS) -lsmoke -lmo5b -lioapi -lnetcdf
LIBS2_DEBUG   = $(LDIRS) -lsmoke.debug -lmo5b.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
##

all: premobl emisfac

debug: premobl.debug emisfac.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

PRE=premobl.o adjsmet.o applumat.o asgnpsi.o chkfulldy.o dyminmax.o \
    getdysav.o openpsiout.o opensmet.o setsrcdy.o \
    tmprinfo.o wrsmet.o genefmet.o setoutdt.o

premobl: $(PRE) $(ALLI) $(LIBPATH)
	echo PREMOBL:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(PRE) $(LIBS)

premobl.debug: $(PRE) $(ALLI) $(LIBPATH_DEBUG)
	echo 'PREMOBL (debug):'
	$(FC) $(FFLAGS) -o $(@) $(PRE) $(LIBS_DEBUG)
#
#
EMS=emisfac.o chkefhdr.o chkemfac.o efcombo.o efsetup.o emfacio.o expndeft.o \
    expndpsi.o mobil5b.o rdefdat.o rdeftmpr.o rdefupd.o seteftmpr.o setefupd.o \
    tmprinfo.o wremfacs.o

emisfac: $(EMS) $(ALLI) $(LIBPATH)
	echo EMISFAC:
	$(FC) $(FFLAGS) -I$(SMKROOT)/mo5b -o $(SMK_BIN)/$(@) $(EMS) $(LIBS2)

emisfac.debug: $(EMS) $(ALLI) $(LIBPATH_DEBUG)
	echo 'EMISFAC (debug):'
	$(FC) $(FFLAGS) -I$(SMKROOT)/mo5b -o $(@) $(EMS) $(LIBS2_DEBUG)
#
#
# Modules
#
MODSOURC = $(MD_SRC)/modsourc.o
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc

adjsmet.o  : $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
asgnpsi.o  : $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
chkefhdr.o: $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
chkemfac.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
efcombo.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
efsetup.o : $(MODEMFAC) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
emfacio.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
emisfac.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
expndeft.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
expndpsi.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
genefmet.o: $(MODXREF) $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
mobil5b.o : $(MODEMFACS) $(MODINFO) $(MODMET) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -I$(SMKROOT)/mo5b -c -o $@ $<
openpsiout.o  : $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
opensmet.o  : $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
premobl.o  : $(MODSOURC) $(MODXREF) $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
rdefdat.o : $(MODEMFAC) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
rdeftmpr.o : $(MODEMFAC) $(MODINFO) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
seteftmpr.o : $(MODEMFAC) $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
setefupd.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<
tmprinfo.o  : $(MODMET)
	$(FC) $(FFLAGS) -c -o $@ $<
wremfacs.o : $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm *.o
