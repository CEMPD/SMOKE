# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBS2         = $(LDIRS) -lsmoke -lmo5b -ledsstools  -lioapi -lnetcdf
LIBS2_DEBUG   = $(LDIRS) -lsmoke.debug -lmo5b.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBS3         = $(LDIRS) -lsmoke -lmo6 -ledsstools -lioapi -lnetcdf
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
FFLAGS_NOOPT  = $(DBGFLAGS)
##

all: mbsetup premobl emisfac

debug: mbsetup.debug premobl.debug emisfac.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

MBSET=mbsetup.o rdmcref.o rdmvref.o wrspdsum.o wrtimegr.o

mbsetup: $(MBSET) $(ALLI) $(LIBPATH)
	echo MBSETUP:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MBSET) $(LIBS)

mbsetup.debug: $(MBSET) $(ALLI) $(LIBPATH_DEBUG)
	echo 'MBSETUP (debug):'
	$(FC) $(FFLAGS) -o $(@) $(MBSET) $(LIBS_DEBUG)

#
#
PRE=premobl.o adjsmet.o chkfulldy.o dyminmax.o \
    getdysav.o openpsiout.o opensmet.o setsrcdy.o \
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
PRE2=premobl_new.o rdspdsrc.o getctynum.o chkfulldy.o \
     openshour.o setoutdt.o hourtemp.o avertemp.o \
     wrshour.o opengroup.o rdgroups.o asgngrps.o

premobl_new: $(PRE2) $(ALLI) $(LIBPATH)
	echo 'PREMOBL (new):'
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(PRE2) $(LIBS)

#
#
EMS2=emisfac_new.o rdm6list.o wrm6input.o rdgrplist.o \
     openscen.o rdspdline.o chkscnyr.o getspdln.o \
     wrspdvmt.o rdhourtemp.o rplctemp.o wrm6header.o \
     opensef.o wremfacs.o

emisfac_new: $(EMS2) $(ALLI) $(LIBPATH)
	echo 'EMISFAC (new):'
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(EMS2) $(LIBS3)

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

adjsmet.o     : $(MODSOURC) $(MODINFO)          # This should not get optimization
adjshour.o    : $(MODSOURC)
chkefhdr.o    : $(INC)/M5CNST3.EXT
chkemfac.o    : $(MODEMFAC)
chkscnyr.o    : $(MODINFO)
efcombo.o     : $(MODEMFAC) $(MODINFO) $(MODMET)
emfacio.o     : $(MODEMFAC) $(MODINFO) $(MODMET)
emisfac.o     : $(MODEMFAC) $(MODINFO) $(MODMET)
emisfac_new.o : $(MODSOURC) $(MODINFO) $(MODMBSET)
expndeft.o    : $(MODEMFAC)
expndpsi.o    : $(MODEMFAC)
genefmet.o    : $(MODXREF) $(MODEMFAC)
mbsetup.o     : $(MODSOURC) $(MODINFO) $(MODMBSET)
mobil5b.o     : $(MODEMFACS) $(MODINFO) $(MODMET) $(INC)/M5CNST3.EXT
	$(FC) $(FFLAGS) -I$(M5_SRC) -c -o $@ mobil5b.f
openscen.o    : $(MODMBSET)
opensef.o     : $(MODINFO)
opensmet.o    : $(MODINFO) $(MODMET)
openshour.o   : $(MODINFO) $(MODMET)
premobl.o     : $(MODSOURC) $(MODXREF) $(MODEMFAC) $(MODINFO) $(MODMET)
premobl_new.o : $(MODSOURC) $(MODINFO) $(MODMET) $(MODMBSET)
rdefdat.o     : $(MODEMFAC) $(INC)/M5CNST3.EXT
rdeftmpr.o    : $(MODEMFAC) $(MODINFO) $(MODMET)
rdm6list.o    : $(MODMBSET)
rdmcref.o     : $(MODSOURC) $(MODINFO) $(MODLISTS) $(MODMBSET)
rdmvref.o     : $(MODSOURC) $(MODMBSET)
rdspdline.o   : $(MODMBSET)
seteftmpr.o   : $(MODEMFAC) $(MODMET)
setefupd.o    : $(MODEMFAC)
#wremfacs.o    : $(MODEMFAC)
wremfacs.o    : $(MODSOURC) $(MODINFO) $(MODMBSET)
wrm6input.o   : $(MODSOURC) $(MODINFO) $(MODMBSET)
wrspdsum.o    : $(MODSOURC) $(MODINFO) $(MODMBSET)
wrspdvmt.o    : $(MODMBSET)
wrtimegr.o    : $(MODSOURC) $(MODINFO) $(MODMBSET)

clean:
	/bin/rm *.o
