# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBS2         = $(LDIRS) $(LIBS_SYS) -lmo6 $(LIBS_COMPILE)
LIBS2_DEBUG   = $(LDIRS) $(LIBS_DEBUG_SYS) -lmo6.debug $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a  $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#

all: mbsetup premobl emisfac

debug: mbsetup.debug premobl.debug emisfac.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

MBSET=mbsetup.o rdmcref.o rdmvref.o wrspdsum.o wrtimegr.o asgnspds.o

mbsetup: $(MBSET) $(ALLI) $(LIBPATH)
	echo MBSETUP:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MBSET) $(LIBS)

mbsetup.debug: $(MBSET) $(ALLI) $(LIBPATH_DEBUG)
	echo 'MBSETUP (debug):'
	$(FC) $(FFLAGS) -o $(@) $(MBSET) $(LIBS_DEBUG)

#
#
PRE=premobl.o rdspdsrc.o openshour.o asgngrps.o \
      hourtemp.o avertemp.o wrshour.o opengroup.o \
      rdgroups.o

premobl: $(PRE) $(ALLI) $(LIBPATH)
	echo 'PREMOBL:'
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(PRE) $(LIBS)

premobl.debug: $(PRE) $(ALLI) $(LIBPATH_DEBUG)
	echo 'PREMOBL (debug):'
	$(FC) $(FFLAGS) -o $(@) $(PRE) $(LIBS_DEBUG)

#
#
EMS=emisfac.o rdm6list.o wrm6input.o rdgrplist.o \
     openscen.o rdspdline.o chkm6scn.o getspdln.o \
     wrspdvmt.o rdhourtemp.o rplctemp.o wrm6header.o \
     opensef.o wremfacs.o rplcmet.o calcrelhum.o

emisfac: $(EMS) $(ALLI) $(LIBPATH)
	echo 'EMISFAC:'
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(EMS) $(LIBS2)

emisfac.debug: $(EMS) $(ALLI) $(LIBPATH_DEBUG)
	echo 'EMISFAC (debug):'
	$(FC) $(FFLAGS) -o $(@) $(EMS) $(LIBS2_DEBUG)

#
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgngrps.o    : $(MODMBSET) $(MODMET) $(MODINFO)
asgnspds.o    : $(MODSOURC) $(MODXREF) $(MODINFO)
avertemp.o    : $(MODMET)
emisfac.o     : $(MODINFO) $(MODLISTS) $(MODMET) $(MODEMFAC) $(MODFILESET)
hourtemp.o    : $(MODMBSET) $(MODSOURC) $(MODINFO) $(MODMET)
mbsetup.o     : $(MODSOURC) $(MODINFO) $(MODMBSET)
openscen.o    : $(MODEMFAC)
opensef.o     : $(MODINFO) $(MODEMFAC) $(MODFILESET)
openshour.o   : $(MODINFO) $(MODMET)
premobl.o     : $(MODSOURC) $(MODINFO) $(MODMET) $(MODGRID) $(MODMBSET)
rdgroups.o    : $(MODINFO)
rdhourtemp.o  : $(MODMET)
rdm6list.o    : $(MODEMFAC)
rdmcref.o     : $(MODMBSET)
rdmvref.o     : $(MODMBSET)
rdspdline.o   : $(MODEMFAC)
rplcmet.o     : $(MODMET)
rplctemp.o    : $(MODMET)
wremfacs.o    : $(MODSOURC) $(MODINFO) $(MODEMFAC)
wrm6header.o  : $(MODEMFAC)
wrspdsum.o    : $(MODSOURC) $(MODINFO) $(MODMBSET) $(MODXREF)
wrspdvmt.o    : $(MODMOBIL)
wrtimegr.o    : $(MODMBSET)

clean:
	/bin/rm -f *.o
ifcclean:
	/bin/rm -f *.o  
	/bin/rm -f *.d  
	/bin/rm -f work.pc*
	/bin/rm -f $(MODALL)
	/bin/rm -f $(IOMODALL)
	/bin/rm -f $(FSMODALL)
