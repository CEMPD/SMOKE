# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBS2         = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG2   = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a

LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a
LIBS_PVM      = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp
LIBS_PVMDB    = $(LDIRS) -liocpl $(LIBS_DEBUG) -lpvm3 -mp

#

all: smkmerge mrggrid mrgelev

debug: smkmerge.debug mrggrid.debug mrgelev.debug

pvm: smkmerge.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
MRG= smkmerge.o allocmrg.o bldmrgidx.o getmrgev.o mrgmult.o mrgonams.o \
     mrgvnams.o openmrgin.o openmrgout.o wmrgemis.o wrmrgrep.o initstcy.o \
     rdmrginv.o mrgelev.o mrgunits.o mrgbio.o wmrgelev.o
#
smkmerge: $(MRG) $(ALLI) $(LIBPATH)
	echo SMKMERGE:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MRG) $(LIBS)
#
smkmerge.debug: $(MRG) $(ALLI) $(LIBPATH_DEBUG)
	echo SMKMERGE \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(MRG) $(LIBS_DEBUG)
#
#
smkmerge.pvm: $(MRG) $(ALLI) $(LIBPATH)
	echo SMKMERGE \(pvm\):
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MRG) $(LIBS_PVM)
#
MGD= mrggrid.o
#
mrggrid: $(MGD) $(ALLI) $(LIBPATH)
	echo MRGGRID:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MGD) $(LIBS2)
#
mrggrid.debug: $(MGD) $(ALLI) $(LIBPATH_DEBUG)
	echo MRGGRID \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(MGD) $(LIBS_DEBUG2)
#
MGE= mrgasciielev.o
#
mrgelev: $(MGE) $(ALLI) $(LIBPATH)
	echo MRGELEV:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MGE) $(LIBS2)
#
mrgelev.debug: $(MGE) $(ALLI) $(LIBPATH_DEBUG)
	echo MRGELEV \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(MGE) $(LIBS_DEBUG2)
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

allocmrg.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODELEV) $(MODGRID) 
bldmrgidx.o : $(MODMERGE)
getmrgev.o  : $(MODMERGE)
initstcy.o  : $(MODMERGE) $(MODSTCY)
mrgelev.o   : $(MODMERGE) $(MODCNTL) $(MODELEV)
mrggrid.o   : $(MODGRID)
mrgmult.o   : $(MODMERGE) $(MODSTCY) $(MODELEV)
mrgonams.o  : $(MODMERGE)
mrgunits.o  : $(MODMERGE)
mrgvnams.o  : $(MODMERGE)
openmrgin.o : $(MODMERGE) $(MODELEV) $(MODGRID) $(MODFILESET) $(MODINFO)
openmrgout.o: $(MODMERGE) $(MODELEV) $(MODGRID) 
rdmrginv.o: $(MODMERGE)
smkmerge.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODSURG) $(MODELEV) $(MODGRID) $(MODLISTS)
wmrgelev.o  : $(MODMERGE) $(MODSOURC) $(MODELEV) $(MODGRID) 
wmrgemis.o  : $(MODMERGE) $(MODELEV)
wrmrgrep.o  : $(MODMERGE) $(MODSTCY) $(MODGRID) 

clean:
	/bin/rm -f *.o
