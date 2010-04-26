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

all: movesmrg

debug: movesmrg.debug

pvm: movesmrg.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
MRG= movesmrg.o bldmrgidx.o bldprocidx.o bldsrccell.o \
     getmrgev.o initstcy.o mrgonams.o mrgunits.o \
     mrgvnams.o openmrgin.o openmrgout.o rdfmref.o \
     rdmetmoves.o rdmrclist.o rdmxref.o \
     rdrpdemfacs.o rdrppemfacs.o rdrpvemfacs.o \
     setrefcnty.o wrmrgrep.o

#
movesmrg: $(MRG) $(ALLI) $(LIBPATH)
	echo MOVESMRG:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MRG) $(LIBS)
#
movesmrg.debug: $(MRG) $(ALLI) $(LIBPATH_DEBUG)
	echo MOVESMRG \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(MRG) $(LIBS_DEBUG)
#
#
movesmrg.pvm: $(MRG) $(ALLI) $(LIBPATH)
	echo MOVESMRG \(pvm\):
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MRG) $(LIBS_PVM)
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

bldmrgidx.o : $(MODMERGE)
bldprocidx.o: $(MODMVSMRG) $(MODMERGE)
bldsrccell.o: $(MODMVSMRG)
getmrgev.o  : $(MODINFO) $(MODMERGE) $(MODMVSMRG)
initstcy.o  : $(MODMERGE) $(MODSTCY) $(MODSOURC)
movesmrg.o  : $(MODMERGE) $(MODMVSMRG) $(MODLISTS) $(MODSTCY) $(MODGRID) $(MODMBSET) $(MODSOURC)
mrgonams.o  : $(MODMERGE)
mrgunits.o  : $(MODMERGE)
mrgvnams.o  : $(MODMERGE) $(MODLISTS)
openmrgin.o : $(MODMERGE) $(MODMVSMRG) $(MODINFO) $(MODSOURC) $(MODEMFAC) $(MODGRID) $(MODFILESET)
openmrgout.o: $(MODMERGE) $(MODGRID) $(MODFILESET)
rdfmref.o   : $(MODMBSET)
rdmetmoves.o: $(MODMVSMRG) $(MODLISTS)
rdmrclist.o : $(MODMBSET) $(MODMVSMRG)
rdmxref.o   : $(MODMBSET)
rdrpdemfacs.o: $(MODMBSET) $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrppemfacs.o: $(MODMBSET) $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrpvemfacs.o: $(MODMBSET) $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
setrefcnty.o: $(MODMVSMRG) $(MODINFO) $(MODSOURC) $(MODLISTS) $(MODMBSET)
wrmrgrep.o  : $(MODMERGE) $(MODSTCY) $(MODGRID) 

clean:
	/bin/rm -f *.o
