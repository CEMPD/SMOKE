# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS	  = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a
LIBS_PVM      = $(LDIRS) -liocpl -lsmoke -ledsstools -lioapi_new -lnetcdf -lpvm3 -mp
LIBS_PVMDB    = $(LDIRS) -liocpl -lsmoke.debug -ledsstools.debug -lioapi_new -lnetcdf -lpvm3 -mp

#

all: smkmerge mrggrid

debug: smkmerge.debug mrggrid.debug

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
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MGD) $(LIBS)
#
mrggrid.debug: $(MGD) $(ALLI) $(LIBPATH_DEBUG)
	echo MRGGRID \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(MGD) $(LIBS_DEBUG)
#
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc

allocmrg.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ allocmrg.f
bldmrgidx.o : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ bldmrgidx.f
getmrgev.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ getmrgev.f
initstcy.o  : $(MODMERGE) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ initstcy.f
mrgelev.o   : $(MODMERGE) $(MODCNTL) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ mrgelev.f
mrggrid.o   : $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ mrggrid.f
mrgmult.o   : $(MODMERGE) $(MODSTCY) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ mrgmult.f
mrgonams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ mrgonams.f
mrgunits.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ mrgunits.f
mrgvnams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ mrgvnams.f
openmrgin.o : $(MODMERGE) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ openmrgin.f
openmrgout.o: $(MODMERGE) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ openmrgout.f
rdmrginv.o: $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ rdmrginv.f
smkmerge.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODSURG) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ smkmerge.f
wmrgelev.o  : $(MODMERGE) $(MODSOURC) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ wmrgelev.f
wmrgemis.o  : $(MODMERGE) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ wmrgemis.f
wrmrgrep.o  : $(MODMERGE) $(MODSTCY) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ wrmrgrep.f


clean:
	/bin/rm *.o
