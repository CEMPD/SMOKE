# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS	  = $(LDIRS) -lsmoke -liotools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -liotools.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libiotools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libiotools.debug.a
#

all: smkmerge mrggrid

debug: smkmerge.debug mrggrid.debug

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
include $(IOMOD_SRC)/modmake.inc

allocmrg.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<
bldmrgidx.o : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
getmrgev.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
initstcy.o  : $(MODMERGE) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgelev.o   : $(MODMERGE) $(MODCNTL) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgmult.o   : $(MODMERGE) $(MODSTCY) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgonams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgunits.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgvnams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
openmrgin.o : $(MODMERGE) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<
openmrgout.o: $(MODMERGE) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<
rdmrginv.o: $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
smkmerge.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODSURG) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<
wmrgelev.o  : $(MODMERGE) $(MODSOURC) $(MODELEV) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<
wmrgemis.o  : $(MODMERGE) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
wrmrgrep.o  : $(MODMERGE) $(MODSTCY) $(MODGRID) 
	$(FC) $(FFLAGS) -c -o $@ $<


clean:
	/bin/rm *.o
