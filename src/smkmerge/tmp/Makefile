# Version @(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: smkmerge.debug 
#
FC=f90
FFLAGS = $(DBGFLAGS)
LDFLAGS = -g
LDIRS = -L/env/proj/ioapi/IRIX6n32f90 \
        -L$(EDSS_LIB) \
        -L$(EL_SRC)/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH= $(EL_SRC)/$(EDSS_OS)/libsmoke.debug.a
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
MRG= smkmerge.o allocmrg.o bldmrgidx.o getmrgev.o mrgmult.o mrgonams.o \
     mrgvnams.o openmrgin.o openmrgout.o wmrgemis.o wrmrgrep.o initstcy.o \
     rdmrginv.o mrgping.o
#
smkmerge.debug: $(MRG) $(ALLI) $(LIBPATH)
	echo SMKMERGE:
	$(FC) $(FFLAGS) -o $(@) $(MRG) $(LIBS)
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc

allocmrg.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
bldmrgidx.o : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
getmrgev.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
initstcy.o  : $(MODMERGE) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgmult.o   : $(MODMERGE) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgonams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgping.o   : $(MODMERGE) $(MODCNTL) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
mrgvnams.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
openmrgin.o : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
openmrgout.o: $(MODMERGE) $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
rdmrginv.o: $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
smkmerge.o  : $(MODMERGE) $(MODCNTRL) $(MODSTCY) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
wmrgemis.o  : $(MODMERGE)
	$(FC) $(FFLAGS) -c -o $@ $<
wrmrgrep.o  : $(MODMERGE) $(MODSTCY)
	$(FC) $(FFLAGS) -c -o $@ $<
