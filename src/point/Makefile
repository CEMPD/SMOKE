# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBS_PVM      = $(LDIRS) -liocpl -lsmoke -lioapi_new -lnetcdf -lpvm3 -mp
LIBS_PVMDB    = $(LDIRS) -liocpl -lsmoke.debug -lioapi_new -lnetcdf -lpvm3 -mp

LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#
all: laypoint elevpoint

debug: laypoint.debug elevpoint.debug

pvm: laypoint.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
LAY=laypoint.o preplm.o plmris.o postplm.o
#
laypoint: $(LAY) $(ALLI) $(LIBPATH)
	echo LAYPOINT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(LAY) $(LIBS)
#
laypoint.debug: $(LAY) $(ALLI) $(LIBPATH_DEBUG)
	echo LAYPOINT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(LAY) $(LIBS_DEBUG)

laypoint.pvm: $(LAY) $(ALLI) $(LIBPATH)
	echo LAYPOINT \(pvm\):
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(LAY) $(LIBS_PVM)
#
laypoint.pvm.debug: $(LAY) $(ALLI) $(LIBPATH_DEBUG)
	echo LAYPOINT \(pvmdb\):
	$(FC) $(FFLAGS) -o $(@) $(LAY) $(LIBS_PVMDB)

#
ELV=elevpoint.o openeout.o plumesti.o wpingstk.o
#
elevpoint: $(ELV) $(ALLI) $(LIBPATH)
	echo ELEVPOINT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(ELV) $(LIBS)
#
elevpoint.debug: $(ELV) $(ALLI) $(LIBPATH_DEBUG)
	echo ELEVPOINT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(ELV) $(LIBS_DEBUG)
#

#include $(MD_SRC)/modmake.inc
MD_OUT = $(MD_SRC)/$(SMOKE_EXE)
MODSOURC = $(MD_OUT)/modsourc.o
MODELEV  = $(MD_OUT)/modelev.o
MODINFO  = $(MD_OUT)/modinfo.o

laypoint.o: $(MODSOURC) $(MODELEV) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
elevpoint.o: $(MODSOURC) $(MODELEV) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openeout.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wpingstk.o: $(MODELEV)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	/bin/rm *.o 
