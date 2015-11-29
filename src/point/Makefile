# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBS_PVM      = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp $(LIBS_COMPILE)
LIBS_PVMDB    = $(LDIRS) -liocpl $(LIBS_DEBUG) -lpvm3 -mp $(LIBS_COMPILE)

LIBPATH       = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a
#
all: laypoint elevpoint

debug: laypoint.debug elevpoint.debug

pvm: laypoint.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
LAY=laypoint.o openlayout.o preplm.o plmris.o postplm.o plsprd.o \
    fire_preplm.o fire_plmris.o fire_postplm.o 

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
ELV=elevpoint.o asgngrps.o openeout.o plumesti.o wpingstk.o mxgrpemis.o \
    rpelvcfg.o 

#
elevpoint: $(ELV) $(ALLI) $(LIBPATH)
	echo ELEVPOINT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(ELV) $(LIBS)
#
elevpoint.debug: $(ELV) $(ALLI) $(LIBPATH_DEBUG)
	echo ELEVPOINT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(ELV) $(LIBS_DEBUG)
#

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgngrps.o: $(MODSOURC) $(MODELEV) $(MODINFO)
laypoint.o: $(MODSOURC) $(MODELEV) $(MODGRID) $(MODINFO)
mxgrpemis.o: $(MODSOURC) $(MODELEV) $(MODINFO) $(MODFILESET)
elevpoint.o: $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID)
openeout.o: $(MODSOURC) $(MODINFO)
openlayout.o: $(MODELEV) $(MODGRID) $(MODINFO)
rpelvcfg.o: $(MODREPRT) $(MODELEV) $(MODINFO)
wpingstk.o: $(MODELEV)
#
clean:
	/bin/rm -f *.o 
