# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS	   = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH    = $(SMK_BIN)/libsmoke.a $(TOOLS_BIN)/libedsstools.a 
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a 
#

all: grdmat

debug: grdmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
GRD= asgnsurg.o genagmat.o genmgmat.o genpgmat.o grdmat.o lnk2grd.o opengmat.o \
     sizgmat.o setfrac.o
#
grdmat: $(GRD) $(ALLI) $(LIBPATH)
	echo GRDMAT
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(GRD) $(LIBS)
#
#
grdmat.debug: $(GRD) $(ALLI) $(LIBPATH_DEBUG)
	echo GRDMAT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(GRD) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

asgnsurg.o: $(MODSOURC) $(MODXREF) $(MODSURG) $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ asgnsurg.f
genagmat.o: $(MODSOURC) $(MODXREF) $(MODSURG) $(MODLISTS) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ genagmat.f
genmgmat.o: $(MODSOURC) $(MODXREF) $(MODSURG) $(MODLISTS) $(MODINFO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ genmgmat.f
grdmat.o: $(MODSOURC) $(MODXREF) $(MODINFO) $(MODGRID) $(MODFILESET)
	$(FC) $(FFLAGS) -c -o $@ grdmat.f
opengmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ opengmat.f
setfrac.o: $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ setfrac.f
sizgmat.o: $(MODSOURC) $(MODXREF) $(MODGRID) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ sizgmat.f
#
clean:
	/bin/rm -f *.o 
ifcclean:
	/bin/rm -f *.o  
	/bin/rm -f *.d  
	/bin/rm -f work.pc*
	/bin/rm -f $(MODALL)
	/bin/rm -f $(IOMODALL)
	/bin/rm -f $(FSMODALL)
