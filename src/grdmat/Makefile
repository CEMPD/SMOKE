# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS       = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
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
GRD= asgnsurg.o genagmat.o genggmat.o genmgmat.o genpgmat.o grdmat.o \
     opengmat.o sizgmat.o setfrac.o rdsrg4grd.o
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
genagmat.o: $(MODSOURC) $(MODSURG) $(MODLISTS) $(MODINFO) $(MODGRID) $(MODXREF)
genggmat.o: $(MODSOURC) $(MODINFO) $(MODGRID) 
genmgmat.o: $(MODSOURC) $(MODSURG) $(MODLISTS) $(MODINFO) $(MODGRID) $(MODXREF)
grdmat.o: $(MODSOURC) $(MODXREF) $(MODINFO) $(MODGRID) $(MODFILESET)
opengmat.o: $(MODINFO)
setfrac.o: $(MODSURG)
sizgmat.o: $(MODSOURC) $(MODGRID) $(MODSURG)
rdsrg4grd.o: $(MODGRID) $(MODSURG)
#
clean:
	/bin/rm -f *.o 
