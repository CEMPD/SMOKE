# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS	   = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH    = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
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

asgnsurg.o: $(MODSOURC) $(MODXREF) $(MODSURG) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genmgmat.o: $(MODSOURC) $(MODXREF) $(MODSURG) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
grdmat.o: $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
opengmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
setfrac.o: $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
sizgmat.o: $(MODSOURC) $(MODXREF) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	/bin/rm *.o
