# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#

all: spcmat

debug: spcmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
SPC= spcmat.o opensmat.o asgnspro.o
     
#
spcmat: $(SPC) $(ALLI) $(LIBPATH)
	echo SPCMAT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(SPC) $(LIBS)
#
#
spcmat.debug: $(SPC) $(ALLI) $(LIBPATH_DEBUG)
	echo SPCMAT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(SPC) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc

asgnspro.o: $(MODSOURC) $(MODXREF) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
opensmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
spcmat.o: $(MODSOURC) $(MODSPRO) $(MODINFO) $(MODEMFAC)
	$(FC) $(FFLAGS) -c -o $@ $<


#
clean:
	/bin/rm *.o
