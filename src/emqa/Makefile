# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a

all: getrecs 
#
debug: getrecs.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
MRG=getrecs.o
#
getrecs: $(MRG) $(ALLI) $(LIBPATH)
	echo SMKINVEN:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(RAW) $(LIBS)

getrecs.debug: $(MRG) $(ALLI) $(LIBPATH_DEBUG)
	echo GETRECS \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(MRG) $(LIBS_DEBUG)
#
#
MODSOURC = $(MD_SRC)/modsourc.o

include $(MD_SRC)/modmake.inc

getrecs.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
#
#
clean:
	/bin/rm *.o
