# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#
all: rawbio tmpbio 

debug: rawbio.debug tmpbio.debug 

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/BIODIMS3.EXT


#
RBIO=rawbio.o  rdbefac.o  
#
rawbio: $(RBIO) $(ALLI) $(LIBPATH)
	echo RAWBIO:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(RBIO) $(LIBS)
#
rawbio.debug: $(RBIO) $(ALLI) $(LIBPATH_DEBUG)
	echo RAWBIO \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(RBIO) $(LIBS_DEBUG)
#
TBIO=tmpbio.o czangle.o cldaten_calc.o hrbio.o hrbios.o prebmet.o rdbpro.o

tmpbio: $(TBIO) $(ALLI) $(LIBPATH)
	echo TMPBIO:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(TBIO) $(LIBS)

tmpbio.debug: $(TBIO) $(ALLI) $(LIBPATH_DEBUG)
	echo TMPBIO \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(TBIO) $(LIBS_DEBUG)


include $(MD_SRC)/modmake.inc

rawbio.o: $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ $<
tmpbio.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdbpro.o:$(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	/bin/rm *.o 
