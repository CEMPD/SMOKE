# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBS_PVM      = $(LDIRS) -liocpl -lsmoke -ledsstools -lioapi_new -lnetcdf -lpvm3 -mp
LIBS_PVMDB    = $(LDIRS) -liocpl -lsmoke.debug -ledsstools.debug -lioapi_new -lnetcdf -lpvm3 -mp

LIBPATH       = $(SMK_BIN)/libsmoke.a   $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

#
all: rawbio tmpbio 

debug: rawbio.debug tmpbio.debug 

pvm: tmpbio.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/BIODIMS3.EXT


#
RBIO=rawbio.o  cybio.o grdbio.o rdbefac.o  
#
rawbio: $(RBIO) $(ALLI) $(LIBPATH)
	echo RAWBIO:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(RBIO) $(LIBS)
#
rawbio.debug: $(RBIO) $(ALLI) $(LIBPATH_DEBUG)
	echo RAWBIO \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(RBIO) $(LIBS_DEBUG)
#
TBIO=tmpbio.o czangle.o cldaten_calc.o hrbio.o hrbios.o prebmet.o rdbpro.o

tmpbio: $(TBIO) $(ALLI) $(LIBPATH)
	echo TMPBIO:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(TBIO) $(LIBS)

tmpbio.debug: $(TBIO) $(ALLI) $(LIBPATH_DEBUG)
	echo TMPBIO \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(TBIO) $(LIBS_DEBUG)

tmpbio.pvm: $(TBIO) $(ALLI) $(LIBPATH)
	echo TMPBIO \(pvm\):
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(TBIO) $(LIBS_PVM)

tmpbio.pvm.debug: $(TBIO) $(ALLI) $(LIBPATH_DEBUG)
	echo TMPBIO \(pvmdb\):
	$(FC) $(DBGFLAGS) -o $(SMK_BIN)/$(@) $(TBIO) $(LIBS_PVMDB)

include $(MD_SRC)/modmake.inc

rawbio.o: $(MODSURG) $(MODBIOG)
	$(FC) $(FFLAGS) -c -o $@ rawbio.f
cybio.o: $(MODSURG) $(MODBIOG)
	$(FC) $(FFLAGS) -c -o $@ cybio.f
grdbio.o: $(MODBIOG)
	$(FC) $(FFLAGS) -c -o $@ grdbio.f
tmpbio.o: $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ tmpbio.f
rdbpro.o:$(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ rdbpro.f
#
clean:
	/bin/rm *.o 
