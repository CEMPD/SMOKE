# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS) 
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG)
LIBS_PVM      = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp
LIBS_PVMDB    = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp

LIBPATH       = $(SMK_BIN)/libsmoke.a   $(TOOLS_BIN)/libedsstools.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

#  BEIS3 executables : beis3 and normbeis3
#  BEIS2 executables : rawbio and tmpbio

all: beis3 normbeis3 rawbio tmpbio 

debug: beis3.debug normbeis3.debug rawbio.debug tmpbio.debug 

pvm: tmpbio.pvm


.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/BIODIMS3.EXT
ALL3=$(INC)/EMCNST3.EXT $(INC)/B3DIMS3.EXT

B3X = czangle.o getpar.o hrbeis3.o prebmet.o rdbpro.o tmpbeis3.o 

beis3: $(B3X) $(ALL3) $(LIBPATH)
	echo BEIS3:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(B3X) $(LIBS)

beis3.debug: $(B3X) $(ALL3) $(LIBPATH_DEBUG)
	echo BEIS3 \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(B3X) $(LIBS_DEBUG)

#
NORM3= normbeis3.o rdb3fac.o

normbeis3: $(NORM3) $(ALL3) $(LIBPATH)
	echo NORMBEIS3:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(NORM3) $(LIBS)

normbeis3.debug: $(NORM3) $(ALL3) $(LIBPATH_DEBUG)
	echo NORMBEIS3 \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(NORM3) $(LIBS_DEBUG)
#
RBIO=rawbio.o cybio.o grdbio.o rdbefac.o  
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


beis3.o:  $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ beis3.f
normbeis3.o:  $(MODBEIS3)
	$(FC) $(FFLAGS) -c -o $@ normbeis3.f    
rawbio.o: $(MODSURG) $(MODBIOG) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ rawbio.f
cybio.o: $(MODSURG) $(MODBIOG)
	$(FC) $(FFLAGS) -c -o $@ cybio.f
grdbio.o: $(MODBIOG) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ grdbio.f
hrbio.o: $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ hrbio.f
hrbios.o: $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ hrbios.f
tmpbio.o: $(MODSPRO) $(MODGRID)
	$(FC) $(FFLAGS) -c -o $@ tmpbio.f
rdbpro.o:$(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ rdbpro.f
#
clean:
	/bin/rm *.o 
