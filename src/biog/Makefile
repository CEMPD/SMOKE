# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBS_PVM      = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp $(LIBS_COMPILE)
LIBS_PVMDB    = $(LDIRS) -liocpl $(LIBS) -lpvm3 -mp $(LIBS_COMPILE)

LIBPATH       = $(SMK_BIN)/libsmoke.a   $(TOOLS_BIN)/libedsstools.a 
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a $(TOOLS_BIN)/libedsstools.debug.a

#  BEIS3 executables : normbeis3 and tmpbeis3
#  BEIS2 executables : rawbio and tmpbio

all: normbeis3 tmpbeis3 rawbio tmpbio 
debug: normbeis3.debug tmpbeis3.debug rawbio.debug tmpbio.debug 
pvm: tmpbeis3.pvm tmpbio.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/BIODIMS3.EXT
ALL3=$(INC)/EMCNST3.EXT $(INC)/B3DIMS3.EXT

B3X = czangle.o getpar.o hrbeis3.o prebmet.o rdbpro.o tmpbeis3.o \
      tmpbeis309.o tmpbeis312.o hrbeis309.o hrno.o

tmpbeis3: $(B3X) $(ALL3) $(LIBPATH)
	echo TMPBEIS3:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(B3X) $(LIBS)

tmpbeis3.debug: $(B3X) $(ALL3) $(LIBPATH_DEBUG)
	echo TMPBEIS3 \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(B3X) $(LIBS_DEBUG)

tmpbeis3.pvm: $(B3X) $(ALLI) $(LIBPATH)
	echo TMPBEIS3 \(pvm\):
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(B3X) $(LIBS_PVM)

#
NORM3= normbeis3.o rdb3fac.o normbeis309.o normbeis312.o

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
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

normbeis3.o:  $(MODBEIS3)
rawbio.o: $(MODSURG) $(MODBIOG) $(MODGRID)
cybio.o: $(MODSURG) $(MODBIOG)
grdbio.o: $(MODBIOG) $(MODGRID)
hrbio.o: $(MODGRID)
hrbios.o: $(MODGRID)
tmpbeis3.o: $(MODSPRO) $(MODBEIS3)
tmpbio.o: $(MODSPRO) $(MODGRID)
rdbpro.o:$(MODSPRO)
#
clean:
	/bin/rm -f *.o  
