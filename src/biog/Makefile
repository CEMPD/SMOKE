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

all: normbeis3 tmpbeis3
debug: normbeis3.debug tmpbeis3.debug
pvm: tmpbeis3.pvm

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/BIODIMS3.EXT
ALL3=$(INC)/EMCNST3.EXT $(INC)/B3DIMS3.EXT

B3X = czangle.o getpar.o getparb.o hrbeis_360.o hrbeis3.o prebmet.o rdbpro.o tmpbeis3.o \
      tmpbeis314.o tmpbeis360.o  hrno.o

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
NORM3= normbeis3.o rdb3fac.o rdb4fac.o rdb4fac_csv.o  normbeis314.o normbeis360.o

normbeis3: $(NORM3) $(ALL3) $(LIBPATH)
	echo NORMBEIS3:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(NORM3) $(LIBS)

normbeis3.debug: $(NORM3) $(ALL3) $(LIBPATH_DEBUG)
	echo NORMBEIS3 \(debug\):
	$(FC) $(DBGFLAGS) -o $(@) $(NORM3) $(LIBS_DEBUG)

include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

normbeis3.o:  $(MODBEIS3)
tmpbeis3.o: $(MODSPRO) $(MODBEIS3)
rdbpro.o:$(MODSPRO)
#
clean:
	/bin/rm -f *.o  
