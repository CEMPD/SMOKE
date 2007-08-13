# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) $(LIBS_SYS) $(LIBS_COMPILE)
LIBS_DEBUG    = $(LDIRS) $(LIBS_DEBUG_SYS) $(LIBS_COMPILE)
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
##

all: aggwndw beld3to2 bluesky2inv cemscan edms2inv extractida geofac invsplit metcombine metscan pktreduc \
     smk2emis surgtool uam2ncf 

debug: aggwndw.debug beld3to2.debug bluesky2inv.debug cemscan.debug edms2inv.debug extractida.debug invsplit.debug \
       geofac.debug metcombine.debug metscan.debug pktreduc.debug smk2emis.debug \
       surgtool.debug uam2ncf.debug  

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	


ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

AGGW=aggwndw.o

aggwndw: $(AGGW) $(ALLI) $(LIBPATH)
	echo AGGWNDW:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(AGGW) $(LIBS)

aggwndw.debug: $(AGGW) $(ALLI) $(LIBPATH_DEBUG)
	echo 'AGGWNDW (debug):'
	$(FC) $(DBGFLAGS) -o $(@) $(AGGW) $(LIBS_DEBUG)

B3TO2=beld3to2.o

beld3to2: $(B3TO2) $(ALLI) $(LIBPATH)
	echo BELD3TO2:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(B3TO2) $(LIBS)

beld3to2.debug: $(B3TO2) $(ALLI) $(LIBPATH_DEBUG)
	echo 'BELD3TO2 (debug):'
	$(FC) $(FFLAGS) -o $(@) $(B3TO2) $(LIBS_DEBUG)

BL2IV=bluesky2inv.o

bluesky2inv: $(BL2IV) $(ALLI) $(LIBPATH)
	echo BlueSky2Inv:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(BL2IV) $(LIBS)

bluesky2inv.debug: $(BL2IV) $(ALLI) $(LIBPATH_DEBUG)
	echo 'BlueSky2Inv (debug):'
	$(FC) $(FFLAGS) -o $(@) $(BL2IV) $(LIBS_DEBUG)

CEM=cemscan.o

cemscan: $(CEM) $(ALLI) $(LIBPATH)
	echo CEMScan:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(CEM) $(LIBS)

cemscan.debug: $(CEM) $(ALLI) $(LIBPATH_DEBUG)
	echo 'CEMScan (debug):'
	$(FC) $(FFLAGS) -o $(@) $(CEM) $(LIBS_DEBUG)

EM2IV=edms2inv.o

edms2inv: $(EM2IV) $(ALLI) $(LIBPATH)
	echo Edms2Inv:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(EM2IV) $(LIBS)

edms2inv.debug: $(EM2IV) $(ALLI) $(LIBPATH_DEBUG)
	echo 'Edms2Inv (debug):'
	$(FC) $(FFLAGS) -o $(@) $(EM2IV) $(LIBS_DEBUG)

EXTI=extractida.o

extractida: $(EXTI) $(ALLI) $(LIBPATH)
	echo EXTRACTIDA:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(EXTI) $(LIBS)

extractida.debug: $(EXTI) $(ALLI) $(LIBPATH_DEBUG)
	echo 'EXTRACTIDA (debug):'
	$(FC) $(FFLAGS) -o $(@) $(EXTI) $(LIBS_DEBUG)

GEOF=geofac.o

geofac: $(GEOF) $(ALLI) $(LIBPATH)
	echo GEOFAC:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(GEOF) $(LIBS)

geofac.debug: $(GEOF) $(ALLI) $(LIBPATH_DEBUG)
	echo 'GEOFAC (debug):'
	$(FC) $(FFLAGS) -o $(@) $(GEOF) $(LIBS_DEBUG)

INVS=invsplit.o

invsplit: $(INVS) $(ALLI) $(LIBPATH)
	echo INVSPLIT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(INVS) $(LIBS)

invsplit.debug: $(INVS) $(ALLI) $(LIBPATH_DEBUG)
	echo 'INVSPLIT (debug):'
	$(FC) $(FFLAGS) -o $(@) $(INVS) $(LIBS_DEBUG)

MCMB=metcombine.o

metcombine: $(MCMB) $(ALLI) $(LIBPATH)
	echo METCOMBINE:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MCMB) $(LIBS)

metcombine.debug: $(MCMB) $(ALLI) $(LIBPATH_DEBUG)
	$(FC) $(DBGFLAGS) -o $(@) $(MCMB) $(LIBS_DEBUG)

MSCN=metscan.o
#
metscan: $(MSCN) $(ALLI) $(LIBPATH)
	echo METSCAN:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MSCN) $(LIBS)
#
metscan.debug: $(MSCN) $(ALLI) $(LIBPATH_DEBUG)
	$(FC) $(DBGFLAGS) -o $(@) $(MSCN) $(LIBS_DEBUG)

PKT=pktreduc.o

pktreduc: $(PKT) $(ALLI) $(LIBPATH)
	echo PKTREDUC:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(PKT) $(LIBS)

pktreduc.debug: $(PKT) $(ALLI) $(LIBPATH_DEBUG)
	echo 'PKTREDUC (debug):'
	$(FC) $(FFLAGS) -o $(@) $(PKT) $(LIBS_DEBUG)

SMKE=smk2emis.o

smk2emis: $(SMKE) $(ALLI) $(LIBPATH)
	echo SMK2EMIS:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(SMKE) $(LIBS)

smk2emis.debug: $(SMKE) $(ALLI) $(LIBPATH_DEBUG)
	echo 'SMK2EMIS (debug):'
	$(FC) $(FFLAGS) -o $(@) $(SMKE) $(LIBS_DEBUG)

SURG=surgtool.o

surgtool: $(SURG) $(ALLI) $(LIBPATH)
	echo SURGTOOL:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(SURG) $(LIBS)

surgtool.debug: $(SURG) $(ALLI) $(LIBPATH_DEBUG)
	echo 'SURGTOOL (debug):'
	$(FC) $(FFLAGS) -o $(@) $(SURG) $(LIBS_DEBUG)

UAMC=uam2ncf.o

uam2ncf: $(UAMC) $(ALLI) $(LIBPATH)
	echo UAM2NCF:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(UAMC) $(LIBS)

uam2ncf.debug: $(UAMC) $(ALLI) $(LIBPATH_DEBUG)
	echo 'UAM2NCF (debug):'
	$(FC) $(FFLAGS) -o $(@) $(UAMC) $(LIBS_DEBUG)

#
## Module dependencies
#
include $(MD_SRC)/modmake.inc
include $(ETDSRC)/modmake.inc
include $(FS_ROOT)/modmake.inc

invsplit.o : $(MODINFO) 
	$(FC) $(FFLAGS) -c -o $@ invsplit.f

surgtool.o : $(MODGRID) $(MODSURG)
	$(FC) $(FFLAGS) -c -o $@ surgtool.f

smk2emis.o  : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ smk2emis.f

clean:
	/bin/rm -f *.o
