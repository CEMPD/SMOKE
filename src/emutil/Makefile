# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -ledsstools -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -ledsstools.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
##

all: beld3to2 geofac mvcondns mvsetup pktreduc smk2emis uam2ncf 

debug: beld3to2.debug geofac.debug mvcondns.debug mvsetup.debug pktreduc.debug smk2emis.debug uam2ncf.debug 

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	


ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

B3TO2=beld3to2.o

beld3to2: $(B3TO2) $(ALLI) $(LIBPATH)
	echo BELD3TO2:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(B3TO2) $(LIBS)

beld3to2.debug: $(B3TO2) $(ALLI) $(LIBPATH_DEBUG)
	echo 'BELD3TO2 (debug):'
	$(FC) $(FFLAGS) -o $(@) $(B3TO2) $(LIBS_DEBUG)


GEOF=geofac.o

geofac: $(GEOF) $(ALLI) $(LIBPATH)
	echo GEOFAC:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(GEOF) $(LIBS)

geofac.debug: $(GEOF) $(ALLI) $(LIBPATH_DEBUG)
	echo 'GEOFAC (debug):'
	$(FC) $(FFLAGS) -o $(@) $(GEOF) $(LIBS_DEBUG)

MVC=mvcondns.o

mvcondns: $(MVC) $(ALLI) $(LIBPATH)
	echo MVCONDNS:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MVC) $(LIBS)

mvcondns.debug: $(MVC) $(ALLI) $(LIBPATH_DEBUG)
	echo 'MVCONDNS (debug):'
	$(FC) $(FFLAGS) -o $(@) $(MVC) $(LIBS_DEBUG)

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

UAMC=uam2ncf.o

uam2ncf: $(UAMC) $(ALLI) $(LIBPATH)
	echo UAM2NCF:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(UAMC) $(LIBS)

uam2ncf.debug: $(UAMC) $(ALLI) $(LIBPATH_DEBUG)
	echo 'UAM2NCF (debug):'
	$(FC) $(FFLAGS) -o $(@) $(UAMC) $(LIBS_DEBUG)

MVS = mvsetup.o

mvsetup: $(MVS) $(ALLI) $(LIBPATH)
	echo MVSETUP:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MVS) $(LIBS)

mvsetup.debug: $(MVS) $(ALLI) $(LIBPATH_DEBUG)
	echo 'MVSETUP (debug):'
	$(FC) $(FFLAGS) -o $(@) $(MVS) $(LIBS_DEBUG)

#
#
# Modules
#
MODSOURC = $(MD_SRC)/modsourc.o
#
## Module dependencies
#
include $(MD_SRC)/modmake.inc

mvcondns.o  : $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ mvcondns.f

smk2emis.o  : $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ smk2emis.f

clean:
	/bin/rm *.o
