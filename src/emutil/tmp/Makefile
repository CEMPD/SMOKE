# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
##

all: mvcondns smk2emis uam2ncf

debug: mvcondns.debug smk2emis.debug uam2ncf.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

MVC=mvcondns.o

mvcondns: $(MVC) $(ALLI) $(LIBPATH)
	echo MVCONDNS:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(MVC) $(LIBS)

mvcondns.debug: $(MVC) $(ALLI) $(LIBPATH_DEBUG)
	echo 'MVCONDNS (debug):'
	$(FC) $(FFLAGS) -o $(@) $(MVC) $(LIBS_DEBUG)

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
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	/bin/rm *.o
