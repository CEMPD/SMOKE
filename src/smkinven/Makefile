# Version "@(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: smkinven.debug 
#
FC=f90
FFLAGS = $(DBGFLAGS)
LDFLAGS = -g
LDIRS = -L/env/proj/ioapi/IRIX6n32f90 \
        -L$(EDSS_LIB) \
        -L$(EL_SRC)/$(EDSS_OS)
LIBS=$(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH= $(EL_SRC)/$(EDSS_OS)/libsmoke.debug.a
#

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT

#
RAW=smkinven.o rdinven.o rdidapt.o rdemspt.o fixstk.o procinven.o openpnts.o \
     srcmem.o initinfo.o chklstfl.o wpntschr.o wpntspol.o wrptref.o \
     openinvin.o 
#
smkinven.debug: $(RAW) $(ALLI) $(LIBPATH)
	echo SMKINVEN:
	$(FC) $(FFLAGS) -o $(@) $(RAW) $(LIBS)

include $(MD_SRC)/modmake.inc

fixstk.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
initinfo.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openpnts.o: $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
procinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
rdemspt.o: $(MODSOURC) 
	$(FC) $(FFLAGS) -c -o $@ $< 
rdidapt.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
rdinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
smkinven.o: $(MODSOURC) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
srcmem.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
wpntschr.o: $(MODSOURC)
	$(FC) $(FFLAGS) -c -o $@ $<
#
clean:
	rm *.o *.s
