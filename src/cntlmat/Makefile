# Version "@(#)$Id$ $Source$ $Date$ 
#
all: debug 
#
debug: cntlmat.debug
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
CTL=cntlmat.o aloccmat.o alocpkts.o asgncntl.o errpkts.o \
    fillcdat.o fillcntl.o  genproj.o genreact.o openpmat.o \
    openrmat.o pktloop.o procpkts.o rdpacket.o wrrmat.o
#
cntlmat.debug: $(CTL) $(ALLI) $(LIBPATH)
	echo CNTLMAT:
	$(FC) $(FFLAGS) -o $(@) $(CTL) $(LIBS)
#
# Module dependencies
#
include $(MD_SRC)/modmake.inc

aloccmat.o: $(MODXREF) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
alocpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
asgncntl.o: $(MODSOURC) $(MODXREF) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
cntlmat.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
fillcdat.o: $(MODCNTRL)
	$(FC) $(FFLAGS) -c -o $@ $< 
fillcntl.o: $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genproj.o:  $(MODSOURC) $(MODCNTRL)
	$(FC) $(FFLAGS) -c -o $@ $< 
genreact.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO)
	$(FC) $(FFLAGS) -c -o $@ $< 
openpmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
openrmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
pktloop.o: $(MODXREF) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
procpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdpacket.o: $(MODCNTRL)
	$(FC) $(FFLAGS) -c -o $@ $< 

#
clean:
	rm *.o
