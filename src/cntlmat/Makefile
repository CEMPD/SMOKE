# Version @(#)$Id$ 
# Path    $Source$
# Date    $Date$ 
#
LIBS          = $(LDIRS) -lsmoke -lioapi -lnetcdf
LIBS_DEBUG    = $(LDIRS) -lsmoke.debug -lioapi.debug -lnetcdf.debug
LIBPATH       = $(SMK_BIN)/libsmoke.a
LIBPATH_DEBUG = $(SMK_BIN)/libsmoke.debug.a
#

all: cntlmat

debug: cntlmat.debug

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<	

ALLI=$(INC)/EMCNST3.EXT $(INC)/EMSTRG3.EXT 

#
CTL=cntlmat.o aloccmat.o alocpkts.o asgncntl.o errpkts.o \
    fillcdat.o fillcntl.o  genmultc.o genproj.o genreact.o \
    openpmat.o opencmat.o openctmp.o openrmat.o pktloop.o \
    procpkts.o rdpacket.o wrctmp.o wrrmat.o wcntlrep.o
#
cntlmat: $(CTL) $(ALLI) $(LIBPATH)
	echo CNTLMAT:
	$(FC) $(FFLAGS) -o $(SMK_BIN)/$(@) $(CTL) $(LIBS)
#
cntlmat.debug: $(CTL) $(ALLI) $(LIBPATH_DEBUG)
	echo CNTLMAT \(debug\):
	$(FC) $(FFLAGS) -o $(@) $(CTL) $(LIBS_DEBUG)
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
fillcdat.o: $(MODCNTRL) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ $< 
fillcntl.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
genmultc.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
genproj.o:  $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
genreact.o: $(MODSOURC) $(MODCNTRL) $(MODSPRO) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
opencmat.o: $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
openpmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
openrmat.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
pktloop.o: $(MODXREF) $(MODLISTS) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ $< 
procpkts.o: $(MODCNTRL) $(MODLISTS) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $< 
rdpacket.o: $(MODCNTRL) $(MODINFO) $(CL_SRC)/CPKTDAT.EXT
	$(FC) $(FFLAGS) -c -o $@ $<
wrctmp.o: $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<
wcntlrep.o: $(MODSOURC) $(MODCNTRL) $(MODINFO)
	$(FC) $(FFLAGS) -c -o $@ $<


#
clean:
	rm *.o
